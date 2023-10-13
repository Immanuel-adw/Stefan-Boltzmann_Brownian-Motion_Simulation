 # -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 12:17:46 2019

@author: ioa15
"""

import numpy as np
import scipy as sp
import pylab as pl
import matplotlib.animation as anim
import math as math
import random as random
from scipy import optimize


class Ball:
    """ 
    a ball
    """
    def __init__(self, m =1, R = 1, r= np.array([0,0]), v = np.array([0,0]), iscont = False):
        #self.patches = []
        self.m = m
        self.R = R
        self._r = np.array(r)
        self.v = np.array(v)
        self.iscont = iscont
        if iscont == False:
            self.patch = pl.Circle(self.get_r(), R, fc='r')
        else:
            self.patch = pl.Circle([0., 0.], R, ec='b', fill=False, ls='solid')
    
    def get_r(self):
        return np.array(self._r)
           
    def get_patch(self):
        return self.patch
    
    def set_r(self, new_r):
        self._r = np.array(new_r)
    
    #np.finfo(np.float32).eps
    def time_to_collision(self, other):
        r_sep = self.get_r() - other.get_r()
        v_sep =  self.v - other.v
        if  self.iscont ==True or other.iscont == True:
            R_vec = self.R - other.R
        else:
            R_vec = self.R + other.R
        a = np.dot(v_sep,v_sep)
        b = 2*np.dot(r_sep,v_sep)
        c = np.dot(r_sep, r_sep) - (R_vec*R_vec)
        
        disc = b**2-(4*a*c)
                #print(disc)
        if disc < 0:
            return np.nan
        elif disc == 0:    # one real solution
            dt = (-b + np.sqrt(disc))/(2*a)
            return dt
        else: 
            dt1 = (-b - np.sqrt(disc)) /(2*a)
            dt2 = (-b + np.sqrt(disc)) /(2*a)
            solutions = [dt1, dt2]
            if dt1 < np.finfo(np.float32).eps and dt2 < np.finfo(np.float32).eps: # if there are two positive solutions return the smaller of the two
                return np.nan 
            else: # if one solution is ngative take the larger one and check whether it is positive
                 return min(i for i in solutions if i > np.finfo(np.float32).eps)
        

        
    def move(self,dt):
        self.set_r(self.get_r() + self.v*dt)
        self.patch.center = self.get_r()
    
    def vel(self):
        return self.v
    
    def pos(self):
        return self.get_r()
    
    
        
    def collide(self, other): #cpresent means container collision
        
        r_sep = self.get_r() - other.get_r()             
        v1 = self.v
        v2 = other.v
        m1 = self.m
        m2 = other.m
        m_sum = m1 + m2
        m_diff = m1 - m2
        r_unit_para = (r_sep)/(abs(np.dot(r_sep, r_sep))**0.5)
        r_unit_perp = np.array([-r_unit_para[1], r_unit_para[0]])
        v1_para = np.dot(v1, r_unit_para)
        v1_perp = np.dot(v1, r_unit_perp)
        v2_para = np.dot(v2, r_unit_para)
        v2_perp = np.dot(v2, r_unit_perp)
        v1_para_new = (m_diff/m_sum) * v1_para + (2*m2/m_sum) * v2_para
        v2_para_new = (-m_diff/m_sum) * v2_para + (2*m1/m_sum) * v1_para
        v1b = v1_para_new * r_unit_para + v1_perp * r_unit_perp
        v2b = v2_para_new * r_unit_para + v2_perp * r_unit_perp
        self.v = v1b
        other.v = v2b
       



class simulation:
       
    def __init__(self, numberballs, cont_R, ball_R, r_max, m, v_max, timeelapsed = 0): 
        self.numberballs = numberballs 
        self.cont = Ball(m=float(10e25), R=cont_R, iscont=True) 
        #self.ball_list = setup(ball_list, 10, .2, 20)
        self.cont_R = cont_R
        self.ball_R = ball_R
        self.ball_m = m
        self.generate_random_r()
        self.v_max = v_max
        self.r_max = r_max
        self.ball_list = []
        self.tot_del_mom = 0
        self.tot_mom = 0
        self.containerhitcount = 0
        self.t_elapsed = 0
        
        self.timelist = []
  
        self.ball_list = []
        
        self.pressures = []
        
        self.vdistx = [] 
        
        self.vdisty = []
        
        self.speedlist = []
        
        self.r_list = []
        
        self.templist = []
        
        #ball = Ball(self.ball_m, self.ball_R, self.r_max, self.v_max)
        for i in range(numberballs):
            self.create_ball(0)
        self.ball_list.append(self.cont)
        self.pressures = []
        
        
    def ballmass(self):
        return self.ball_m
    
    def balls_listed(self):
        return self.ball_list
        
    def create_ball(self, gas_type = 0):
        """
        Creates a 'ball' in 'gas_system' object. This method is called during instantiation and should not be called
        manually by the user.
        """
        #print(self.__possible_r_remaining)
        r = random.choice(self.__possible_r_remaining)
        v = np.array([self.v_max*np.random.uniform(-1.0,1.0),self.v_max*np.random.uniform(-1.0,1.0)])
        #print('v is', v)
        self.__possible_r_remaining.remove(r)       
        b = Ball(self.ball_m, self.ball_R, r, v, iscont = False)       
        self.ball_list.append(b)
            
    def generate_random_r(self):
        """
        Generates a list of positions which a gas particle can be instantiated with. This ensures that particles
        will not overlap.
        """
        def generate(R, T):
            for j in range(len(R)):
                #print('len(R) is', len(R))
                #print('R[j] is', R[j] )
                #print('rad_calc is', rad_calc)
                
                for i in range(T[j]):
                    #print('T[j] is', T[j] )
                    r = R[j]
                    t = i*2*np.pi/T[j]
                    #print('len(T[j]) is', len(T[j]))
                    #print('t is', t)
                    yield r,t                
            
        count = 0
        R = []
        while count < self.cont_R - self.ball_R:
            R.append(count)
            count += 2.5*self.ball_R
            
        T=[]
        for i in range(len(R)):
            T.append(int(math.trunc(R[i]*2*np.pi/(self.ball_R*2.00000000001))))
            #print('T is', T)
        T[0] = 1
                
        self.__possible_r_remaining = []
        for r, t in generate(R, T):
            self.__possible_r_remaining.append([r*np.cos(t), r*np.sin(t)])
       

    """"
    def init_anim(self):
        patches = []
        for i in range(self.numberballs):
            b = self.ball_list[i]
            patch = b.get_patch()
            addpatch = ax.add_patch(patch)
            patches.append(addpatch)
            """
    def ball_speeds(self):
        """
        A list containing the speeds of all balls.
        """
        self.v_all_balls = []
        
        for i in self.ball_list:
            self.v_all_balls.append(np.linalg.norm((i.ball__v)))
        
        return self.v_all_balls
    
    def total_kinetic_energy(self):
        totKE = 0
        for a in range(0, self.numberballs):
            temp_ball = self.ball_list[a]
            totKE += 0.5* self.ball_m * np.dot(temp_ball.v, temp_ball.v)
        #print(totKE)
        return totKE
    
    def av_KE(self):
        self.avKE = self.total_kinetic_energy()/ self.numberballs
        return self.avKE
    
    def Temperature(self):
        kb = 1.38 * 10e-23
        Temp = (self.av_KE())/kb
        self.templist.append(Temp)
        #print('temp is', Temp)
        return Temp
    
    
    def check_KE(self):
        EKlist = []
        for i in range(500):
             #print(i)
             self.run(10)
             EKlist.append(self.av_KE())
        #print('EKlist', EKlist)
        return EKlist
    
    
    def next_collision(self):
        conttimelist = []
        dtlist = []
        for i in self.ball_list:
            #dtjlist = []
            #KEb = 0.5 * m * (v**2)
            #print('KE before is', KEb)
            for j in self.ball_list:
                if self.ball_list.index(j) <= self.ball_list.index(i):
                    dtlist.append(1000000)
                else:
                    dt = i.time_to_collision(j)
                    if j is self.cont:
                        conttimelist.append(dt)
                        dtlist.append(dt)
                    #dtjlist.append(dt)
                    #dt = np.nanmin(dtjlist)
                    else:
                        dtlist.append(dt)
                    
        dt = np.nanmin(dtlist)
        
        for i in self.ball_list:
            #print(dt)
            i.move(dt)
            
        
        #if animate:
            #for i in range(self.numberballs):
                #bi = self.ball_list[i]
                #patches.append(bi.get_patch())
        if dt in conttimelist:
            dt_index = np.nanargmin(dtlist)
            index = dt_index // (len(self.ball_list))
            
            vbefore = self.ball_list[index].v 
            
            self.ball_list[index].collide(self.cont)
            vafter = self.ball_list[index].v
            
            del_v = vafter - vbefore
            
            delmomentum = self.ball_list[index].m * abs(del_v)
            self.tot_del_mom += delmomentum
            
            if self.containerhitcount < 1:
                self.containerhitcount += 1
                self.t_elapsed += dt                
            else:                    
                self.containerhitcount += 1
                self.t_elapsed += dt
                self.timelist.append(self.t_elapsed)  
                #self.delmomlist.append(delmomentum)
                Force = ((self.tot_del_mom)/self.t_elapsed)
                ForceMag = np.sqrt(np.dot(Force,Force))
                #print(Force)
                pressure = ForceMag /(2*np.pi*self.cont_R)
                #print(pressure)
                self.pressures.append(pressure)
                self.containerhitcount = 0
                self.tot_del_mom = 0
                self.t_elapsed = 0
            
           
        else:
            dt_index = np.nanargmin(dtlist)
            ball1_index = dt_index//(len(self.ball_list))
            ball2_index = dt_index%(len(self.ball_list))
            b1 = self.ball_list[ball1_index]
            b2 = self.ball_list[ball2_index]
            #print("colliding", b1,"and",b2)
            b1.collide(b2)
            self.t_elapsed += dt
            
        return dt
        

        
        self.delmomlist = []
    
    
    
    def find_average_pressure(self):
        ps = self.pressures[len(self.pressures)-3:-1]
        av_p = sum(ps)/len(ps)
        return av_p
            
        """ change next collision back to normal"""
        

    
    #vdist = np.sqrt((np.array(vdistx)**2)+np.array(vdisty)**2)
    
             
    def run(self, num_frames, animate=False):
        x = self.cont.R
        if animate:
            pl.figure(figsize=(5,5))
            ax = pl.axes(xlim=(-x, x), ylim=(-x, x))
            ax.add_artist(self.cont.get_patch())
        if animate:
            patches = []
            for i in self.ball_list:
                patch = i.get_patch()
                addpatch = ax.add_patch(patch)
                patches.append(addpatch)
                #ax.add_artist(pl.Circle([0,0], self.ball_list[-1].R, ec='b', fill=False, ls='solid'))
                #self.init_anim()                
                #pl.pause(0.0001)
        for frame in range(num_frames):
            #print(frame)
            for i in self.ball_list:
                self.vdistx.append(i.v[0])                
                self.vdisty.append(i.v[1])
                if self.ball_list.index(i) < len(self.ball_list) - 1: 
                    self.speedlist.append(np.sqrt((i.v[0]**2)+(i.v[1]**2)))
                    self.r_list.append(np.sqrt((i.get_r()[0]**2)+(i.get_r()[1]**2)))
                else:
                    pass
            self.next_collision()
            self.total_kinetic_energy()
            if animate:
                pl.pause(0.0001)
        if animate:
            pl.figure(1)
            pl.show()
        #print(self.pressures)
        #print(timelist
        
    def plot(self):
        mass = self.ball_m
        avKE = self.av_KE()
        kT = avKE
        A = np.sqrt((mass)/(2*np.pi*kT))
        C = 2*np.pi*(A*A) 
        alpha = (mass)/(2*kT)
        
        
        pl.figure(2)
        pl.suptitle('Pressure versus time')
        pl.plot(self.timelist, self.pressures, 'bx')
        pl.xlabel('Time (s)')
        pl.ylabel('Pressure (Pa)')
        pl.show()
        
        
        pl.figure(5)
        z = np.linspace(0, 18, 3000)
        pl.suptitle('Maxwell Boltzmann distribution')
        pl.ylabel('Probability density')
        pl.xlabel('velocity2d')
        pl.hist(self.speedlist, bins=75, normed=1)
        pl.plot(z, C * z * np.exp(-alpha*(z)**2))
        pl.show()
        
        pl.figure(6)
        z = np.linspace(-10, 10, 3000)
        pl.suptitle('pd of Particle Separation from Center')
        pl.ylabel('Probability density')
        pl.xlabel('separation from Centre')
        pl.hist(self.r_list, bins=75, normed=1)
        #pl.plot(z, C * np.exp(-alpha*(z)**2))
        pl.show()
        
        pl.figure(7)
        pl.suptitle('Ek conservation checker')
        pl.ylabel('EK conservation checker')
        pl.xlabel('count')
        pl.plot(np.arange(0,5000,10), self.check_KE(), 'bo')
        pl.show()   

l = 100
n = 100
v = 1
mb = .01
ballr = 0.05

sim = simulation(numberballs= n,  cont_R=10, ball_R= ballr, r_max=6, m= mb, v_max=v)
sim.run(l)
sim.plot()
 


sim_pV = []
average_p_for_pV = []
Cont_vol = []
TemperatureV = []

sim_pT = []
average_p_for_pT = []
Temperature1 = []

for i in range(0, 9):
    sim_pV.append(simulation(numberballs= n,  cont_R = (i+1) * 10, ball_R=ballr, r_max=6, m= mb, v_max=v))


  
for i in range(0, 9):
    spV = sim_pV[i]
    spV.run(l)
    Avogadros = 6.02 * 10e23
    b1 = np.pi * ((spV.ball_R)**2) * (spV.numberballs)
    Cont_vol.append(1/((np.pi * ((spV.cont_R)**2))-b1))
    #print('spvballr is', spV.ball_R)
    #print('spvnumberballs is', spV.numberballs)
    #print((np.pi * ((spV.cont_R)**2))-b1)
    average_p_for_pV.append(spV.find_average_pressure())
    TemperatureV.append(spV.Temperature())
    #print('Temp V is', TemperatureV)
    
    
    
for i in range(10, 19):
    sim_pT.append(simulation(numberballs= n,  cont_R = 10, ball_R = ballr, r_max=6, m= mb, v_max=i *v))
   
for i in range(0, 9):
    spT = sim_pT[i]
    spT.run(l)
    Temperature1.append(spT.Temperature())
    average_p_for_pT.append(spT.find_average_pressure())
    #print(Temperature1)
       
#print('avtemp is ' + str(avtemp))
    
    
Cont_vol = np.array(Cont_vol)
average_p_for_pV = np.array(average_p_for_pV)
Temperature1 = np.array(Temperature1)
average_p_for_pT = np.array(average_p_for_pT)



def lineS(x,a,b):
    linear = (a * x) + b
    return linear

popt1, pcov1 = optimize.curve_fit(lineS, Cont_vol, average_p_for_pV)
a = popt1[0]
b = popt1[1]


def line2(x, c, d):
    linear2 = (c * x) + d
    return linear2

popt2, pcov2 = optimize.curve_fit(line2, Temperature1, average_p_for_pT)
c = popt2[0]
d = popt2[1]

avtemp = sum(np.array(TemperatureV))/len(TemperatureV)
print('avtemp is', avtemp)
s = sim_pV[-1]
kb = 1.38 * 10e-23
grad1 =   (spV.numberballs) * kb * avtemp
Err1 = np.sqrt(pcov1[0,0])


pl.figure(8)
pl.suptitle('Pressure vs 1/(V-b)' )
pl.xlabel('1/(V-b)') 
pl.ylabel('Pressure')   
pl.plot(Cont_vol, average_p_for_pV, 'b+', label= ('temperature = ' + str(avtemp), 'grad = ' + str(grad1)))
pl.plot(Cont_vol, lineS(Cont_vol, a, b))
pl.legend(loc='upper right')
pl.show()

print('grad1 is '+ str(a)+' +/- '+ str(Err1))


Avogadros = 6.02 * 10e23
V = np.pi * (spV.cont_R**2)
b1 = np.pi * (spV.ball_R**2) * spV.numberballs

s = sim_pT[-1]
gradcalc2 =   (spV.numberballs * kb)/(V-b1)
Err2 =  np.sqrt(pcov2[0,0])
 
  
pl.figure(9)
pl.suptitle('Pressure vs Temperature' )
pl.xlabel('Temperature') 
pl.ylabel('Pressure')   
pl.plot(Temperature1, average_p_for_pT,  'b+', label= 'grad = ' + str(gradcalc2))
pl.plot(Temperature1, line2(Temperature1, c, d))
pl.legend(loc='upper right')
pl.show()              
print('grad2 is '+ str(c)+' +/- '+ str(Err2))
        

            
            
