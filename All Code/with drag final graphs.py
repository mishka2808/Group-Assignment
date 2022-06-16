
"""
Created on Wed Jun 15 20:00:52 2022

@author: Mish
"""
import numpy as np
import matplotlib.pyplot as plt

Ax = -20    #Point A where the mass starts
Ay = 0
Bx = 20     #Point B for other attached point
By = 0
g  = 9.8    #gravity
L  = 60     #length of bungee cord
m  = 90     #mass
mu = 0.77   #air friction
k  = 10     #spring constant for bungee cord
N  = 480    #sim steps
dt = 0.25   #time step
t  = 0      #initial time
x  = Ax     #initial x pos of mass
y  = Ay     #initial y pos of mass
vx = 0      #initial x velocity
vy = 0      #initial y velocity

def dist(x1,y1,x2,y2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)
 
def normalize_vec(v):
    v_mag = np.sqrt(v[0]*v[0] + v[1]*v[1])
    if v_mag>0:
        return (v[0]/v_mag,v[1]/v_mag)
    return (0,0)
    
def length_vec(v):
    return np.sqrt(v[0]*v[0] + v[1]*v[1])
    
def unit_vec(x1,y1,x2,y2):  #unit vector from (x1,y1) to (x2,y2)
    dx = x2-x1
    dy = y2-y1
    return normalize_vec((dx,dy))
    
def F_A():              #F_A is spring force in bungee cord towards point A
    d = dist(x,y,Ax,Ay)
    if d<= L:
        return (0,0)    #No force if bungee cord is slack
    else:
        dL = d-L                        #extension length
        F_mag = k*dL                    #magnitude of spring force
        u = unit_vec(x,y,Ax,Ay)         #unit vector in directtion of force
        F_vec = (F_mag*u[0],F_mag*u[1]) #spring force vector
        return F_vec

def F_B():              #F_B is spring force in bungee cord towards point B
    d = dist(x,y,Bx,By)
    if d<= L:
        return (0,0) #No force if bungee cord is slack
    else:
        dL = d-L                        #extension length
        F_mag = k*dL                    #magnitude of spring force
        u = unit_vec(x,y,Bx,By)         #unit vector in directtion of force
        F_vec = (F_mag*u[0],F_mag*u[1]) #spring force vector   
        return F_vec

def F_D():           #F_D is drag force in opposite direction to motion
    v_hat = normalize_vec((vx,vy))
    v = length_vec((vx,vy))
    F_mag = mu*v*v
    F_vec = (-F_mag*v_hat[0],-F_mag*v_hat[1])
    return F_vec
   
def draw_vec(p,v,plt):
    plt.plot([p[0],p[0]+v[0]],[p[1],p[1]+v[1]],color="Blue")
    
#F_net = F_g + F_A + F_B + F_D

#F_net is net force acting on mass m at time t
#F_g is force of gravity going down

#Runge Kutta integration functions
#dx/dt = f1 = vx
#d(vx)/dt = g1 = F_net_x

#dy/dt = f2 = vy
#d(vy)/dt = g2 = F_net_y

def f1(t,x,vx):
    return vx

def g1(t,x,vx):
    return (F_A()[0] + F_B()[0] + F_D()[0])/m

def f2(t,y,vy):
    return vy
    
def g2(t,y,vy):
    return (F_A()[1] + F_B()[1] + F_D()[1])/m - g
    
x_values = []
y_values = []
E_kin = []
E_pot = []
E_tot = [] 
stepper = []

for step in np.arange(N):
    k0 = dt*f1(t,x,vx)
    l0 = dt*g1(t,x,vx)
    
    k1 = dt*f1(t+0.5*dt,x+0.5*k0,vx+0.5*l0)
    l1 = dt*g1(t+0.5*dt,x+0.5*k0,vx+0.5*l0)
    
    k2 = dt*f1(t+0.5*dt,x+0.5*k1,vx+0.5*l1)
    l2 = dt*g1(t+0.5*dt,x+0.5*k1,vx+0.5*l1)
    
    k3 = dt*f1(t+dt,x+k2,vx+l2)
    l3 = dt*g1(t+dt,x+k2,vx+l2)
    
    x += (1/6)*(k0+2*k1+2*k2+k3)
    vx += (1/6)*(l0+2*l1+2*l2+l3)
    
    k0 = dt*f2(t,y,vy)
    l0 = dt*g2(t,y,vy)
    
    k1 = dt*f2(t+0.5*dt,y+0.5*k0,vy+0.5*l0)
    l1 = dt*g2(t+0.5*dt,y+0.5*k0,vy+0.5*l0)
    
    k2 = dt*f2(t+0.5*dt,y+0.5*k1,vy+0.5*l1)
    l2 = dt*g2(t+0.5*dt,y+0.5*k1,vy+0.5*l1)
    
    k3 = dt*f2(t+dt,y+k2,vy+l2)
    l3 = dt*g2(t+dt,y+k2,vy+l2)
    
    y += (1/6)*(k0+2*k1+2*k2+k3)
    vy += (1/6)*(l0+2*l1+2*l2+l3)

    x_values.append(x)
    y_values.append(y)
    
    E_pot.append(m*g*y)
    E_kin.append(0.5*m*np.sqrt(vx**2 + vy**2) )
    E_tot.append(m*g*y + 0.5*m*np.sqrt(vx**2 + vy**2) )
    
    stepper.append(t)
    
    t += dt
    
plt.plot(x_values, y_values, '.-k') 
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("Trajectory")
plt.legend()  
plt.show()

plt.plot(stepper, E_pot, '--b', label = "Potential Energy")
plt.plot(stepper, E_kin, '--r', label = "Kinetic Energy")
plt.plot(stepper, E_tot, '-k', label = "Equilibruim Energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Kinetic, Potential and Equilibruim Energy")
plt.legend()  
plt.show()