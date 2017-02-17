# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:21:46 2017

@author: Сергей
"""

print("dasda")

import numpy 
from matplotlib import pyplot
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
#############################################################
#Basic initial condition parameters
#defining grid size, time steps, CFL condition, etc...
nx = 681
T = 0.01
dx = 20./(nx-1)
sigma1 = 0.0008
dt = sigma1*dx
nt = T/dt
gamma = 1.4

#############################################################
x = numpy.linspace(-10,10,nx)
#############################################################

def u_initial(nx):
    
    Rho = numpy.ones(nx)
    Rho[int((nx-1)/2.):] = 0.125
    v = numpy.zeros(nx)
    p = 100000*numpy.ones(nx)
    p[int((nx-1)/2.):] = 10000
    return numpy.array([Rho, v, p])
    
#############################################################    
z = u_initial(nx)        #initial conditions
print(z)
#############################################################
def f(u, gamma):
    
    """
    Functia potokov konservativnih peremenih
    """
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    
    
    return numpy.array([u2, u2*u2/u1+(gamma-1)*(u3-0.5*u2*u2/u1),
                      (u3+(gamma-1)*(u3-0.5*u2*u2/u1))*u2/u1])
                      
 #############################################################                     
                      
def f1(u, gamma):
    "Functia perehoda ot primitivnih k konservativnim peremennim"
    Rho = u[0]
    v = u[1]
    p = u[2]
    e = p/((gamma-1)*Rho)
    eT = e+0.5*v*v
    return numpy.array([Rho, Rho*v, Rho*eT]) 

#############################################################

u = f1(z, gamma) #flux variables 


u_star = numpy.zeros((len(u),len(u[0])))
u_n = numpy.zeros_like(u)      
    #copy the initial u array into each row of our new array
print(u.shape,u_star.shape)


#############################################################

def Godunov(u,nt,dt,dx,gamma):
    u_plus = numpy.zeros_like(u)
    u_minus = numpy.zeros_like(u)
    flux = numpy.zeros_like(u)
    for n in range(0,int(nt)):
    #print('next step')
    
        u_n = u.copy() 
    #print(n)
    
    
   
    
    
    
        u_plus[:,:-1] = u[:,1:] # Can't do i+1/2 indices, so cell boundary
        u_minus = u.copy() # arrays at index i are at location i+1/2
        flux = 0.5 * (f(u_minus, gamma) + 
                      f(u_plus, gamma) + 
                      dx / dt * (u_minus - u_plus))
        u_n[:,1:-1] = u[:,1:-1] + dt/dx*(flux[:,:-2]-flux[:,1:-1])
        u_n[:,0] = u[:,0]
        u_n[:,-1] = u[:,-1]
        u = u_n.copy()
    return u

#############################################################

def f3(u, gamma):
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
   
    return numpy.array([u1, u2/u1, (u3/u1-0.5*u2*u2/u1/u1)*((gamma-1)*u1)])

#############################################################  

 
"polucenie resenia"
U2 = Godunov(u,nt,dt,dx,gamma)

z3 = f3(U2, gamma)
"zapolnenie massiva"

p3 = numpy.arange(len(u[0]))
p3 = z3[2,:]

rho3 = numpy.arange(len(u[0]))
rho3 = z3[0,:]
v3 = numpy.arange(len(u[0]))
v3 = z3[1,:]



#############################################################  

def func1(P,p,rho):
    gamma = 1.4 
    s = (gamma*p/rho)**(1/2)
    for P in P:
        if p > P:
            return 2*s/(gamma-1)*((P/p)**((gamma-1)/(2*gamma))-1)
        else:
            return (P-p)/(rho*s)*((gamma+1)/(2*gamma)*P/p+(gamma-1)/(2*gamma))**(-1/2)

############################################################# 
"""
def func2(P,p,rho):
    gamma = 1.4 
    s = (gamma*p/rho)**(1/2)
    
    return (P-p)/(rho*s)*((gamma+1)/(2*gamma)*P/p+(gamma-1)/(2*gamma))**(-1/2)
"""
############################################################# 

def func1_shtrih(P,p,rho):
    gamma = 1.4 
    s = (gamma*p/rho)**(1/2)
    for P in P:
        if p > P:
            return s/(gamma*P)*(P/p)**((gamma-1)/(2*gamma))
        else:
            return ((gamma+1)*P/p+(3*gamma-1))/(4*gamma*p*s)*((gamma+1)/(2*gamma)*P/p+(gamma-1)/(2*gamma))**(-3/2)

############################################################# 
"""
def func2_shtrih(P,p,rho):
    gamma = 1.4 
    s = (gamma*p/rho)**(1/2)
    
    return ((gamma+1)*P/p+(3*gamma-1))/(4*gamma*p*s)*((gamma+1)/(2*gamma)*P/p+(gamma-1)/(2*gamma))**(-3/2)
"""
############################################################# 

"vichislenie tocinogo reshenia" 
P_exact1 = numpy.zeros(len(u[0]))

P_exact1[:] = p3[:] - (func1(p3[:],100000,1)+func1(p3[:],10000,0.125))/(func1_shtrih(p3[:],100000,1)+func1_shtrih(p3[:],10000,0.125))   


P_exact2 = numpy.zeros(len(u[0]))
    
P_exact2[:] = P_exact1[:] - (func1(P_exact1[:],100000,1)+func1(P_exact1[:],10000,0.125))/(func1_shtrih(P_exact1[:],100000,1)+func1_shtrih(P_exact1[:],10000,0.125))

print(P_exact1)  
print(P_exact2)
pyplot.plot(x, P_exact1, color='green', ls='-', lw=3)
#pyplot.plot(x, p3, color='red', ls='-', lw=3)
pyplot.ylabel('pressure')
pyplot.xlabel('Distance')



























