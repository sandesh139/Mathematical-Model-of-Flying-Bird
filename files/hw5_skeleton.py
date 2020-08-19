from scipy import *
from matplotlib import pyplot
import matplotlib.animation as manimation
import os, sys

# HW5 Skeleton 

def RK4(f, y, t, dt, food_flag, alpha, gamma_1, gamma_2, kappa, rho, delta):
    '''
    Carry out a step of RK4 from time t using dt
    
    Input
    -----
    f:  Right hand side value function (this is the RHS function)
    y:  state vector
    t:  current time
    dt: time step size
    
    food_flag:  0 or 1, depending on the food location
    alpha:      Parameter from homework PDF
    gamma_1:    Parameter from homework PDF
    gamma_2:    Parameter from homework PDF
    kappa:      Parameter from homework PDF
    rho:        Parameter from homework PDF
    delta:      Parameter from homework PDF
    

    Output
    ------
    Return updated vector y, according to RK4 formula
    '''
    
    # Task: Fill in the RK4 formula
    k1 = f(y,t,food_flag,alpha,gamma_1,gamma_2,kappa,rho,delta)
    k2 = f(y+(dt/2)*k1,t+(dt/2),food_flag,alpha,gamma_1,gamma_2,kappa,rho,delta)
    k3 = f(y+(dt/2)*k2,t+(dt/2),food_flag,alpha,gamma_1,gamma_2,kappa,rho,delta)
    k4 = f(y+dt*k3,t+dt,food_flag,alpha,gamma_1,gamma_2,kappa,rho,delta) 

    return y + (dt/6)*(k1+2*k2+2*k3+k4)


def RHS(y, t, food_flag, alpha, gamma_1, gamma_2, kappa, rho, delta):
    '''
    Define the right hand side of the ODE

    '''
    N = y.shape[0] #The number of bird in the flock
    f = zeros_like(y)
    # Task:  Fill this in by assigning values to f
    
    # Set f_food
    f_food = zeros_like(y)
    if food_flag == 1:
        C = array([sin(alpha*t),cos(alpha*t)],dtype='float')
    elif food_flag == 0:
        C = array([0,0],dtype='float')
    f_food[0,] = gamma_1*(C-y[0,])

    # Set f_follow
    f_follow = zeros_like(y)
    for i in range(1,N):
        f_follow[i,] = y[0,]
    f_follow[1:,] = gamma_2*(f_follow[1:,]-y[1:,])

    #set_f_flock
    f_flock = zeros_like(y)
    B_bar = flockCenter(y)
    f_flock[1:,] = kappa*(B_bar[1:,]-y[1:,])

    # Set f_rep
    f_rep = zeros_like(y)
    for j in range(2):
        for k in range(1,N):
            neighbors = []
            findNeighbors(y,k,neighbors)
            summ = 0
            for x in neighbors:
                summ += rho*(y[k,j]-y[x,j])/((y[k,j]-y[x,j])**2 + delta)
            f_rep[k,j] = summ

    return f_flock + f_food + f_follow + f_rep

#Find the distance between two birds
def distanceFun(first, second):
    return math.sqrt((first[0]-second[0])**2+(first[1]-second[1])**2)

#find the 5 closest neighbors
def findNeighbors(y,k,neighbors):
    distances = zeros(len(y))
    closest = 100

    for i in range(len(y)):
        distances[i] = distanceFun(y[k],y[i])
        if distances[i] < closest and i != k:
            closest = distances[i]
            closestIndex = i
    
    distances[closestIndex] = 100
    neighbors.append(closestIndex)

    while len(neighbors) < 5:
        closest = 100
        for j,d in enumerate(distances):
            if d < closest and j != k:
                closest = d
                closestIndex = j
        neighbors.append(closestIndex)
        distances[closestIndex] = 100

# This function returns an N by 2 matrix whos entries are all the coordinants of the center of the flock
def flockCenter(y):
    B_bar = ones((N,2))
    B_bar[:,0] = B_bar[:,0]*(sum(y[:,0])/N)
    B_bar[:,1] = B_bar[:,1]*(sum(y[:,1])/N)
    return B_bar


# We'll define flock diamter as the doubled distance from the flock center to the furthest bird
def flockDia(y):
    maxDist = 0
    center = flockCenter(y)
    for i in range(1,len(y)-1):
        dist = distanceFun(center[i],y[i])
        if dist > maxDist:
            maxDist = dist

    return 2*maxDist

##
# Set up problem domain
t0 = 0.0        # start time
T = 10.0        # end time
nsteps = 50     # number of time steps

# Task:  Experiment with N, number of birds
N = 10

# Task:  Experiment with the problem parameters, and understand how each parameter affects the system
dt = (T - t0) / (nsteps-1.0)
gamma_1 = 2.0
gamma_2 = 8.0
alpha = 0.4
kappa = 4.0
rho = 2.0
delta = 0.5
food_flag = 1   # food_flag == 0: C(x,y) = (0.0, 0.0)
                # food_flag == 1: C(x,y) = (sin(alpha*t), cos(alpha*t))

# Intialize problem
y = rand(N,2)  # This is the state vector of each Bird's position.  The k-th bird's position is (y[k,0], y[k,1])
flock_diam = zeros((nsteps,))


# Initialize the Movie Writer
# --> The movie writing code has been done for you
FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=6)

fig = pyplot.figure(0)
pp, = pyplot.plot([],[], 'k+') 
rr, = pyplot.plot([],[], 'r+')
pyplot.xlabel(r'$X$', fontsize='large')
pyplot.ylabel(r'$Y$', fontsize='large')
pyplot.xlim(-3,3)       # you may need to adjust this, if your birds fly outside of this box!
pyplot.ylim(-3,3)       # you may need to adjust this, if your birds fly outside of this box!


# Begin writing movie frames
with writer.saving(fig, "movie.mp4", dpi=1000):

    # First frame
    pp.set_data(y[1:,0], y[1:,1]) 
    rr.set_data(y[0,0], y[0,1])
    writer.grab_frame()

    t = t0
    for step in range(nsteps):
        
        # Task: Fill in the code for the next two lines 
        y = RK4(RHS,y,t,dt,food_flag,alpha,gamma_1,gamma_2,kappa,rho,delta) 
        flock_diam[step] = flockDia(y)
        t += dt
        
        # Movie frame
        pp.set_data(y[:,0], y[:,1]) 
        rr.set_data(y[0,0], y[0,1])
        writer.grab_frame()
 
# Task: Plot flock diameter
pyplot.figure(1)
pyplot.plot(linspace(0,10,nsteps),flock_diam)
pyplot.xlabel('Time', fontsize='large')
pyplot.ylabel('Flock Diameter', fontsize='large')
pyplot.title('Flock Diameter over time')
pyplot.savefig('FlockDiam'+str(N)+'v3.png',dpi=500, format='png', bbox_inches='tight', pad_inches=0.0)
pyplot.show(1)
