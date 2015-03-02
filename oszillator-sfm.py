"""
Implementation of SFM to investigate oscillations. For dertails see
http://arxiv.org/abs/1412.1133
"""
#-------------------------------------------------------------
# from numpy import *

import numpy as np
import logging
import time
#----------------------- Parameter ---------------------------------------

fps = 16      # frames per second
dt = 0.001    # [s] integrator step length
t_end = 10    # [s] integration time
N_ped = 2     # number of pedestrians
Length = 10   # [m] length of corridor. *Closed boundary conditions*
RK4 = 0       # 1 --> RK4. 0 --> Euler. Euler with very small dt converges to RK4.


#======================================================
# init (x_n, dx_n) for n=0 --> N-1
def init(N):

    shift = float(Length)/N
    x_n = 0.5*shift + shift*np.arange(N)
    dx_n = np.zeros(N)
     
    Fd = np.zeros(N)
    Fr = np.zeros(N)
    
    return np.vstack( [x_n, dx_n] ), Fd, Fr    
#======================================================
# runge-kutta solver
def rk4(x, h, y, f):
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0
#======================================================
def euler(x, h, y, f):
    y_new, Fd, Fr = f(x, y)
    return x + h, y + h * y_new, Fd, Fr
#======================================================
def get_state_vars(state):
    """
    state variables and dist 
    """
    x_n = state[0,:]       # x_n  
    x_m = np.roll(x_n, -1)    # x_{n+1}
    x_m[ -1 ] += Length    # the first one goes ahead by Length
    dx_n = state[1,:]      # dx_n
    dx_m = np.roll(dx_n, -1)  # dx_{n+1}
    dist = x_m - x_n  
    return x_n, x_m, dx_n, dx_m, dist
#======================================================
def model(t, state):
    """
    exp-distance model nature 2000
    """
    x_n, x_m, dx_n, dx_m, dist = get_state_vars(state)

    #get acceleration => sum of forces
    f_drv = (v0 - dx_n)/tau
    f_rep = -a*np.exp(-dist/b)

    Fd = f_drv
    Fr = f_rep

    ########################
    a_n = f_rep + f_drv
    ########################
    x_n_new = dx_n
    dx_n_new = a_n

    x_n_new[-1] = 0   #ped 2 should not move
    dx_n_new[-1] = 0  #ped 2 should not move

    new_state = np.vstack([x_n_new, dx_n_new])
    return new_state, Fd, Fr
#======================================================
def simulation(N, t_end, state, Fd, Fr):
    t = 0
    frame = 0
    ids = np.arange(N_ped)
    while t <= t_end:
        if frame%(1/(dt*fps)) < 1e-8: # one frame per second

            T = t*np.ones(N_ped)
            output = np.vstack([ids, T, state, Fd, Fr]).transpose()
            np.savetxt(f, output, fmt="%d\t %f\t %f\t %f\t %f\t %f")
            f.flush()

        if RK4: # Runge-Kutta 4
            t, state = rk4(t, dt, state, model)
        else: # Euler
            t, state, Fd, Fr = euler(t, dt, state, model)

        frame += 1


#===================== MAIN =================================
if __name__=="__main__":
    b = 0.08 # 0.5 #0.08
    a = 2000 # 300 #2000
    tau = 0.5
    v0 = 1.2

    filename = "traj_N%d_a%.2f_b%.2f_v0%.2f_tau%.2f.txt"%(N_ped, a, b, v0, tau) #trajectory file

    f = open(filename, 'w+')

    state, Fd, Fr = init(N_ped)

    t1 = time.clock()

    print 'filename ', filename
    constant = float(a)/b*(tau**2)
    if constant < 0.25:
        print 'No oscillations: constant = %.3f'%constant
    else:
        print 'Oscillations: constant = %.3f'%constant
     ######################################################
    simulation(N_ped, t_end, state, Fd, Fr)
    ######################################################
    t2 = time.clock()
    f.close()
