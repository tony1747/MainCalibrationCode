import numpy as np
import matplotlib.pyplot as plt
from random import random
from math import log
from copy import copy
import seaborn as sns
from scipy.integrate import solve_ivp

import pandas


class Error(Exception):
    """Base class for other exceptions"""
    pass


class NegativeRate(Error):
    """Negative Birth or Death Rate"""
    pass


data = []

replicates = 1
cutoff = 2


# varFactor is birth + death for all cell types in all drugs:
# try for K=1, 2, 3
d = 0.0598
a = 0.3953
b = 1.33e-5
p = 0.4862
m = 1.784e-6
k_n1 = 0
k_2 = 1
k_3 = 8.8577e-4*k_2  #k_3/k_2=8.8577e-4
g = 1.7531e5
s = 0
k_1 = m*(k_2+k_3+k_n1)/k_3 # try a K that makes sense
d_1 = d

# try for T_initial=10^7 or 8, E_initial = 2*10^7 or any reasonable number
relativePopInit = [2500, 2500*5];

# Check for negative rates:
for parameter in [a, b, d, p, k_1, k_n1, k_2, k_3]:
    if parameter < 0:
        raise NegativeRate






# T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1



# -----------------Stoch----------------------------------------
def first_n_sum(p,n):  # define a function that can sum up the first n terms of a list
    sum = 0
    for i in range(0,n):
        sum += p[i]
    return sum

def update(state):

    # state is in form of {E(t), T(t), C(t), t(time)}
    rt, re = random(), random()
    alpha = d*state[0]+a*state[1]*(1-b*state[1])+(p*state[0])/(g+state[1])*state[1]+k_1*state[0]*state[1]+k_n1*state[2]
    +k_2*state[2]+k_3*state[2]
    dt = -np.log(rt)/alpha
    state[3] = state[3]+dt
    list = np.divide([d*state[0],a*state[1]*(1-b*state[1]),(p*state[0])/(g+state[1])*state[1],k_1*state[0]*state[1],
                  k_n1*state[2],k_2*state[2],k_3*state[2]], alpha)
# add a while/for loop
    if re < first_n_sum(list,1):
        state[0] -= 1
    elif re < first_n_sum(list,2):
        state[1] += 1
    elif re < first_n_sum(list,3):
        state[0] += 1
    elif re < first_n_sum(list,4):
        state[0] -= 1
        state[1] -= 1
        state[2] += 1
    elif re < first_n_sum(list,5):
        state[0] += 1
        state[1] += 1
        state[2] -= 1
    elif re < first_n_sum(list,6):
        state[0] += 1
        state[2] -= 1
    else:
        state[2] -= 1
        state[1] += 1
    print(state)

    return state


replicate = 0

extinctionTimes = []

data = []

while (replicate < replicates):

    extinct = False;

    state = [relativePopInit[0], relativePopInit[1],
             0.0, 0.0]  # initial conditions for state vector(0=CAR-T, 1=Tumor, 2=Conjugated cells, 3=time)
    life_history = []
    life_history.append(state)

    while (extinct == False) & (state[3] < cutoff):# wanted time length is 100 days, why do we set up time<cutoff=2

        # print "Time: ", state[3]

        new_state = copy(state)

        state = update(new_state)

        life_history.append(state)

        if ((state[0] == 0) & (state[1] == 0) & (state[2]==0)):
            extinct = True
            extinctionTimes.append(state[3])

   # life_history.append([0, 0, cutoff])

    data.append(life_history)

    replicate += 1

# -----------------ODE----------------------------------------
print(life_history)
print(k_1)
def column(matrix, i):
    return [row[i] for row in matrix]
plt.plot(column(life_history,3),column(life_history,0),'b',column(life_history,1),'r',column(life_history,2),'g')
plt.show()
def ode_system(t, z):
    E, T, C=z
    return [p*C/(g+T)-d_1*E-k_1*E*T+(k_n1+k_2)*C, #check the value of s
            a*T*(1-b*T)-k_1*E*T+(k_n1+k_3)*C,
            k_1*E*T-(k_n1+k_2+k_3)*C,
            ]

sol=solve_ivp(ode_system, [1,100],[2500,5*2500,0])
plt.plot(sol.t, sol.y[0],'b', sol.t, sol.y[1],'r', sol.t, sol.y[2],'g')
         # consider not use log
plt.show()
#sol=odeint(ode_system, [1e6,1e6,0],[1,100])

# ------------------Plots--------------------------


# Bin the data into replicates:


# ax.set_ylim(0,500)


# print(len(df))