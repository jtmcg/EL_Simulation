import math

import numpy as np
import scipy as sp
from scipy import optimize

#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab


A = 20                              #Dimension of Nodes
R = 2.0                            #R Defines the resistance of the Molly
R_t = 10.0                          #R_T Defines the resistance of the TCO
M = 1000                            #Number of iterations to calculated the steady state

#Diode Properties
n = 1.5                             #Diode Quality Factor
J_o = 10.0 ** (-11.0)               #Diode Dark Current
Rsh = 1000.                         #Diode Conductance
r = 1.0                             #Diode Series Resistance
T = 298.0                           #Ambient Temperature
k = 1.381 * 10.0 ** (-23.0)         #Boltzmann's Constant
q = 1.602 * 10.0 ** (-19.0)         #Elementary Charge
I_o = 100.0

Jsc = 0.0                           #AM 1.5 Current (mA/cm^2)

V_th = 0.0256796821 * 1000

w = [J_o, r, Rsh, n, 0.0, 0.0,
     0.0]         #w[J_o, Series Resistance, Shunt Resistance, Diode Quality, Node[X-1], Node[X+1], Node[above]]
#     0 ,1, 2 ,3, 4 , 5 , 6

Input = 5.0
Ground = 0.0

TCO = np.zeros(A)
Molly = np.zeros(A)

Molly[0] = Input

count = 0

G = lambda V: (w[4] + w[5] - V) / R - w[0] * (-1 + math.exp(((V - w[6]) - w[1] * ((w[4] + w[5] - V) / R)) / (V_th * w[3]))) - (V - w[6]) / w[2]

K = lambda V: (w[4] + w[5] - V) / R_t + w[0] * (-1 + math.exp(((w[6] - V) - w[1] * ((w[4] + w[5] - V) / R)) / (V_th * w[3]))) + (w[6] - V) / w[2]

def F(left, right, top):
    return lambda V: (left + right - 2 * V) / R - J_o * (-1 + math.exp(((V - top) - r * ((left + right - 2 * V) / R)) / (V_th * n))) - (V - top) / Rsh

def L(left, right, bottom):
    return lambda V: (left + right - 2 * V) / R_t + J_o * (-1 + math.exp(((bottom - V) - r * ((left + right - 2 * V) / R)) / (V_th * n))) + (bottom - V) / Rsh

def G(left, right, top):
    return lambda V: (left + right - V) / R - J_o * (-1 + math.exp(((V - top) - r * ((left + right - V) / R)) / (V_th * n))) - (V - top) / Rsh

def K(left, right, bottom):
    return lambda V: (left + right - V) / R_t + J_o * (-1 + math.exp(((bottom - V) - r * ((left + right - V) / R)) / (V_th * n))) + (bottom - V) / Rsh


counter = 0

while counter < M:
    count = 1
    while count < A - 1:
        Molly[count] = optimize.fsolve(F(Molly[count - 1], Molly[count + 1], TCO[count - 1]), Molly[count])[0]
        TCO[count] = optimize.fsolve(L(TCO[count - 1], TCO[count + 1], Molly[count + 1]), TCO[count])[0]
        count += 1



    Molly[A - 1] = optimize.fsolve(G(Molly[A-2], 0, TCO[A-2]), Molly[A - 1])[0]
    TCO[0] = optimize.fsolve(K(0, TCO[1], Molly[1]), TCO[0])[0]

    counter += 1

x = np.arange(0, A, 1)
#plt.plot(x, Molly, "b-")
#plt.plot(x, TCO, "r--")
print Molly
print TCO

for i in range(0, 19):
    print Molly[i+1] - TCO[i]



#plt.savefig("../Initialization.png", dpi=48)
#plt.show()