## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
##def stress_in_active_state(a, b, c, X, Y, UPPER, LOWER, HEIGHT):

import numpy as np
from Janssen_Equation import Janssen_Equation
from Motzkus_Equation import Motzkus_Equation
import matplotlib.pyplot as plt

Y = [1.4, 0.7, 0]
X = [0.518, 0.518, 0.2]

HEIGHT = 1.1879
RADIUS = 0.518

XX = [x for x, y in zip(X, Y) if y <= HEIGHT]
YY = [y for y in Y if y <= HEIGHT]
#print(XX)
#print(YY)

YY = [HEIGHT] + YY
XX = [RADIUS] + XX
print(YY)
print(XX)

# 0 -> rhob
# 1 -> PHIE
# 2 -> FC
# 3 -> PHILIN
# 4 -> wall friction

N = 10000                                                               # number of mesh in the vertical section
number = len(XX) - 1                                                    # number of sections in hopper (a hopper may have several cylinderical and cone parts)
sigmav = np.zeros(number*N)                                             # We initialized sigma0 (consolidation stress) for the top of the powder as 0.
z = np.zeros(number*N)                                                  # From the top to the bottom direction (m)
delY = np.zeros(number*N)                                               # increment length, which varies at each section of the hopper (m)
for i in range(number):
    X1 = XX[i]
    X2 = XX[i + 1]
    Y1 = YY[i]
    Y2 = YY[i + 1]
    if (X1 / X2 == 1):                                                  # If the ratio of the points remains as 1, it means that we are dealing with a cylinderical part and we will use Janssen equation to estimate the consolidation stress vs height
        sigmav_init = sigmav[i*N-i]
        print(sigmav[i*N-1])
        sigmav_out=Janssen_Equation(X1,X2,Y1,Y2,N,sigmav_init)
        sigmav[i*N:(i+1)*N] = sigmav_out[:N]
        print(sigmav_out[N-1])
    else:
        sigmav_init = sigmav[i*N-i]
        print(sigmav[i*N-1])
        sigmav_out=Motzkus_Equation(X1,X2,Y1,Y2,N,sigmav_init)
        sigmav[i*N:(i+1)*N] = sigmav_out[:N]
        print(sigmav_out[N - 1])


plt.plot(sigmav,'--')
#plt.xlabel("vertical stress (pa)",fontsize=16)
#plt.ylabel("from the top to the bottom of the hopper (m)",fontsize=16)
plt.show()

