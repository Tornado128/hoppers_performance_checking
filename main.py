## (1) This code is based on Leung, Lap Yin, et al. "A proposed complete methodology to predict gravity flow obstruction of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475."
## (2) The purpose of this piece of code is to assess the flow regime in a conical hopper in active and passive modes

import sys
import matplotlib.pyplot as plt
from HopperSpecification import *                                                                                       #This is a constructor where several conical hoppers are defined as objects
from is_number import *                                                                                                 #This function checks if the inputs are numbers and if they are processable
from vessel_volume import *                                                                                             #This function estimates the volume of the hoppers
from height_position import *                                                                                           #This function estimates the height ("HEIGHT") of the powder in the hopper. It also gives the radius (RADIUS) corresponds to the height of the powder in the hopper

print('Choose one of the following systems by their number!')
print('1. 1000L_IBC')
print('2. 100L_IBC')
print("3. Piccola")
print('4. Courtoy')
print('5. Korsch_XL100_poor_flow')
print('6. Korsch_XL100_gravity')
k = input("Enter the number! ")                                                                                         # A number between 1 to 6 must be selected

fill_percent = input("Enter the filling percent of the hopper system! ")                                                # How much is the filling percent of the hopper vessel?

# This function checks that k and fill_percent are numbers and they are within the range.
# It stops the code if these conditions are not satisfied
[S, k, fill_percent]=is_number(k,fill_percent)
if S == 1:
    print("One of or both the numbers you entered is not processable! I have to stop the code!")
    sys.exit(1)     # stop running the code

# This function estimates the volume of the vessel
# percent is the volume percent associated with each section of the hopper
X = r[k-1].x                    # X values for the position of the vessel
Y = r[k-1].y                    # Y values for the position of the vessel
[volume, percent, vol] = vessel_volume(X, Y)

# This function estimates the HEIGHT of the material in the hopper (m)
# It also gives the radius (or x-location) associated with the height of the powder
[HEIGHT, RADIUS] = height_position(X, Y, fill_percent)

print(HEIGHT)
print(RADIUS)

## showing the dimensions of the hopper
percent = percent[::-1]                                                                                                 # reversing the percent of filling order for convinience
plt.plot(r[k-1].x, r[k-1].y,'b-')
plt.plot([-x for x in r[k-1].x], r[k-1].y,'b-')
## putting label of percent of filling
for i in range(len(X) - 1):
    plt.text(X[i],Y[i],'filling%: {}'.format(round(percent[i],1)),fontsize=16,color='g')
    plt.axhline(Y[i], color='g', linestyle='--')
plt.axhline(HEIGHT, color='r', linestyle='--')
plt.xlabel("(m)",fontsize=16)
plt.ylabel("(m)",fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
volume_liter = round(1000*volume,2)                                             # m3 to liter
plt.title("The volume of %s" %r[k-1].name + f" is {volume_liter:0.2f} liter", fontsize=18)                                                            # liter to m3
plt.show()






