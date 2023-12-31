## (1) This code is based on Leung, Lap Yin, et al. "A proposed complete methodology to predict gravity flow obstruction
## of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475."
## (2) The purpose of this piece of code is to assess the vertical load (sigmav), major principal stress (sigma1) and
## unconfined yeild strength (UYS) in a GIVEN conical hopper in both active (filling) and passive modes (emptying).
## (3) The stress profile later is used to determine if there is a mass flow or funnel flow in passive mode.
## (4) If there is a mass flow regime in the given hopper, the code determines if an arch forms or not
## (5) If there is a funnel flow regime in the given hopper, the code determines if a rathole forms in the hopper
## Note: Remember that Jenike method is very conservative for the estimation of arch diameter. It is even more conservative
## for the prediction of rathole diameter
import math

import matplotlib.pyplot as plt
import sys
import numpy as np
from HopperSpecification import *                                                                                       #This is a constructor where several conical hoppers are defined as objects
from is_number import *                                                                                                 #This function checks if the inputs are numbers and if they are processable
from vessel_volume import *                                                                                             #This function estimates the volume of the hoppers
from height_position import *                                                                                           #This function estimates the height ("HEIGHT") of the powder in the hopper. It also gives the radius (RADIUS) corresponds to the height of the powder in the hopper
from stress_profile import *                                                                                            #This function estimates the stress profile in both active and passive modes
from curve_fitting import *                                                                                             #This function does linear and power curve fitting on effective angle of internal friction, bulk density,

KK = 0.4                                                                                                                #Horizontal stress to vertical stress, which is typically equal to 0.3 to 0.6 for powders
print("I am assuming the horizontal to vertical stress ratio is equal to ", KK)

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
# It stops the code if the above conditions are not satisfied
[S, k, fill_percent]=is_number(k,fill_percent)
if S == 1:
    print("One of or both the numbers you entered is not processable! I have to stop the code!")
    sys.exit(1)     # stop running the code

# This function estimates the volume of the vessel
# percent is the volume percent associated with each section of the hopper
X = r[k-1].x                    # X values for the position of the vessel
Z = r[k-1].z                    # Z values (height) for the position of the vessel
[volume, percent, vol] = vessel_volume(X, Z)
# This function estimates the HEIGHT of the material in the hopper (m)
# It also gives the radius (or x-location) associated with the height of the powder
[HEIGHT, RADIUS] = height_position(X, Z, fill_percent)

# Vertical stress by the powder in both active and passive mode
# Passive state: F=1 is funnel flow and F=0 is mass flow in the passive state
# Passive state: If F=0, P=0 is no-arch and P=1 is equivalent to arch formation.
# Passive state: If F=1 and P=2, we have a funnel flow with rathole formation
# Passive state: If F=1 and P=-2, we have a funnel flow but no rathole forms
[M, F, P, theta, theta_critical, z, sigmav, sigma1_active, sigma1_passive, UYS_active, UYS_passive, sigmaf, RH_diameter] \
    = stress_profile(KK, HEIGHT, RADIUS, X, Z)

if (F==1 and P==2):
    output_passive = "We have funnel flow with a rathole in the passive state"
if (F==1 and P==-2):
    output_passive = "We have funnel flow but without any rathole in the passive state"
if (F==0 and P==0):
    output_passive = "We have mass flow without any arch formation in the passive state"
if (F==0 and P==1):
    output_passive = "We have mass flow with arch formation in the passive state"
if (M==0):
    output_active = "No arch in the active state"
if (M==1):
    output_active = "An arch forms in the active state"

print("1. "+output_active)
print("2. "+output_passive)


# Chance of arch formation in the active state
# We will call curve fitting function to obtain fc (or UYS:unconfined yield strength)
# and rhob (bulk density) as a function of MPS. We need a[2], b[2] (because

[a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()
rhob_active_outlet =a[0]*sigma1_active[-1]+b[0]                                                                         # Bulk density as a function of MPS in the active mode
rhob_passive_outlet = a[0]*sigma1_passive[-1]+b[0]                                                                      # Bulk density as a function of MPS in the passive mode
WFA = a[4]*sigma1_active[-1]**b[4]+c[4]                                                                                 # wall friction angle in the active state
theta = 90 - np.arctan((Z[-2]-Z[-1])/(X[-2]-X[-1]+0.000000000000001))*180/np.pi                                         # angle of the outlet of the hopper from a vertical line
H = (130 +theta)/65                                                                                                     # Eq. (3) of the reference: Leung et al
D_arching_active = H * UYS_active[-1] / (rhob_active_outlet*9.8)                                                        # arching diameter (m) in the active mode: calculation based on Eq. (2) and Eq. (3)
D_arching_active = D_arching_active * 1000                                                                              # arching diameter (mm) in the active mode
D_arching_passive = H * UYS_passive[-1] / (rhob_passive_outlet*9.8)                                                     # arching diameter (m) in the passive mode: calculation based on Eq. (2) and Eq. (3)
D_arching_passive = D_arching_passive * 1000                                                                            # arching diameter (mm) in the passive mode
print("3. The outlet diameter is ", round(2*X[-1]*1000,1), " mm")
print("4. The vertical load at the outlet of the hopper is ", round(sigmav[-1],1), " Pa")
print("5. Bulk density at the outlet in the active state is ", round(rhob_active_outlet,1), " kg/m3")
print("6. Bulk density at the outlet in the passive state is ", round(rhob_passive_outlet,1), " kg/m3")
print("7. Critical mass flow angle is ", round(theta_critical,1), " degrees")
print("8. Arching diameter for the active state is ", round(D_arching_active,1), " mm")
print("9. Arching diameter for the passive state is ", round(D_arching_passive,1), " mm")
print("10. Rathole diameter is (for passive state) ", round(1000*RH_diameter,1), " mm")

## showing the dimensions of the hopper
percent = percent[::-1]                                                                                                 # reversing the percent of filling order for convinience
plt.plot(r[k-1].x, r[k-1].z,'b-')
plt.plot([-x for x in r[k-1].x], r[k-1].z,'b-')
## putting label of percent of filling
for i in range(len(X) - 1):
    plt.text(X[i],Z[i],'filling%: {}'.format(round(percent[i],1)),fontsize=16,color='g')
    plt.axhline(Z[i], color='g', linestyle='--')
plt.axhline(HEIGHT, color='r', linestyle='--')
plt.xlabel("x-axis (m)",fontsize=16)
plt.ylabel("y-axis (m)",fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
volume_liter = round(1000*volume,2)                                             # m3 to liter
plt.title("The volume of %s" %r[k-1].name + f" is {volume_liter:0.2f} liter."+"\n %s" %output_active +"\n %s" %output_passive \
          +"\n The critical mass flow angle is %0.1f" %theta_critical + " while the hopper angle from the vertical is %0.1f" %theta, \
          fontsize=18)
plt.show()

L=100
plt.plot(sigmav[::L],z[::L],'b-.',sigma1_active[::L],z[::L],'go',sigma1_passive[::L],z[::L],'m*', UYS_active[::L], z[::L], 'Pr', UYS_passive[::L],z[::L],'pk', sigmaf[::L],z[::L],'cv',markersize=10, markerfacecolor='none', linewidth=2)
plt.xlabel("stress (Pa)",fontsize=22)
plt.ylabel("height (m)",fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(["vertical load", "major principal stress in the active state","major principal stress in the passive state",
            "unconfined yield strength in the active state", "unconfined yield strength in the passive state","circumferential stress"], fontsize=22)
plt.show()

'''
if (sigma > UYS[-1]):
    print("No arching in the active state because, at the outlet, external stress (%0.2f"%sigma+" pa) is larger than UYS (%0.2f"%UYS[-1]+" pa)" )
    output_active = "No arching in the active state because, at the outlet, external stress (%0.2f"%sigma+" pa) is larger than UYS (%0.2f"%UYS[-1]+" pa)"                                                                      # It will be shown in the title of the plots
else:
    print("Arching takes place in the active state because, at the outlet, external stress (%0.2f" % sigma + " pa) is smaller than UYS (%0.2f" %UYS[-1] + " pa)")
    output_active = "Arching takes place in the active state because, at the outlet, external stress (%0.2f" % sigma + " pa) is smaller than UYS (%0.2f" %UYS[-1] + " pa)"
    print("Arching diameter is " + str("{:.2f}".format(D_arching)) + " m")


print("Vertical load at the outlet is " + str("{:.2f}".format(sigmav[-1])) + " Pa")
print("MPS at the outlet is " + str("{:.2f}".format(sigma1[-1])) + " kg/m3")
print("Bulk density at the outlet is " + str("{:.2f}".format(rhob[-1])) + " kg/m3")
print("UYS at the outlet is " + str("{:.2f}".format(UYS[-1])) + " kg/m3")
print("The critical mass flow angle is " + str("{:.2f}".format(theta_critical)) + " while the hopper angle from the vertical is " + str("{:.2f}".format(theta)))

## showing the dimensions of the hopper
percent = percent[::-1]                                                                                                 # reversing the percent of filling order for convinience
plt.plot(r[k-1].x, r[k-1].z,'b-')
plt.plot([-x for x in r[k-1].x], r[k-1].z,'b-')
## putting label of percent of filling
for i in range(len(X) - 1):
    plt.text(X[i],Z[i],'filling%: {}'.format(round(percent[i],1)),fontsize=16,color='g')
    plt.axhline(Z[i], color='g', linestyle='--')
plt.axhline(HEIGHT, color='r', linestyle='--')
plt.xlabel("x-axis (m)",fontsize=16)
plt.ylabel("y-axis (m)",fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
volume_liter = round(1000*volume,2)                                             # m3 to liter
plt.title("The volume of %s" %r[k-1].name + f" is {volume_liter:0.2f} liter."+"\n %s" %output_active +"\n %s" %output_passive \
          +"\n The critical mass flow angle is %0.1f" %theta_critical + " while the hopper angle from the vertical is %0.1f" %theta, \
          fontsize=18)
plt.show()


#estimation of vertical wall stress
sigmaw = np.zeros(len(WFA))
j = 0
for i in WFA:
    sigmaw[j] = KK * np.tan(math.radians(WFA[j]-0.001)) * sigmav[j]
    j = j + 1


'''





