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
print('7. 300L_IBC')

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


[M, F, P, ANGLE, theta_critical, z, sigmav_active, sigma1_active, sigma1_passive, UYS_active, UYS_passive, sigmaf, RH_diameter, D_arching_active, D_arching_passive,
 rhob_active_conical_outlet, rhob_passive_conical_outlet, sigma1_active_conical_outlet, sigma1_passive_conical_outlet] = stress_profile(KK, HEIGHT, RADIUS, X, Z)
# M = 0                         # No arch is predicted in the active state
# M = 1                         # An arch formation is predicted in the active state
# F = 0                         # Mass flow regime in the passive state is predicted
# F = 1                         # Funnel flow regime in the passive state is predicted
# F = 0 and P = 0               # No arch formation is predicted in the passive state
# F = 0 and P = 1               # An arch formation is predicted in the passive state
# F = 1 and P = 2               # Rathole formation is predicted in the funnel flow
# F = 1 and P = -2              # No rathole formation is predicted in the funnel flow
if (M==0):
    output_active = "No arch in the active state is predicted"
if (M==1):
    output_active = "An arch in the active state is predicted"
if (F==1 and P==2):
    output_passive = "A funnel flow with a rathole in the passive state is predicted"
if (F==1 and P==-2):
    output_passive = "A funnel flow without any rathole in the passive state is predicted"
if (F==0 and P==0):
    output_passive = "Mass flow without any arch in the passive state is predicted"
if (F==0 and P==1):
    output_passive = "Arch formation in the passive state is predicted"

#if (Q==2):
#    output_passive2 = "An arch in the passive state is predicted"
#if (Q==-2):
#    output_passive2 = "No arch forms in the passive state is predicted"
#if (Q==-1):
#    output_passive2 = "No funnel flow is predicted"

print("1. "+output_active)
print("2. "+output_passive)


# Chance of arch formation in the active state
# We will call curve fitting function to obtain fc (or UYS:unconfined yield strength)
# and rhob (bulk density) as a function of MPS. We need a[2], b[2] (because
[a, b, c, min_rhob, max_PHIE, max_WFA] = curve_fitting()
rhob_active_outlet =a[0]*sigma1_active[-1]+b[0]                                                                         # Bulk density as a function of MPS in the active mode
rhob_passive_outlet = a[0]*sigma1_passive[-1]+b[0]                                                                      # Bulk density as a function of MPS in the passive mode
WFA = a[4]*sigma1_active[-1]**b[4]+c[4]                                                                                 # wall friction angle in the active state

print("3. The outlet diameter is", round(2*X[-1]*1000,1), "mm and the angle from the vertical for the outlet is", round(ANGLE,1), "degrees")
print("4. The vertical load at the end of the last conical part of the hopper is", round(sigma1_active_conical_outlet,1), "Pa")
print("5. Bulk density at the end of the last conical part of the hopper in the active state is", round(rhob_active_conical_outlet,1), "kg/m3")
print("6. Bulk density at the end of the last conical part of the hopper in the passive state is", round(rhob_passive_conical_outlet,1), "kg/m3")
print("7. Critical mass flow angle is ", round(theta_critical,1),"degrees")
print("8. Arching diameter for the active state is", round(D_arching_active,1), "mm")
print("9. Arching diameter for the passive state is", round(D_arching_passive,1), "mm")
print("10. Rathole diameter is (relevant for funnel flow scenario only)", round(1000*RH_diameter,1), "mm")
print("11. Major principal stress in the passive state at the last conical part of the hopper is", round(sigma1_passive_conical_outlet,1), "Pa")

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
          +"\n The critical mass flow angle is %0.1f" %theta_critical + " while the hopper angle from the vertical is %0.1f" %ANGLE, \
          fontsize=14)
plt.show()

L=10
plt.plot(sigmav_active[::L],z[::L],'b-.',sigma1_active[::L],z[::L],'go', UYS_active[::L], z[::L], 'Pr',sigmaf[::L],z[::L],'cv',markersize=10, markerfacecolor='none', linewidth=2)
plt.xlabel("stress (Pa)",fontsize=22)
plt.ylabel("height (m)",fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(["vertical load", "major principal stress in the active state","unconfined yield strength in the active state","circumferential stress"], fontsize=22)
plt.show()






