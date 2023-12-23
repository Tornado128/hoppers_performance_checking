# This function estimates the volume of the hopper by dividing the hopper into cylinderical and cone parts.
# It also estimates the percent of the volume corresponding to each point
# It adds the volume of all parts to estimate the total volume of the hopper.

def vessel_volume(X, Y):
    import numpy as np
    import math
    number = len(X) - 1                                 # number of sections of the vessel. Typically a hopper have two sections: one is cylinder and one is the cone

    volume = 0                                          # volume of the vessel (cylinderical parts + cone parts): some vessels have several stright and cone parts
    vol=np.ones(number)                                 # defining the volume of each section
    for i in range(number):                             # We estimate the volume of each section of the hopper individually. Then we add them altohther.
        if (X[i]/X[i+1] == 1):                          # If the ratio of the points remains as 1, it means that we are dealing with a cylinder
            vol[i] = (math.pi*X[i]**2.0)*(Y[i]-Y[i+1])     # volume of a cylinder parts
        else:
            Y_aux = -(Y[i]-Y[i+1])*X[i+1]/(X[i]-X[i+1])+Y[i+1]
            vol_aux = (math.pi*X[i]**2)*(abs(Y[i]-Y_aux))/3             # volume the big imaginary cone (apex is imaginary)
            vol[i] = vol_aux - (math.pi*X[i+1]**2)*abs(Y[i+1]-Y_aux)/3     # volume of the cone is the volume of the big imaginary cone minus the volume of the small imagninary cone

        volume = volume + vol[i]                                           # volume of the hopper is the summation of cone parts and the cylinderical parts

    vol = vol[::-1]                                                         # reverse the order of the array
    percent = np.zeros(len(X))                                              # percent of the total volume at each location in the hopper
    V = 0
    for i in range(number+1):
        if i != 0:
            V = V + vol[i-1]
            percent[i] = 100*(V)/volume
        else:
            percent[number] = 0

    return volume, percent, vol



