## The goal of this function is to determine the height of the material in the hopper at a given filling percent

#Y = [0.518, 0.478, 0.313, 0.125, 0]
#X = [0.125, 0.125, 0.025, 0.025, 0.018]
#fill_percent = 25
def height_position(X, Y, fill_percent):
    import numpy as np
    import math
    number = len(X) - 1                                                                         # number of sections of the vessel. Typically a hopper have two sections: one is cylinder and one is the cone. However, some hoppers may have several cones and cylinderical parts
    N = 1000                                                                                    # dividing the height of the hopper into N elements
    m = 0                                                                                       # overall numerator
    vol = np.ones(number*N)
    height = np.zeros(number*N)                                                                 # height of each mesh element from the top to the bottom (m)
    radius = np.zeros(number*N)                                                                 # radius of at each mesh element from the top to the bottom (m)
    for i in range(number):                                                                     # We estimate the volume of each section of the hopper individually. Then we add them altohther.
        delY = abs((Y[i + 1] - Y[i]) / N)                                                       # divide this section of the hopper into N elements
        if (X[i]/X[i+1] == 1):                                                                  # If the ratio of the points remains as 1, it means that we are dealing with a cylinder
            vol_inc = (math.pi * X[i] ** 2.0) * delY                                            # incremental volume in the straight part of the hopper (m^3)
            j = 0                                                                               # local numerator for this for loop
            while (j < N):
                vol[m] = vol_inc                                                                # volume at each volume element in the straight part of the hopper
                height[m] = np.max([Y[i+1], Y[i]])-delY*(j+1)
                radius[m] = X[i]
                j = j + 1                                                                       # updating the numerator
                m = m + 1

        else:                                                                                   # Cone part of the hopper
            Y_aux = -(Y[i]-Y[i+1])*X[i+1]/(X[i]-X[i+1])+Y[i+1]                                  # apex of the bigger imaginary cone
            X_p_before = X[i]
            Y_p_before = Y[i]
            j = 0                                                                               # local numerator
            while (j < N):
                X_p_after = X[i] - (1+j) * delY * (X_p_before - X[i + 1]) / (Y_p_before - Y[i + 1])
                radius[m] = X_p_after                                                                  # X coordinate corresponds to the height of the material
                Y_p_after = Y_p_before + ((Y_p_before -Y[i + 1]) / (X_p_before - X[i + 1])) * (X_p_after-X_p_before)
                vol_aux_big = (math.pi * X_p_before ** 2) * (abs(Y_p_before - Y_aux)) / 3               # volume the bigger imaginary cone (apex is imaginary): m^3
                vol_aux_small = (math.pi * X_p_after ** 2) * (abs(Y_p_after - Y_aux)) / 3               # volume the smaller imaginary cone (apex is imaginary): m^3
                vol[m] = vol_aux_big - vol_aux_small
                X_p_before = X_p_after
                Y_p_before = Y_p_after
                height[m] = np.max([Y[i+1], Y[i]])-delY*(j+1)
                j = j + 1
                m = m + 1


    ## the goal is to find the height associated to the filling percent of the hopper
    volume = sum(vol)                                               #volume of the hopper! (m^3)                                                #reversing the list for the sake of convinience
    vol_cumulative = 100                                                #percent of the material in the hopper (it wil be updated)
    i = 0                                                               #numerator for the list
    while (i < number*N):
        vol_cumulative = vol_cumulative - 100*vol[i]/volume
        HEIGHT = height[i]
        RADIUS = radius[i]
        if vol_cumulative < fill_percent:                           # stop the while loop once we reach to the filling percent of the hopper
            break
        i = i + 1

    return HEIGHT, RADIUS
