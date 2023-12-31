#(1) We implemented implicit Euler method to obtain major principal stress for the conical part of the hopper in the active state
#(2) We are assuming that (likewise Jenike) wall normal stress in the conical hopper is equal to major principal stress
#(3) Major principal stress in the  needed to estimate unconfined yield strength
#(4) Major principal stress in the active mode in obtained in another function in the code

def MPS_in_passive_mode(X1,X2,Z1,Z2,N,sigmav_init):
    import numpy as np
    from curve_fitting import curve_fitting
    import math

    #"curve_fitting" function fits bulk density, effective angle of internal friction, FC and FFC vs sigma1.
    #It also does a power fit for wall friction angle vs normal stress
    #a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()

    delZ = abs(Z1 - Z2) / N                                                                                             # increment size in z direction (m)
    theta = (np.pi / 2 - np.arctan((Z2 - Z1) / (X2 - X1))) * 180 / np.pi                                                # hopper angle form vertical (degree)

    sigma1_passive = sigmav_init * np.ones(N)                                                                           # major principal stress in the passive mode (Pa)
    sigma1_passive[0] = sigmav_init                                                                                     # major principal stress in the active mode at the first node
    z_loc = np.zeros(N)                                                                                                 # increments in the vertical direction (m)

    rhob_passive = average_rhob*np.ones(N)                                                                              # bulk density (kg/m3) in the passive mode
    UYS_passive = np.ones(N)                                                                                            # unconfined yield strength or "FC" (Pa) in the passive mode

    # Eq. (6) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
    beta = (average_WFA + (180 / np.pi) * np.arcsin(
        np.sin(average_WFA * np.pi / 180) / np.sin(average_PHIE * np.pi / 180))) / 2

    # Eq. 20 and 21 in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
    # X and Z are just two parameters needed to estimate the major principal stress at the outlet of the hopper
    # in the passive state ("sigma1_passive")
    X = (2 * np.sin(math.radians(average_PHIE)) / (1 - np.sin(math.radians(average_PHIE)))) * (
            abs(np.sin(math.radians(2 * beta + theta))) / (np.sin(math.radians(theta)) + 1))
    Z = ((2 * (1 - np.cos(math.radians(beta + theta)))) * np.sin(math.radians(theta)) + np.sin(
        math.radians(beta)) * np.sin(math.radians(beta + theta)) ** 2) / (
                1 - np.sin(math.radians(average_PHIE)) * np.sin(math.radians(beta + theta)) ** 3)

    RHO = average_rhob
    PHI = average_PHIE
    for i in range(N):

        ## B is the diameter of the cone at this specific element
        B = X1 - i*(X1-X2)/N
        B = 2.0 * B



        #Eq. (19) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475.
        #sigma1_passive is the major principal stress at the outlet of the hopper in the passive state
        sigma1_passive[i] = (1 + np.sin(math.radians(PHI))) * Z * (RHO) * 9.8 * B / (
                    2 * (X - 1) * np.sin(math.radians(theta)))

        rhob_passive[i] = a[0] * sigma1_passive[i]+b[0]
        RHO = rhob_passive[i]
        PHI = a[1] * sigma1_passive[i] + b[1]
        UYS_passive[i] = a[2] * sigma1_passive[i] + b[2]

        z_loc[i] = (i + 1) * delZ


    sigma1_passive[0] = sigma1_passive[1]
    rhob_passive[0] = rhob_passive[1]
    UYS_passive[0] = UYS_passive[1]

    return sigma1_passive, UYS_passive
