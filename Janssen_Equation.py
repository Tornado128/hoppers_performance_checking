# (1) We implemented implicit Euler method to solve Janssen equation numerically in the cylinderical part of the hopper.
# (2) This function gives us vertical stress distribution (sigmav) and major principal stress (sigma1) in the cylinderical
# sections of the hoppers
# (3) The details for Janssen equation is given in "Dietmar Schulze. Flow properties of bulk solids. Powders and Bulk solids: Behavior, characterization, storage and flow (2021): 57-100. Page 259"
# (4) Implicit methods are always stable although they are often slower in covergence compared to explicit methods
def Janssen_Equation(KK,X1,X2,Z1,Z2,N,sigmav_init):
    import numpy as np
    from curve_fitting import curve_fitting                                                                                             #This function fits bulk density, effective angle of internal friction, FC and FFC vs sigma1

    D = 2*X1                                                                            # diameter of this section of the hopper
    g = 9.8                                                                             # gravity (m/s2)

    #This function fits bulk density, effective angle of internal
    # friction, FC and FFC vs sigma1. It also does a power fit for WFA vs normal stress
    # a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()

    delZ = (Z1 - Z2)/N                                                  # increment size in z direction (m)
    sigmav=sigmav_init*np.ones(N)                                       # vertical load in the vertical section of the hopper (pa) (it varies with z)
    sigmav[0] = sigmav_init                                             # vertical load at the top of the cylinderical part of the hopper
    sigma1 = 0.1*np.ones(N)                                             # major principal stress (pa)
    sigmaf_o = 0.1*np.ones(N)                                           # stress in the abutment (only for the case of funnel flow: Eq. (23) of the reference)
    z_loc = np.zeros(N)                                                 # vertical direction from the top to the bottom (m)

    rhob=np.zeros(N)                                                    # bulk density (kg/m3)
    PHIE=np.zeros(N)                                                    # effective angle of internal friction (degree)
    UYS=np.zeros(N)                                                     # unconfined yield strength or "FC" (pa)
    PHILIN=np.zeros(N)                                                  # linearized angle of internal friction
    WFA = np.zeros(N)                                                   # wall friction angle
    for i in range(N):

        ## B is the diameter of the cone at this specific element
        B = X1 - i*(X1-X2)/N
        B = 2.0 * B

        ## parameters for the Janssen equation: parameters needed to estimate an initial guess for sigmav for implementation of implicit Euler method
        if i == 0:
            rhob[i] = average_rhob
            PHIE[i] = average_PHIE
            WFA[i] = average_WFA
        else:
            rhob[i] = a[0] * sigmav[i-1] + b[0]
            PHIE[i] = a[1] * sigmav[i-1] + b[1]
            WFA[i] = a[4] * (sigmav_guess*KK) ** b[4] + c[4]


        #UYS[i] = a[2]*sigmav[i] + b[2]
        #PHILIN[i] = a[3]*sigmav[i] + b[3]
        if (WFA[i]>85):
            WFA[i] = 85                                                                                                 # Wall friction angle can not be above 85 degrees; otherwise it is a big number

        sigmav_guess = (g * rhob[i] - (4 * KK / D) * np.tan(WFA[i] * np.pi / 180) * sigmav[i]) * delZ + sigmav[i]       # See Eq. (9.7) in Dietmar Schulze. Flow properties of bulk solids. Powders and Bulk solids: Behavior, characterization, storage and flow (2021)
                                                                                                                        # We use explicit euler method to estimate sigma0 provide an initial guess for the implicit euler method!
        j = 0                                                                                                           # numerator for the while loop
        error_percent = 100                                                                                             # the goal is to bring the error percent to below 0.1 percent for the convergence
        # while loop checks that error percent (between our initial guess and the estimation) is below 0.1%
        while (error_percent>0.1 and j<100 and i<N-1):

            sigmav[i+1] = (g * rhob[i] - (4 * KK / D) * np.tan(WFA[i] * np.pi / 180) * sigmav_guess) * delZ + sigmav[i]
            error_percent = 100*abs( (sigmav_guess - sigmav[i + 1]) / sigmav_guess)
            sigmav_guess = sigmav[i+1]

            rhob[i] = a[0] * sigmav_guess + b[0]
            PHIE[i] = a[1] * sigmav_guess + b[1]
            UYS[i] = a[2] * sigmav_guess + b[2]
            PHILIN[i] = a[3] * sigmav_guess + b[3]
            WFA[i] = a[4] * (sigmav_guess*KK) ** b[4] + c[4]

            z_loc[i + 1] = (i + 1) * delZ

            j = j + 1

    sigmav[N-1] = sigmav[N-1]
    z_loc[N-1] = (N-1) * delZ
    sigma1[:N] = sigmav[:N]                                                                                             #For cylinderical part of a hopper, major principal stress (sigma1) is equal to the vertical stress

    sigma1[0] = sigma1[1]
    rhob[0] = rhob[1]
    UYS[0] = UYS[1]

    sigma1[-1] = sigma1[-2]
    rhob[-1] = rhob[-2]
    UYS[-1] = UYS[-2]

    return sigmav, sigma1, UYS
    #for i in range(N):
    #    # These three lines are used only for evaluation of rathole for the case of the formation of funnel flow in the passive state if the outlet section is vertical (like in Piccola)
    #    PHILIN_p = a[3] * sigma1[i] + b[3]
    #    # Jenike Bulletin 123 P67
    #    G = -6.86712 + 0.58911*PHILIN_p-0.012966*PHILIN_p**2.0+0.00011939*PHILIN_p**3.0

    #    # stress in the abutment (only for the case of funnel flow: Eq. (23) of the reference)
    #    # Some hoppers like Piccola have cylinderical outlets. That is why we put RADIUS of the outlet here
    #    # remeber B (diameter) is constant in the cylinderical part of the hopper
    #    sigmaf_o[i] = rhob[i] * g * B / G

    #RH_diameter = G * UYS[-1] / (rhob[-1] * g)                                                                          # Rathole diameter (m)
    #return sigmav, sigma1, sigmaf_o, UYS, RH_diameter







