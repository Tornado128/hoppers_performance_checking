## We implemented implicit Euler method to solve Janssen equation numerically
## The details for Janssen equation is given in "Schulze, Dietmar, and Dietmar Schulze. "Flow properties of bulk solids." Powders and Bulk solids: Behavior, characterization, storage and flow (2021): 57-100. Page 259"
## Implicit methods are always stable although they are often slower in coverging compared to explicit methods
def Janssen_Equation(X1,X2,Y1,Y2,N,sigmav_init):
    import numpy as np
    from curve_fitting import curve_fitting                                                                                             #This function fits bulk density, effective angle of internal friction, FC and FFC vs sigma1
    import matplotlib.pyplot as plt

    #X1 = 0.2845; X2 = 0.2845; Y1 = 0.4372; Y2 = 0.325
    D = 2*X1                                                                            # diameter of this section of the hopper
    g = 9.8                                                                             # gravity (m/s2)
    K = 0.4                                                                             # horizontal stress to vertical stress ratio (it is typically between 0.3 to 0.6 in hoppers and silos)
    #This function fits bulk density, effective angle of internal
    # friction, FC and FFC vs sigma1.It also does a power fit for WFA vs normal stress
    # a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c] = curve_fitting()

                                                              # diameter of the cylinderical part of the hopper (m): used in Janssen equation


    delY = (Y1 - Y2)/N                                                  # increment size in z direction (m)
    sigmav=0.01*np.ones(N+1)                                            # powder load in the vertical section of the hopper (pa) (it varies with z)
    z_loc = np.zeros(N+1)                                               # increments in the vertical direction (m)
    rhob=np.zeros(N)                                                    # bulk density (kg/m3)
    PHIE=np.zeros(N)                                                    # effective angle of internal friction (degree)
    UYS=np.zeros(N)                                                     # unconfined yield strength or "FC" (pa)
    PHILIN=np.zeros(N)                                                  # linearized angle of internal friction
    WFA = np.zeros(N)                                                   # wall friction angle
    for i in range(N):

        ## parameters for the Janssen equation: parameters needed to estimate an initial guess for sigmav for implementation of implicit Euler method
        rhob[i] = a[0]*sigmav[i]**b[0]+c[0]
        PHIE[i] = a[1] * sigmav[i] ** b[1] + c[1]
        UYS[i] = a[2]*sigmav[i] + b[2]
        PHILIN[i] = a[3]+sigmav[i] + b[3]
        WFA[i] = a[4]*sigmav[i] ** b[4] + c[4]
        if (WFA[i]>85):
            WFA[i] = 85                                                                                                 # Wall friction angle can not be above 85 degrees; otherwise it is a big number
        sigmav_guess = (g * rhob[i] - (4 * K / D) * np.tan(WFA[i] * np.pi / 180) * sigmav[i]) * delY + sigmav[i]        # We use explicit euler method to estimate sigma0 provide an initial guess for the implicit euler method!
        j = 0                                                                                                           # numerator for the while loop
        error_percent = 100                                                                                             # the goal is to bring the error percent to below 0.1 percent for the convergence
        # while loop checks that error percent (between our initial guess and the estimation) is below 0.1%
        while (error_percent>0.1 and j<100):
            rhob[i] = a[0] * sigmav_guess ** b[0] + c[0]
            PHIE[i] = a[1] * sigmav_guess ** b[1] + c[1]
            UYS[i] = a[2] * sigmav_guess + b[2]
            PHILIN[i] = a[3] + sigmav_guess + b[3]
            WFA[i] = a[4] * sigmav_guess ** b[4] + c[4]
            sigmav[i+1] = (g * rhob[i] - (4 * K / D) * np.tan(WFA[i] * np.pi / 180) * sigmav_guess) * delY + sigmav[i]
            error_percent = 100*abs( (sigmav_guess - sigmav[i + 1]) / sigmav_guess)
            sigmav_guess = sigmav[i+1]

            j = j + 1
        z_loc[i+1] = (i+1) * delY
    #plt.plot(z_loc, sigmav, 'b-')
    #plt.xlabel("(m)", fontsize=16)
    #plt.ylabel("v", fontsize=16)
    #plt.show()
    return sigmav







