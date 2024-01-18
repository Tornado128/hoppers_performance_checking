##(1) Motzkus equation is used to obtain the vertical stress distribution (sigmav) and then major principal stress (sigma1)
## in the conical part of the hopper in the active mode.
##(2) The equations in this function are mostly obtain from "A proposed complete methodology to predict gravity flow obstruction
# of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475. However, some equations
# are also obtained from "Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Handling, 1994, 14, 497-503"
##(3) We implemented implicit Euler method to solve Motzkus equation (Eq. (7) of Leung et al) numerically
##(4) Implicit methods are always stable although they are often slower in coverging compared to explicit methods
##(5) The reason for picking the numerical approach to solve Eq. (7) of the above reference is that the parameters like
##wall friction angle, bulk density, ... are changing in the vertical direction because of the increasing load

def properties_in_the_conical_part_in_active_mode(X1,X2,Z1,Z2,N,sigmav_init):
    import numpy as np
    from curve_fitting import curve_fitting
    import math
    from D_arch_active_estimation import D_arch_active_estimation

    # imaginary position of the apex of the cone (see Figure 4 of the reference)
    Z_apex = ((Z2-Z1)/(X2-X1))*(0-X1)+Z1

    # gravity (m/s2)
    g = 9.8

    #"curve_fitting" function fits bulk density, effective angle of internal friction, FC and FFC vs sigma1.
    #It also does a power fit for wall friction angle vs normal stress
    #a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()


    theta = (np.pi/2 - np.arctan((Z2-Z1)/(X2-X1)))*180/np.pi                                                            # hopper vertical angle (degree)
    delZ = (Z1 - Z2) / N                                                                                                # increment size in z direction (m)
    sigmav = sigmav_init*np.ones(N)                                                                                     # vertical stress in the cone section of the hopper (pa) (it varies with z)
    sigmav[0] = sigmav_init                                                                                             # vertical stress (pa) at the top of the cone section of the hopper.
    sigma1_active = 0.1 * np.ones(N)                                                                                    # major principal stress (sigma1) in the active mode (Pa)
    z_loc = np.zeros(N)                                                                                                 # increments in the vertical direction (m)
    sigmaf_o = 0.1 * np.ones(N)                                                                                         # stress in the to have rathole (only for the case of funnel flow: Eq. (23) of the reference)

    rhob = average_rhob * np.ones(N)                                                                                    # bulk density (kg/m3)
    PHIE = average_PHIE * np.ones(N)                                                                                    # effective angle of internal friction (degree)
    WFA = average_WFA * np.ones(N)                                                                                      # wall friction angle (degree)
    UYS = np.zeros(N)                                                                                                   # unconfined yield strength or "FC" (pa)
    PHILIN = np.zeros(N)                                                                                                # linearized angle of internal friction (degree)
    for i in range(N):

        ## B is the diameter of the cone at this specific element
        B = X1 - i*(X1-X2)/N
        B = 2.0 * B

        # parameters for the Motzkus equation: parameters needed to estimate an initial guess for sigmav for
        # implementation of implicit Euler method
        if (i > 0):
            rhob[i] = a[0]*sigma1_active[i-1]+b[0]
            PHIE[i] = a[1] * sigma1_active[i-1] + b[1]
            WFA[i] = a[4]*sigma1_active[i-1] ** b[4] + c[4]

        if (WFA[i] > 85):
            WFA[i] = 85                                 # Wall friction angle can not be above 85 degrees; otherwise it is a big number

        if (WFA[i]>PHIE[i]):
            WFA[i] = PHIE[i]                            # Wall friction angle may be unphysically larger than effective angle of internal friction at very low stresses (close to the top of the powder)

        ## Eq. (12) of Leung et al "A proposed complete methodology to predict gravity flow obstruction
        # of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475.
        lambda_if = (1 - np.sin(math.radians(WFA[i])) ** 2.0 - np.sqrt(
            (1 - np.sin(math.radians(WFA[i])) ** 2) * (np.sin(math.radians(PHIE[i])) ** 2 - np.sin(math.radians(WFA[i])) ** 2))) / (
                            1 + np.sin(math.radians(WFA[i])) ** 2.0 + np.sqrt((1 - np.sin(math.radians(WFA[i])) ** 2) * (
                            np.sin(math.radians(PHIE[i])) ** 2 - np.sin(math.radians(WFA[i])) ** 2)))

        ## Eq. (15) of Leung et al "A proposed complete methodology to predict gravity flow obstruction
        # of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475.
        lambda_f = (np.tan(math.radians(90 - theta)) - np.tan(math.radians(WFA[i])) - lambda_if * np.tan(math.radians(WFA[i])) * (1 + np.tan(math.radians(WFA[i])) ** 2.0 -(np.tan(math.radians(90 - theta)) - np.tan(math.radians(WFA[i]))) ** 2)) / (np.tan(math.radians(90 - theta)) * (1 + np.tan(math.radians(WFA[i])) * np.tan(math.radians(90 - theta))))


        #lambda_f = (np.tan(math.radians(90 - theta)) - mu_if) - lambda_if * np.tan(math.radians(WFA[i])) * (1 + np.tan(math.radians(WFA[i])) ** 2.0 -(np.tan(math.radians(90 - theta)) - np.tan(math.radians(WFA[i]))) ** 2)) / (np.tan(math.radians(90 - theta)) * (1 + np.tan(math.radians(WFA[i])) * np.tan(math.radians(90 - theta))))
        ## see Eq. (23) and Eq. (21) of Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Hand;img, 1994, 14, 497-503
        mu_if = np.tan(math.radians(WFA[i]))
        mu_f = mu_if * lambda_if/lambda_f

        if (theta > 90 - math.degrees(np.arcsin(np.sin(math.radians(WFA[i])) / np.sin(math.radians(PHIE[i]))))):        # material failure
            ## Eq. (18) in Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Handling, 1994, 14, 497-503
            #K = 0.5 * (1 + lambda_if) - 0.5 * (1 - lambda_if) * np.cos(2 * math.radians(theta)) + mu_if*lambda_if *np.sin(2*math.radians(theta))

            # Eq. (19) in Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Hand;img, 1994, 14, 497-503
            n = 2 * mu_if * lambda_if * np.tan(math.radians(90 - theta))

        else:                                                                                                           # wall failure

            ## Eq. (16) in Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Handling, 1994, 14, 497-503
            #K = 0.5 * (1 + lambda_f) - 0.5 * (1 - lambda_f) * np.cos(2 * math.radians(theta)) + mu_f * lambda_f * np.sin(2*math.radians(theta))

            # Eq. (17) in Dietmar Schulze, "The prediction of initial stresses in hoppers", Bulk Solids Hand;img, 1994, 14, 497-503
            n = 2 * mu_f * lambda_f * np.tan(math.radians(90 - theta))

        sigmav_guess = (g * rhob[i] - n*sigmav[i]/((Z1-Z_apex-z_loc[i])))*delZ + sigmav[i]                              # We use explicit euler method to estimate sigma0 provide an initial guess for the implicit euler method!
        if (i>0):
            sigmav_guess = (g * rhob[i] - n * sigmav[i] / ((Z1 - Z_apex - z_loc[i]))) * delZ + sigmav[i]                # We use explicit euler method to estimate sigma0 provide an initial guess for the implicit euler method!

        j = 0                                                                                                           # numerator for the while loop
        error_percent = 100                                                                                             # the goal is to bring the error percent to below 0.1 percent for the convergence
        # while loop checks that error percent (between our initial guess and the estimation) is below 0.1%
        while (error_percent>0.1 and j<100 and i<N-1):

            sigmav[i+1] = (g * rhob[i] - n*sigmav_guess/((Z1-Z_apex-z_loc[i])))*delZ + sigmav[i]
            error_percent = 100*abs( (sigmav_guess - sigmav[i + 1]) / sigmav_guess)

            sigmav_guess = sigmav[i+1]
            z_loc[i + 1] = (i + 1) * delZ

            j = j + 1

        # It is the same as Eq. (14) of Leung et al "A proposed complete methodology to predict gravity flow obstruction
        # of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475.
        K_for_MPS_calculation = 0.5 * (1 + lambda_f) - 0.5 * (1 - lambda_f) * np.cos(2 * math.radians(theta)) + mu_f * lambda_f * np.sin(2 * math.radians(theta))

        # Eq. (18) of the reference Leung et al "A proposed complete methodology to predict gravity flow obstruction
        # of pharmaceutical powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475.
        beta_p = math.degrees(np.arctan(((1 + np.cos(2 * math.radians(theta))) * np.tan(math.radians(WFA[i]))) / (1 + np.sin(2 * math.radians(theta)) * np.tan(math.radians(WFA[i])) - 1 / K_for_MPS_calculation)))

        #print( ((1 + np.cos(2 * math.radians(theta))) * np.tan(math.radians(WFA[i]))) / (1 + np.sin(2 * math.radians(theta)) * np.tan(math.radians(WFA[i])) - 1 / K_for_MPS_calculation))
        # # major principal stress in the conical part in the active state: Eq. (17) of the appendix in the reference
        sigma1_active[i] = 0.5 * sigmav[i] * (1 + lambda_f) + abs(K_for_MPS_calculation * sigmav[i] - sigmav[i] * (1 + lambda_f) / 2) / np.cos(math.radians(beta_p))

        rhob[i] = a[0] * sigma1_active[i] + b[0]
        PHIE[i] = a[1] * sigma1_active[i] + b[1]
        UYS[i] = a[2] * sigma1_active[i] + b[2]
        PHILIN[i] = a[3] * sigma1_active[i] + b[3]
        WFA[i] = a[4] * sigma1_active[i] ** b[4] + c[4]

        # sigma1_active may become negative at very low consolidation stresses at the top of the powder.
        if sigma1_active[i]<0:
            sigma1_active[i] = None

    sigma1_active[0] = sigma1_active[1]
    rhob[0] = rhob[1]
    UYS[0] = UYS[1]

    sigma1_active[-1] = sigma1_active[-2]
    rhob[-1] = rhob[-2]
    UYS[-1] = UYS[-2]

    sigmav[N-1] = sigmav[N-2]
    z_loc[N-1] = (N-1) * delZ

    # arch diameter in the active state
    D_arching_active = D_arch_active_estimation(UYS, rhob, theta)

    return sigmav, sigma1_active, UYS, D_arching_active