
def MassFlow_or_FunnelFlow(X1, X2, Z1, Z2, sigmav, sigma1, PHIE, rhob, WFA, UYS, sigmaf, N, number,RADIUS):

    import numpy as np
    from curve_fitting import curve_fitting
    import math

    # Eq. (5) and Eq. (6) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
    if (WFA[-1]>PHIE[-1]):
        WFA[-1] = PHIE[-1]          #To avoid nan in the estimation of beta in the next line
    beta = (WFA[-1] + (180 / np.pi) * np.arcsin(np.sin(WFA[-1] * np.pi / 180) / np.sin(PHIE[-1] * np.pi / 180))) / 2
    theta_critical = 90 - 0.5 * (180 / np.pi) * np.arccos((1 - np.sin(PHIE[-1] * np.pi / 180)) / (2 * np.sin(PHIE[-1] * np.pi / 180))) - beta

    ## angle of the outlet of the hopper from a vertical line
    # (I added 0.000000000000001 to the denominator to avoid having the denominator equal to zero)
    theta = 90 - abs(np.arctan((Z2 - Z1) / (X2 - X1 + 0.000000000000001))) * 180 / np.pi

    F = -1                                          # Initialization of F: Later, F=1 is funnel flow and F=0 is mass flow in the passive state
    P = -1                                          # Initialization of P: Later, If F=0, P=0 means no arch forms while P=1 is equivalent to arch formation
    if (theta>theta_critical):
        F = 1                                       # funnel flow
    else:
        F = 0                                       # mass flow

    if (F == 0):                                    # if there is a mass flow, determine if there is a risk of arch formation

        # B is the outlet diameter of the hopper (m)
        B = 2 * X2
        #B = 2*RADIUS

        # Eq. 20 and 21 in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
        # X and Z are just two parameters needed to estimate the major principal stress at the outlet of the hopper
        # in the passive state ("sigma1_out")
        X = (2 * np.sin(math.radians(PHIE[-1])) / (1 - np.sin(math.radians(PHIE[-1])))) * (
                    np.sin(math.radians(2 * beta + theta)) / (np.sin(math.radians(theta)) + 1))
        Z = ((2 * (1 - np.cos(math.radians(beta + theta)))) * np.sin(math.radians(theta)) + np.sin(
            math.radians(beta)) * np.sin(math.radians(beta + theta)) ** 2) / (
                        1 - np.sin(math.radians(PHIE[-1])) * np.sin(math.radians(beta + theta)) ** 3)

        #Eq. (19) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475.
        #sigma1_out is the major principal stress at the outlet of the hopper in the passive state
        sigma1_out = (1 + np.sin(math.radians(PHIE[-1]))) * Z * (rhob[-1]) * 9.8 * B / (
                    2 * (X - 1) * np.sin(math.radians(theta)))

        # UYS_out is the unconfined yield strength (or cohesive strength) of the powder at the outlet of the hopper
        [a, b, c] = curve_fitting()
        UYS_out = a[2]*sigma1_out + b[2]

        # Eq. (2) and Eq. (3): sigma is the stress on the abutment of the hopper caused by gravity in the case of mass flow
        H = (130+theta)/65
        sigma = 9.8*rhob[-1]*B/H

        if (sigma>UYS_out):
            P = 0                           # no arch
        else:
            P = 1                           # arch formation

    ## If there is a funnel flow, determine if there is a risk of rathole formation in the ENTIRE depth of the powder bed
    # in the hopper
    if (F == 1):
        for i in range(number*N):
            if (sigmaf[i]<UYS[i]):
                P = 2                   #rathole formation
                break
            else:
                P = -2                  # no rathole formation

    return F, P, theta, theta_critical