## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
def stress_profile(HEIGHT, RADIUS, X, Z):

    import numpy as np
    from Janssen_Equation import Janssen_Equation
    from Motzkus_Equation import Motzkus_Equation
    from MassFlow_or_FunnelFlow import MassFlow_or_FunnelFlow

    XX = [x for x, z in zip(X, Z) if z < HEIGHT]                            # We eliminate the x points above the powder level
    ZZ = [z for z in Z if z < HEIGHT]                                       # We eliminate the z points above the powder level

    ZZ = [HEIGHT] + ZZ                                                      # Height of the powder in the hopper
    XX = [RADIUS] + XX                                                      # x-location corresponds to the height of the powder in the hopper

    N = 10000                                                               # number of mesh in the vertical section
    number = len(XX) - 1                                                    # number of sections in hopper (a hopper may have several cylinderical and cone parts)
    sigmav = np.zeros(number*N)                                             # We initialized sigma0 (vertical stress) for the top of the powder as 0.
    sigma1 = np.zeros(number*N)                                             # We initialized sigma1 (major principal stress) for the top of the powder as 0.
    z = np.zeros(number*N)                                                  # (Hieght) From the top to the bottom direction (m)
    sigmaf = np.zeros(number*N)                                             ## sigmaf is used ONLY for the case of estimation of stress on the abutment in the case of funnel flow
                                                                            # for passive state. sigmaf is given in Eq. (22) of the reference.
    UYS = np.zeros(number*N)                                                # unconfined yield strength (pa)

    numerator = 0                                                           # numerator, which is equal to the total number of elements in z direction
    # This loop goes through all the cylinderical and cone parts of the hopper.
    # If the section is cylinderical, it uses Janssen equation to obtain vertical stress distribution.
    # if the section is cone shaped, it uses Motzkus equarion to obtain vertical stress distribution.
    for i in range(number):
        X1 = XX[i]
        X2 = XX[i + 1]
        Z1 = ZZ[i]
        Z2 = ZZ[i + 1]
        for j in range(N):
            z[numerator] = Z1 - (Z1 - Z2) * j / N
            numerator = numerator + 1
        if (X1 / X2 == 1):                                                  # If the ratio of the points remains as 1, it means that we are dealing with a cylinderical part and we will use Janssen equation to estimate the consolidation stress vs height
            sigmav_init = sigmav[i*N-i]
            [sigmav_o, sigma1_o, WFA, PHIE, rhob, sigmaf_o, UYSf_o]=Janssen_Equation(X1,X2,Z1,Z2,N,sigmav_init, RADIUS)
            #print(sigmav_o[-1])
            sigmav[i*N:(i+1)*N] = sigmav_o[:N]
            sigma1[i*N:(i+1)*N] = sigma1_o[:N]
            sigmaf[i * N:(i + 1) * N] = sigmaf_o[:N]                                  # sigmaf is useful only for the funnel flow case to evaluate the possibility of ratholing (Eq. (23) of the reference)
            UYS[i * N:(i + 1) * N] = UYSf_o[:N]
        else:                                                               # We are in the cone part of the hopper
            sigmav_init = sigmav[i*N-i]
            [sigmav_o, sigma1_o, WFA, PHIE, rhob, sigmaf_o, UYSf_o]=Motzkus_Equation(X1,X2,Z1,Z2,N,sigmav_init, RADIUS)
            #print(sigmav_o[-1])
            sigmav[i*N:(i+1)*N] = sigmav_o[:N]
            sigma1[i*N:(i+1)*N] = sigma1_o[:N]
            sigmaf[i * N:(i + 1) * N] = sigmaf_o[:N]                                  # sigmaf is useful only for the funnel flow case to evaluate the possibility of ratholing (Eq. (23) of the reference)
            UYS[i * N:(i + 1) * N] = UYSf_o[:N]

    ##WFA_out = WFA[-1]                                                       # wall friction angle at the outlet
    ##PHIE_out = PHIE[-1]                                                     # effective angle of internal friction at the outlet
    ##rhob_out = rhob[-1]
    # We want to determine if we are dealing with a mass flow or a funnel flow
    [F, P] = MassFlow_or_FunnelFlow(X1, X2, Z1, Z2, sigmav, sigma1, PHIE, rhob, WFA, UYS, sigmaf, N, number, RADIUS)
    return z, sigmav, sigma1, F, P


