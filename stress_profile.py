## This function estimates the stress profile in the hopper in active (initial discharging and filling) mode
def stress_profile(KK, HEIGHT, RADIUS, X, Z):

    import numpy as np
    from Janssen_Equation import Janssen_Equation
    from properties_in_the_conical_part_in_active_mode import properties_in_the_conical_part_in_active_mode
    from properties_in_the_conical_part_in_passive_mode import properties_in_the_conical_part_in_passive_mode
    from MassFlow_or_FunnelFlow import MassFlow_or_FunnelFlow



    XX = [x for x, z in zip(X, Z) if z < HEIGHT]                            # We eliminate the x points above the powder level
    ZZ = [z for z in Z if z < HEIGHT]                                       # We eliminate the z points above the powder level

    ZZ = [HEIGHT] + ZZ                                                      # Height of the powder in the hopper
    XX = [RADIUS] + XX                                                      # x-location corresponds to the height of the powder in the hopper

    N = 10000                                                               # number of mesh in the vertical section
    number = len(XX) - 1                                                    # number of sections in hopper (a hopper may have several cylinderical and cone parts)
    sigmav = np.zeros(number*N)                                             # We initialized sigmav (vertical stress load), Pa
    sigma1_active = np.zeros(number*N)                                      # We initialized major principal stress in active mode (vertical stress load)
    sigma1_passive = np.zeros(number*N)                                     # We initialized major principal stress in passive mode (vertical stress load)
    z = np.zeros(number*N)                                                  # (Hieght) From the top to the bottom direction (m)
    sigmaf = np.zeros(number*N)                                             ## sigmaf is used ONLY for the case of estimation of stress
                                                                            # on the abutment in the case of funnel flow for passive state.
                                                                            # sigmaf is given in Eq. (22) of the reference.
    UYS_active = np.zeros(number*N)                                         # unconfined yield strength in the active mode  (Pa)
    UYS_passive = np.zeros(number*N)                                        # unconfined yield strength in the passive mode  (Pa)

    numerator = 0                                                           # numerator, which is equal to the total number of elements in z direction
    # This loop goes through all the sections of the hopper sequentially (from the top to the bottom).
    # If the section is cylinderical, it uses Janssen equation to obtain vertical stress distribution.
    # if the section is cone shaped, it uses Motzkus equation to obtain vertical stress distribution.
    for i in range(number):
        X1 = XX[i]
        X2 = XX[i + 1]
        Z1 = ZZ[i]
        Z2 = ZZ[i + 1]
        for j in range(N):
            z[numerator] = Z1 - (Z1 - Z2) * j / N
            numerator = numerator + 1
        if (X1 / X2 == 1):                                                  # If the ratio of the points remains as 1, it means that we are dealing with a cylinderical part and we will use Janssen equation to estimate the consolidation stress vs height

            ## Janssen equation is valid for both active and passive modes to obtain MPS, vertical stress and UYS
            sigmav_init = sigmav[i*N-i]
            [sigmav_o, sigma1_o, sigmaf_o, UYSf_o, RH_diameter]=Janssen_Equation(KK,X1,X2,Z1,Z2,N,sigmav_init)
            sigmav[i*N:(i+1)*N] = sigmav_o[:N]
            sigma1_active[i*N:(i+1)*N] = sigma1_o[:N]                       # MPS
            sigma1_passive[i*N:(i+1)*N] = sigma1_o[:N]                      # MPS for the passive mode which is equal to active mode for the Janssen equation
            sigmaf[i * N:(i + 1) * N] = sigmaf_o[:N]                        # sigmaf is useful only for the funnel flow case to evaluate the possibility of ratholing (Eq. (22) of the reference)
            UYS_active[i * N:(i + 1) * N] = UYSf_o[:N]                      # UYS for active mode
            UYS_passive[i * N:(i + 1) * N] = UYSf_o[:N]                     # UYS for the passive mode which is equal to active mode for the Janssen equation

        else:                                                               # otherwise we are in the conical part of the hopper

            ## We use Motzkus equation to obtain vertical stress and MPS in the active mode
            sigmav_init = sigmav[i*N-i]
            [sigmav_o, sigma1_o, UYSf_o, sigmaf_o, RH_diameter]=properties_in_the_conical_part_in_active_mode(X1,X2,Z1,Z2,N,sigmav_init)
            sigmav[i*N:(i+1)*N] = sigmav_o[:N]                                                                          #vertical stress in the active mode (Pa) in the conical part of the hopper
            sigma1_active[i*N:(i+1)*N] = sigma1_o[:N]                                                                   #MPS in the active mode (Pa) in the conical part of the hopper
            UYS_active[i * N:(i + 1) * N] = UYSf_o[:N]                                                                  #unconfined yield strength in the active mode (Pa)
            sigmaf[i * N:(i + 1) * N]=sigmaf_o[:N]                                                                      #sigmaf is useful only for the funnel flow case to evaluate the possibility of ratholing (Eq. (22) of the reference)

            ## We use radial stress theory to estimate major principal stress in the passive mode
            sigmav_init = sigmav[i*N-i]
            [sigma1_o, UYSf_o]=properties_in_the_conical_part_in_passive_mode(X1,X2,Z1,Z2,N,sigmav_init)
            sigma1_passive[i*N:(i+1)*N] = sigma1_o[:N]                                                                  #MPS in the passive mode (Pa) in the conical part of the hopper
            UYS_passive[i * N:(i + 1) * N] = UYSf_o[:N]                                                                 #unconfined yield strength in the passive mode (Pa)

    [Q, M, F, P, theta, theta_critical] = MassFlow_or_FunnelFlow(X1, X2, Z1, Z2, sigma1_active, sigma1_passive, sigmaf, N, number)

    return Q, M, F, P, theta, theta_critical, z, sigmav, sigma1_active, sigma1_passive, UYS_active, UYS_passive, sigmaf, RH_diameter, theta


