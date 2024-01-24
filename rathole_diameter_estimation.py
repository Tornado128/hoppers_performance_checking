def rathole_diameter_estimation(N, number, sigma1_active, X2):
    import numpy as np
    from curve_fitting import curve_fitting
    import math

    sigmaf_o = np.zeros(N*number)
    rhob = np.zeros(N*number)
    UYS = np.zeros(N*number)

    B = 2 * X2                                                                                                          # diameter of the outlet of the hopper (m)
    g = 9.8
    # This function fits bulk density, effective angle of internal
    # friction, FC and FFC vs sigma1. It also does a power fit for WFA vs normal stress
    # a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()

    # stress in the abutment (only for the case of funnel flow: Eq. (23) of Leung et al, J. Pharmaceuticaals Sciences, 2019, 108, 264)
    for i in range(number*N):

        PHILIN_p = a[3] * sigma1_active[i] + b[3]                                                                       # Eq. (23) of the eference
        G = -6.86712 + 0.58911*PHILIN_p-0.012966*PHILIN_p**2.0+0.00011939*PHILIN_p**3.0                                 # Jenike Bulletin 123 P67

        rhob[i] = a[0] * sigma1_active[i] + b[0]
        sigmaf_o[i] = rhob[i] * g * B / G
        UYS[i] = a[2]*sigma1_active[i]**2 + b[2]*sigma1_active[i] + c[2]

    RH_diameter = G * UYS[-1] / (rhob[-1] * g)                                                                          # Eq. (22) of the reference: Rathole diameter (m)

    return RH_diameter, sigmaf_o
