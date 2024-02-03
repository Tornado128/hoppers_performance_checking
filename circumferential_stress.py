def circumferential_stress(N, number, sigma1_active, X2):
    import numpy as np
    from curve_fitting import curve_fitting

    sigmaf_o = np.zeros(N*number)
    rhob = np.zeros(N*number)

    B = 2 * X2                                                                                                          # diameter of the outlet of the hopper (m)
    g = 9.8
    # This function fits bulk density, effective angle of internal
    # friction, FC and FFC vs sigma1. It also does a power fit for WFA vs normal stress
    # a, b and c are the coefficients for power and linear curve fittings.
    [a, b, c, min_rhob, max_PHIE, max_WFA] = curve_fitting()

    for j in range(number*N):

        PHILIN_p = a[3] * sigma1_active[j] + b[3]                                                                       # Eq. (23) of the eference
        G = -6.86712 + 0.58911*PHILIN_p-0.012966*PHILIN_p**2.0+0.00011939*PHILIN_p**3.0                                 # Jenike Bulletin 123 P67

        rhob[j] = a[0] * sigma1_active[j] + b[0]
        sigmaf_o[j] = rhob[j] * g * B / G

    # circumferential stress
    return sigmaf_o
