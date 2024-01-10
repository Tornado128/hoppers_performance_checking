
def MassFlow_or_FunnelFlow(X1, X2, Z1, Z2, sigma1_active, sigma1_passive, sigmaf, N, number, RADIUS):

    import numpy as np
    from curve_fitting import curve_fitting

    [a, b, c, average_rhob, average_PHIE, average_WFA] = curve_fitting()

    theta = (np.pi / 2 - np.arctan((Z2 - Z1) / (X2 - X1+ 0.0000000001))) * 180 / np.pi                                  # hopper angle form vertical (degree)
    rhob_active_outlet = a[0] * sigma1_active[-1] + b[0]                                                                # bulk density at the outlet of the hopper in the active mode (Pa)
    UYS_active_outlet = a[2] * sigma1_active[-1] + b[2]                                                                 # UYS at the outlet of the hopper in the active mode(Pa)

    H = (130+theta)/65                                                                                                  # Eq. (3) of Leung et al
    sigma_active_outlet = rhob_active_outlet*9.8* 2 * X2 /H                                                         # Eq. (2) of Leung et al to estimate the stress in the outlet (Pa)

    ## ACTIVE STATE
    ## If the stress on the abutment (caused by gravity) is more than UYS, mass flow occurs in the active state
    if (sigma_active_outlet > UYS_active_outlet):                                                                       # see Eq. (1) in the reference
        M = 0  # no arch in the active state
    else:
        M = 1  # arch formation in the active state

    ## Passive State
    PHIE_passive_outlet = a[1] * sigma1_passive[-1] + b[1]                                                              # bulk density at the outlet of the hopper in the active mode (Pa)
    WFA_passive_outlet = a[4] * sigma1_passive[-1] ** b[4] + c[4]                                                        # wall friction angle at the outlet of the hopper in the active mode(Pa)

    # Eq. (5) and Eq. (6) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
    if (WFA_passive_outlet>PHIE_passive_outlet):
         PHIE_passive_outlet = WFA_passive_outlet                                                                         # To avoid nan in the estimation of beta in the next line
    beta = (WFA_passive_outlet + (180 / np.pi) * np.arcsin(np.sin(WFA_passive_outlet * np.pi / 180)
                                                           / np.sin(PHIE_passive_outlet * np.pi / 180))) / 2
    theta_critical = 90 - 0.5 * (180 / np.pi) * np.arccos((1 - np.sin(PHIE_passive_outlet * np.pi / 180))
                                                          / (2 * np.sin(PHIE_passive_outlet * np.pi / 180))) - beta

    ## angle of the outlet of the hopper from a vertical line
    # (I added 0.000000000000001 to the denominator to avoid having the denominator equal to zero)
    theta = 90 - abs(np.arctan((Z2 - Z1) / (X2 - X1 + 0.000000000000001))) * 180 / np.pi

    F = -1                                          # Initialization of F: Later, F=1 is funnel flow and F=0 is mass flow in the passive state
    P = -1                                          # Initialization of P: Later, If F=0, P=0 means no arch forms while P=1 is equivalent to arch formation
    Q = -1                                          # Initialization of P: Later, if there is a funnel flow (F=1), the value of Q changes into 2 in the case of arch formation. It changes to -2 if no arch forms
    if (theta>theta_critical):
        F = 1                                       # funnel flow in the passive state
    else:
        F = 0                                       # mass flow in the passive state

    if (F == 0):                                    # if there is a mass flow, determine if there is a risk of arch formation

        UYS_passive_outlet = a[2]*sigma1_passive[-1]+b[2]
        if (sigma1_passive[-1]>UYS_passive_outlet):
            P = 0                           # no arch in the passive state
        else:
            P = 1                           # arch formation in the passive state

    ## If there is a funnel flow, determine if there is a risk of rathole formation in the ENTIRE depth of the powder bed
    # in the hopper
    a2 = a[2]
    b2 = b[2]
    UYS_active = np.zeros(number*N)
    if (F == 1):
        for i in range(number*N):
            UYS_active[i] = a2 * sigma1_active[i] + b2
            if (sigmaf[i]<UYS_active[i]):
                P = 2                   #rathole formation
            else:
                P = -2                  #no rathole formation

    ## If there is a funnel flow, determine if there is a risk of arch formation
    # in the hopper
    a0 = a[0]
    b0 = b[0]
    if (F == 1):
        UYS_passive_outlet = a2 * sigma1_passive[-1] + b2
        rhob_passive_outlet = a0 * sigma1_passive[-1] + b0

        ## Eq. (2) and Eq. (3)
        H = (130 + theta)/65
        sigma = rhob_passive_outlet*9.8*X2*2/H

        if (sigma<UYS_passive_outlet):
            Q = 2                   #arch formation
        else:
            Q = -2                  # no arch formation
    return Q, M, F, P, theta, theta_critical