
def MassFlow_or_FunnelFlow( X2, sigma1_active, sigmaf, N, number, ANGLE, UYS_passive_conical_outlet, sigma1_passive_conical_outlet, WFA_passive_conical_outlet,  PHIE_passive_conical_outlet,UYS_active_conical_outlet, rhob_active_conical_outlet, RH_diameter):

    import numpy as np
    from curve_fitting import curve_fitting

    [a, b, c, min_rhob, max_PHIE, max_WFA] = curve_fitting()

    H = (130+ANGLE)/65                                                                                                  # Eq. (3) of Leung et al
    sigma_active_outlet = rhob_active_conical_outlet * 9.8 * 2 * X2 /H                                                  # Eq. (2) of Leung et al to estimate the stress in the outlet (Pa)

    ## ACTIVE STATE
    ## If the stress on the abutment (caused by gravity) is more than UYS, mass flow occurs in the active state
    if (sigma_active_outlet > UYS_active_conical_outlet):                                                                       # see Eq. (1) in the reference
        M = 0                       # no arch in the active state
    else:
        M = 1                       # arch formation in the active state

    ## Passive State
    # Eq. (5) and Eq. (6) in Leung et al, J. pharmaceutical sciences, 108 (2019) 464-475
    if (WFA_passive_conical_outlet>PHIE_passive_conical_outlet):
         PHIE_passive_conical_outlet = WFA_passive_conical_outlet                                                                         # To avoid nan in the estimation of beta in the next line
    beta = (WFA_passive_conical_outlet + (180 / np.pi) * np.arcsin(np.sin(WFA_passive_conical_outlet * np.pi / 180)
                                                           / np.sin(PHIE_passive_conical_outlet * np.pi / 180))) / 2
    theta_critical = 90 - 0.5 * (180 / np.pi) * np.arccos((1 - np.sin(PHIE_passive_conical_outlet * np.pi / 180))
                                                          / (2 * np.sin(PHIE_passive_conical_outlet * np.pi / 180))) - beta


    if (ANGLE>theta_critical):
        F = 1                                       # funnel flow in the passive state
    else:
        F = 0                                       # mass flow in the passive state

    if (F == 0):                                    # if there is a mass flow in the passive state, determine if there is a risk of arch formation

        #UYS_passive_outlet = a[2]*sigma1_passive[-1]**2 +b[2]*sigma1_passive[-1] + c[2]
        if (sigma1_passive_conical_outlet>UYS_passive_conical_outlet):
            P = 0                           # no arch in the passive state
        else:
            P = 1                           # arch formation in the passive state

    ## If there is a funnel flow, determine if there is a risk of rathole formation in the ENTIRE depth of the powder bed
    ## in the hopper
    #a2 = a[2]
    #b2 = b[2]
    #c2 = c[2]
    #UYS_active = np.zeros(number*N)
    #if (F == 1):                ## if it is a funnel flow in passive state
    #    for i in range(number*N):
    #        UYS_active[i] = a2 * sigma1_active[i]**2 + b2*sigma1_active[i] + c2
    #        if (sigmaf[i]<UYS_active[i]):
    #            P = 2                   #rathole formation
    #        else:
    #            P = -2                  #no rathole formation

    if ( F== 1):                                    ## if it is a funnel flow in passive state
        if RH_diameter > (2*X2):                    ## if the rathole diameter is larger than the outlet diameter
            P = 2                                   #rathole formation is predicted
        else:
            P = -2                                  #no rathole formation is predicted

    return M, F, P, ANGLE, theta_critical