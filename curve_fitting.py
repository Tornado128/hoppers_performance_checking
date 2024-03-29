## This function reads stress walk data (FFC as a function of normal stress) and wall friction \
## data (shear stress as a function of normal stress)
## This function also perform curve fitting

def curve_fitting():
	import numpy as np
	from scipy.optimize import curve_fit
	from input_files import input_files

	# This function returns the stress walk data and wall friction data obtained from ring shear cell tester.
	[df_FFC, df_wall] = input_files()

	sigma1 = df_FFC.iloc[:,2]               # sigma1 is the major principal stress (MPa)
	FC = df_FFC.iloc[:,3]                   # FC is the cohesive strength (MPa)
	PHILIN = df_FFC.iloc[:,9]               # linearized angle of internal friction (see: https://www.mdpi.com/2075-4701/8/4/255)
	rhob = df_FFC.iloc[:,7]                 # bulk density (kg/m3)
	PHIE = df_FFC.iloc[:,8]                 # effective angle of internal friction (see: https://www.researchgate.net/figure/Illustration-of-a-powder-yield-locus-j-kinetic-angle-of-internal-friction-t_fig1_26551784)
	SIGMA = df_wall.iloc[:,1]				# wall normal stress in the wall friction data (pa)
	TAU = df_wall.iloc[:,2]					# shear stress in the wall friction data (pa)

	WFA = np.arctan(TAU/SIGMA)*180/np.pi	# wall friction angle (degree)

	## Average wall friction angle and average angle of effective wall friction angle
	## we need maximum wall friction angle, maximum effective angle of friction and minimum density for passive state.
	max_WFA = np.max(WFA)		#51
	max_PHIE = np.max(PHIE)		#62.85
	min_rhob = np.min(rhob)		#264.1

	## defining arrays of zero. a, b and c are the fitting coefficients for linear and powder curve fittings.
	## The goal of fitting is to have a continuous function for FC, FFC, rhob, PHIE and PHILIN as a function of sigma1 (major principal stress)
	a = np.zeros(6)
	b = np.zeros(6)
	c = np.zeros(6)

	## Linear objective function
	def objective_linear(x, a, b):
		return a * x + b

	## Power objective function
	def objective_2nd_oder_polynomial(x, a, b, c):
		return a * x ** 2 + b * x + c

	## Power objective function for WFA (the intercept is defined as 90 degrees because we expect that wall friction angle approaches 90 degrees at low sigma)
	def objective_power_wall(x, a, b):
		return a * x**b + 90
	#

	popt, _ = curve_fit(objective_linear, sigma1, rhob)     			## curve fitting for rhob vs sigma1 (power equation)
	a[0], b[0] = popt                                             		## summarize the parameter values
	popt, _ = curve_fit(objective_linear, sigma1, PHIE, maxfev=5000)    ## curve fitting for PHIE vs sigma1 (power equation)
	a[1], b[1] = popt                                             		## summarize the parameter values
	popt, _ = curve_fit(objective_2nd_oder_polynomial, sigma1, FC)      ## curve fitting for FC vs sigma1 (linear equation)
	a[2], b[2], c[2] = popt                                             ## summarize the parameter values
	popt, _ = curve_fit(objective_linear, sigma1, PHILIN)               ## curve fitting for PHILIN vs sigma1 (linear equation)
	a[3], b[3] = popt                                                   ## summarize the parameter values
	popt, _ = curve_fit(objective_power_wall, SIGMA, WFA, maxfev=5000)  ## curve fitting for WFA vs SIGMA (power equation)
	a[4], b[4] = popt                                                   ## summarize the parameter values
	c[4] = 90
	popt, _ = curve_fit(objective_linear, SIGMA, TAU)  					## curve fitting for shear stress vs normal stress in the wall friction data
	a[5], b[5] = popt

	return (a, b, c, min_rhob, max_PHIE, max_WFA)													## returning a, b and c arrays, which are the coefficients for curve fitting. Remember c is zero for some!

'''
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
y_pred_rhob = objective_power(sigma1, a[0], b[0], c[0])
print('y = %.8f * x^ %.8f + %.8f' % (a[0], b[0], c[0]), '; R2 = ', r2_score(y_pred_rhob, rhob))

y_pred_PHIE = objective_power(sigma1, a[1], b[1], c[1])
print('y = %.8f * x^ %.8f + %.8f' % (a[1], b[1], c[1]), '; R2 = ', r2_score(y_pred_PHIE, PHIE))

y_pred_FC = objective_linear(sigma1, a[2], b[2])
print('y = %.8f * x + %.8f' % (a[2], b[2]), '; R2 = ', r2_score(y_pred_FC, FC))

y_pred_PHILIN = objective_linear(sigma1, a[3], b[3])
print('y = %.8f * x + %.8f' % (a[3], b[3]), '; R2 = ', r2_score(y_pred_PHILIN, PHILIN))

y_pred_WFA = objective_power_wall(SIGMA, a[4], b[4])
print('y = %.8f * x^ %.8f + 90' % (a[4], b[4]), '; R2 = ', r2_score(y_pred_WFA, WFA))

plt.plot(sigma1, PHIE, 'o',sigma1, y_pred_PHIE, '--')
plt.xlabel("major principal stress (pa)",fontsize=16)
plt.ylabel("effective angle of internal friction (degree)",fontsize=16)
plt.xlim((0, max(sigma1)))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

plt.plot(sigma1, rhob, 'o',sigma1, y_pred_rhob, '--')
plt.xlabel("major principal stress (pa)",fontsize=16)
plt.ylabel("bulk density (g/ml)",fontsize=16)
plt.xlim((0, max(sigma1)))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

plt.plot(sigma1, FC, 'o',sigma1, y_pred_FC, '--')
plt.xlabel("major principal stress (pa)",fontsize=16)
plt.ylabel("cohesive strength (Pa)",fontsize=16)
plt.xlim((0, max(sigma1)))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

plt.plot(sigma1, PHILIN, 'o',sigma1, y_pred_PHILIN, '--')
plt.xlabel("major principal stress (pa)",fontsize=16)
plt.ylabel("linearized angle of internal friction (degree)",fontsize=16)
plt.xlim((0, max(sigma1)))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

plt.plot(SIGMA, WFA, 'o',SIGMA, y_pred_WFA, '*')
plt.xlabel("consolidation pressure (pa)",fontsize=16)
plt.ylabel("wall friction angle",fontsize=16)
plt.xlim((0, max(SIGMA)))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
##plt.show()
'''