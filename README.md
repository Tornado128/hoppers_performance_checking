(1) These codes here are translated into Python as a personal side project and just as a demo of a more advanced version, which was developed in MATLAB. The input and dimensions provided in this code are not real purposely.
---
(2) The code is mostly based on the methodology developed by Leung et al, "A proposed complete methodology to predict gravity flow obstruction of pharmaceutical 
powders in drug product manufacturing." Journal of pharmaceutical sciences 108.1 (2019): 464-475. However, some of the equations in "Dietmar Schulze, 
The prediction of initial stresses in hoppers, Bulk Solids Handling, 1994, 14, 497-503" are also used. 
-----
(3) The vertical stress profile predicted by the code matches with the predictions by Dietmar online software: https://dietmar-schulze.de/sstool_e.html
-----
(4) Remember that Jenike method is very conservative for the estimation of arch diameter. It is even more conservative
for the prediction of rathole diameter! As a rule of thumb, if Jenike method says that a system doesn't rathole, it will definitely NOT rathole! 
-----
(5) If the hopper outlet is not fully activated (as in most tablet presses), Jenike method is even not directionally correct!! The effective angle of internal friction is more predictive of ratholing in such scenarios, based on my experiences. Some expert claims that average angle of friction higher than 46 degree cause ratholing.
-----
(6) The purpose of this code is to evaluate the risks involved in gravity-driven flow of pharmaceutical powders in EXISTING hoppers, including mass flow/funnel 
pattern, arch formation under active stress (filling, initial discharging) and passive state (emptying) and rathole formation.
-----
(7) In the cylinderical part of the hopper, the stress (vertical stress and MPS) profiles are obtained using Janssen equation. The stress profile is the active state in the conical part is obtained using Motzkus equation
-----
