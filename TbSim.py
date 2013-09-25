# Tool for plotting various tau-omega model curves

from __future__ import division
from numpy import *
from matplotlib.pyplot import *

from n_soil import *
from trb_lossy import *

# -------------- inputs

print('\nPlease input parameter values, empty for [default values]:')
Ts  = float(raw_input('Temperature, K              [300] --> ') or 300) # temperature [K]
mv  = float(raw_input('Vol. Soil Moisture, m3/m3  [0.35] --> ') or 0.35) # soil moisture [m3/m3]
tau = float(raw_input('Optical Thickness           [0.1] --> ') or 0.1) # optical thickness
hp  = float(raw_input('Roughness Parameter         [0.1] --> ') or 0.1) # roughness parameter

# -------------- constants

freq = 1.4e9 # frequency, Hz

theta = linspace(0,60,61) # incidence angles, deg
mu = cos(deg2rad(theta))  # direction cosine

er_a = 1 # relative permittivity of air

# IVS site 703 0-20 cm particle size analysis: 38% sand, 33% silt, 29% clay
# IVS site 709 0-20 cm particle size analysis: 36% sand, 35% silt, 28% clay.
clay = 30 # soil clay content, %

# --------------- computations

# soil electrical properties (from Mironov model)
n = n_soil(freq,clay,mv)
erprime = (real(n))**2 - (abs(imag(n)))**2
erdoubleprime = 2*real(n)*abs(imag(n))
er_s = erprime + 1j*erdoubleprime

[dum,Rs_h,duh] = trb_lossy(mu,er_a,er_s,2,freq)  # smooth soil h-pol reflectivity
[dum,Rs_v,duh] = trb_lossy(mu,er_a,er_s,1,freq)  # smooth soil v-pol reflectivity

Rr_h = Rs_h*exp(-hp*mu**2)  # rough soil h-pol reflectivity
Rr_v = Rs_v*exp(-hp*mu**0)  # rough soil v-pol reflectivity
    
# use simplified tau-omega equation (omega = 0)
Tb_h = Ts*(1 - Rr_h * exp(-2*tau/mu))
Tb_v = Ts*(1 - Rr_v * exp(-2*tau/mu))
        
# --------------- plotting
figure(1)
plot(theta, Tb_v, 'b-', lw=2); hold(True)
plot(theta, Tb_h, 'b--', lw=2); grid(True)
xlabel('Incidence Angle [$^{\\circ}$]')
ylabel('Brightness Temperature [K]')
legend(['v-pol', 'h-pol'], loc=2)
title('$T =$ %d K, $\\theta_v =$ %.2f m$^3$ m$^{-3}$, $\\tau =$ %.2f, $h =$ %.2f' % (Ts, mv, tau, hp)) 
tight_layout()
show()
