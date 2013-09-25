from numpy import *

from refractive_index import *

# function to computer soil dielectric constant using Mironov et al 2009
# Brian Hornbuckle, May 25, 2010.
# Converted to Python by Jason Patton, January 10, 2012
# uses function "refractive_index.m"
# output:  n = n' - j*n" = moist soil index of refraction, ejwt time dependence
# inputs:
# f = frequency, Hz
# c = clay content, %
# mv = volumetric water content, m^3 m^-3

def n_soil(f,c,mv):

# ------ questions
# eou is set to 100, 
# but Debye water model has eou = 88.045 - Tc.*0.4147 + Tc**2.*6.295e-4 + Tc**3.*1.075e-5, where Tc = Celsius temperature
# and this yields eou = 76 to 86 for Tc = 5 to 30 deg C
# relaxation time tu agrees with water model

# ------- notes
# temperature dependence was not tested, 
# so best if used near 293 to 295 K, 
# the temperature at which data was collected

    eo = 8.854e-12 # permittivity of free space, F m^-1

# ----- empirical models

    nd = 1.634 - 0.539e-2*c + 0.2748e-4*c**2 # Re(dry soil refractive index) = n'_dry
    kd = 0.03952 - 0.04038e-2*c # Im(dry soil refractive index) = n"_dry
    mvt = 0.02863 + 0.30673e-2*c # max bound water volumetric water content, m^3 m^-3 
    eob = 79.8 - 85.4e-2*c + 32.7e-4*c**2 # low frequency limit of bound water dielectric constant
    tb = 1.062e-11 + 3.450e-12*c*1e-2 # relaxation time for bound water
    sb = 0.3112 + 0.467e-2*c # bound water conductivity
    su = 0.3631 + 1.217e-2*c # free (unbound) water conductivity
    eou = 100 # low frequency limit of free (unbound) water dielectric constant
    tu = 8.5e-12 # relaxation time for free (unbound) water, s

# -------- bound water Debye model

    einf = 4.9 # high-frequency (optical) limit of relative permittivity

    erprime = einf + (eob-einf)/(1+((2*pi)*f*tb)**2) # real part of relative permittivity

    erdoubleprime = (eob-einf)/(1+((2*pi)*f*tb)**2)*(2*pi)*f*tb + sb/((2*pi*eo)*f) # imaginary part of relative permittivity

    erb = erprime + 1j*erdoubleprime # bound water relative permittivity or dielectric constant

    bound = refractive_index(erb)
    nb = real(bound)
    kb = abs(imag(bound)) # bound water index of refraction

# -------- free (unbound) water Debye model

    erprime = einf + (eou-einf)/(1+((2*pi)*f*tu)**2) # real part of relative permittivity

    erdoubleprime = (eou-einf)/(1+((2*pi)*f*tu)**2)*(2*pi)*f*tu + su/((2*pi*eo)*f) # imaginary part of relative permittivity

    eru = erprime + 1j*erdoubleprime # bound water relative permittivity or dielectric constant

    free = refractive_index(eru)
    nu = real(free)
    ku = abs(imag(free)) # free (unbound) water index of refraction

# ------- moist soil index of refraction

    if(mv < mvt): # all soil water is bound
        nm = nd + (nb-1)*mv # real part of refractive index of moist soil
        km = kd + kb*mv # imaginary part of moist soil refractive index
    else: # some free (unbound) water exists
        nm = nd + (nb-1)*mvt + (nu-1)*(mv-mvt) # real part of refractive index of moist soil
        km = kd + kb*mvt + ku*(mv-mvt) # imaginary part of moist soil refractive index

    return nm - 1j*km #  moist soil refractive index, ejwt time dependence
