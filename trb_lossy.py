from numpy import *

# Brian Hornbuckle, 3-07-98.  Modified November 27, 2001.  Modified May 25, 2010.
# Converted to Python by Jason Patton, January 10, 2012

def trb_lossy(mu1, er1, er2, p, f):

  """
Given an direction cosine incident on (lossy) medium 2 from medium 1,
this function returns the transmission direction cosine
and the reflectivity and transmissivity for brightness temperatures.
These expressions happen to be the same as the power transmissivity and reflectivity.

Parameters
----------
mu1 : array_like
    Incident direction cosine, i.e. cos(theta_incident)
er1 : scalar
    Relative dielectric constant of medium 1
er2 : scalar
    Relative dielectric constant of medium 2
p : int
    Polarization:
      1 = vertical
      2 = horizontal
f : scalar
    Frequency [Hz]

Returns
-------
mu2 : ndarray
    Transmission direction cosine, i.e. cos(theta_transmitted)
R : ndarray
    Reflectivity for brightness temperature
T : ndarray
    Transmissivity for brightness temperature

Notes
_____
Polarization orientation follows FSA convention, i.e. h=zxk, v=hxk, where
k is the direction of propagation.
Only the real part of medium 1's relative dielectric constant is considered,
so medium 1 must be low loss for good results.  Medium 2 can be lossy,
in which case the real angle of transmission is used to calculate the
outputs.  The real angle of transmission approach is used when
when tan(delta2) > 0.1.
"""

  er1 = real(er1)
  n1 = sqrt(er1)         # index of refraction in medium 1
  er2p = real(er2)       # er2': relative permittivity
  er2dp = abs(imag(er2)) # er2'': loss

#------------------------- medium 2 is lossy --------------------------------

  n2 = sqrt(er2)                      # index of refraction in 2
  ko = 2*pi*f*sqrt(pi*4e-7*8.854e-12) # free space wave number
  alpha = ko*abs(imag(sqrt(er2)))     # attenuation constant
  beta = ko*real(sqrt(er2))           # phase constant

  # chi = real angle of transmission
  s = 2*alpha*beta
  q = (beta**2)-(alpha**2)-(ko**2)*(n1**2)*(1-mu1**2)
  chi = arctan(sqrt(2)*sqrt((ko**2)*(n1**2)*(1-mu1**2))/sqrt(sqrt((s**2)+(q**2))+q))
  mu2 = cos(chi)

  # Fresnel coefficients
  if(p == 2): # horizontal polarization
    gamma = (n1*mu1-n2*mu2)/(n1*mu1+n2*mu2)
    tau = 1+gamma
  elif(p == 1):	# vertical polarization
    gamma = (n2*mu1-n1*mu2)/(n2*mu1+n1*mu2)
    tau = mu1/mu2*(1-gamma)
  else:
    error('v or h not entered for polarization')

  R = abs(gamma)**2
  T = 1-R	

  return [mu2, R, T]
