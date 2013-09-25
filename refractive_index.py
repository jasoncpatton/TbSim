from numpy import *

# Brian K. Hornbuckle 2-06-2006
# Converted to Python by Jason Patton, January 10, 2012
# This matlab m-file calculates the index of refraction
# when given the relative permittivity.
#
# Input:   e_r = relative permittivity (complex in general) (can be a vector)
# Output:  n = index of refraction (complex in general) (will be a vector if e_r is a vector)

def refractive_index(e_r):

    e_r_prime = abs(real(e_r))       # real part of the relative permittivity
    e_r_doubleprime = abs(imag(e_r)) # imaginary part

    # real part of refractive index
    n_prime = sqrt((sqrt(e_r_prime**2 + e_r_doubleprime**2) + e_r_prime)/2)
    # imaginary part
    n_doubleprime = sqrt((sqrt(e_r_prime**2 + e_r_doubleprime**2) - e_r_prime)/2)

    return n_prime - 1j*n_doubleprime # return the index of refraction
