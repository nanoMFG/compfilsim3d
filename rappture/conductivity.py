import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from numpy import sqrt, sin, cos, pi, exp, log

# cond_g is the conductivity in presence of both GBS and BS
# lambda_0 is background mean free path
# a is nanowire radius
# p is proportion of elastically scattered electrons
# v is the velocity vector magnitude
# r is the radius vector
# alpha is the 
# d = D is the avg grain diameter
# R is the reflection probability of an electron at a single plane.

def total_conductivity(a, lambda_0, R, d, theta, cond_0, p):
    # C = D/(2*a) = 1
    k = 2*a/lambda_0
    t = 1/sin(theta)
    alpha = (lambda_0/d)*(R/(1-R))
    cond_g = 3*cond_0*(1/3 - .5*alpha + alpha**2 - alpha**3*log(1 + 1/alpha))
    H = 1 + alpha/(1-1/t**2)**.5
    f1 = lambda t : (-1/(H)**2)*(t**-5)*sqrt(t**2 - 1)
    integral_1 = integrate.quad(f1, 1, np.inf) 
    # Approximate the geometric series
    def sum_of_series(max):
        total = 0
        for v in range(max):
            f2 = lambda psi : (exp(-(v+1)*H*k*t*sin(psi)) - exp(-v*H*k*t*sin(psi)))*sin(psi)
            integral_2 = integrate.quad(f2, 0, pi/2)
            total = total + p**v * integral_2[0]
        return total
    total_cond = cond_g - (12/(pi*k))*(1-p)*integral_1[0]*sum_of_series(10)
    return total_cond