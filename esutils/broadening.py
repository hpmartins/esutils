import numpy as np
from scipy.special import wofz

################
################
## BROADENING ##
################
################

def voigt(x, y, new_x, psigma, pgamma, sticks=False, aL=0.0):
    if np.all(y == 0):
        return np.zeros_like(new_x)

    out = np.zeros_like(new_x)
    if sticks:
        old = np.sum(y)
    else:
        old = np.trapz(y, x)
    if np.size(psigma) == 1:
        psigma = psigma*np.ones_like(y)
    if np.size(pgamma) == 1:
        pgamma = pgamma*np.ones_like(y)
    if aL > 0:
        pgamma = pgamma*(1 + aL*(x - x[0]))
    for i in np.arange(np.size(y)):
        out += voigt_prof(new_x, amplitude=y[i], center=x[i], sigma=psigma[i], gamma=pgamma[i])
    return old*out/np.trapz(out, new_x)

def gaussian(x, y, new_x, sigma, sticks=False):
    if np.all(y == 0):
        return np.zeros_like(new_x)

    out = np.zeros_like(new_x)
    if sticks:
        old = np.sum(y)
    else:
        old = np.trapz(y, x)
    if np.size(sigma) == 1:
        sigma = sigma*np.ones_like(y)
    for i in np.arange(np.size(y)):
        out += gaussian_prof(new_x, amplitude=y[i], center=x[i], sigma=sigma[i])
    return old*out/np.trapz(out, new_x)

def lorentzian(x, y, new_x, gamma, sticks=False, aL=0.0):
    out = np.zeros_like(new_x)
    if sticks:
        old = np.sum(y)
    else:
        old = np.trapz(y, x)
    if np.size(gamma) == 1:
        gamma = gamma*np.ones_like(y)
    if aL > 0:
        gamma = [gamma*(1 + aL*(a - x[0])) for a in x]
    for i in np.arange(np.size(y)):
        out += lorentzian_prof(new_x, amplitude=y[i], center=x[i], gamma=gamma[i])
    return old*out/np.trapz(out, new_x)

def doniachsunjic(x, y, new_x, sigma_d, gamma):
    out = np.zeros_like(new_x)
    for i in np.arange(np.size(y)):
        out += doniachsunjic_prof(new_x, y[i], x[i], sigma_d, gamma)
    return out

def gaussian_prof(x, amplitude=1.0, center=0.0, sigma=1.0):
    return (amplitude/(np.sqrt(2*np.pi)*sigma)) * np.exp(-(1.0*x-center)**2 / (2*sigma**2))

def lorentzian_prof(x, amplitude=1.0, center=0.0, gamma=1.0):
    return (amplitude/(1 + ((1.0*x-center)/gamma)**2)) / (np.pi*gamma)

def doniachsunjic_prof(x, amplitude=1.0, center=0.0, sigma=1.0, gamma=0.0):
    arg = (x-center)/sigma
    gm1 = (1.0 - gamma)
    scale = amplitude/(sigma**gm1)
    return scale*np.cos(np.pi*gamma/2 + gm1*np.arctan(-arg))/(1 + arg**2)**(gm1/2)

def voigt_prof(x, amplitude=1.0, center=0.0, sigma=0.5, gamma=0.5):
    alpha = sigma/np.sqrt(2*np.log(2))
    z = ((x - center) + 1j*gamma) / (alpha*np.sqrt(2))
    return amplitude*np.real(wofz(z))/(alpha*np.sqrt(2*np.pi))
