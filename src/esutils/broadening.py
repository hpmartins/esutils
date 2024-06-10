import numpy as np
from scipy.special import wofz

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
        out += voigt_k(new_x, I=y[i], x0=x[i], sigma=psigma[i], gamma=pgamma[i])
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
        out += gaussian_k(new_x, I=y[i], x0=x[i], sigma=sigma[i])
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
        out += lorentzian_k(new_x, I=y[i], x0=x[i], gamma=gamma[i])
    return old*out/np.trapz(out, new_x)

def doniachsunjic(x, y, new_x, sigma_d, gamma):
    out = np.zeros_like(new_x)
    for i in np.arange(np.size(y)):
        out += doniachsunjic_k(new_x, y[i], x[i], sigma_d, gamma)
    return out

def gaussian_k(x, I=1.0, x0=0.0, sigma=1.0):
    return (I/(np.sqrt(2*np.pi)*sigma)) * np.exp(-(1.0*x-x0)**2 / (2*sigma**2))

def lorentzian_k(x, I=1.0, x0=0.0, gamma=1.0):
    return (I/(1 + ((1.0*x-x0)/gamma)**2)) / (np.pi*gamma)

def doniachsunjic_k(x, I=1.0, x0=0.0, sigma=1.0, gamma=0.0):
    arg = (x-x0)/sigma
    gm1 = (1.0 - gamma)
    I0 = I/(sigma**gm1)
    return I0*np.cos(np.pi*gamma/2 + gm1*np.arctan(-arg))/(1 + arg**2)**(gm1/2)

def voigt_k(x, I=1.0, x0=0.0, sigma=0.5, gamma=0.5):
    alpha = sigma/np.sqrt(2*np.log(2))
    z = ((x - x0) + 1j*gamma) / (alpha*np.sqrt(2))
    return I*np.real(wofz(z))/(alpha*np.sqrt(2*np.pi))
