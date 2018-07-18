import numpy as np
import pandas as pd
from scipy.special import wofz
from scipy.signal import argrelmax

##########################
##########################
######### QUANTY #########
##########################
##########################

def normalize(A, nmin=True, norm=1.0):
    vmin = 0
    if nmin:
        vmin = np.amin(A)
    return norm*(A-vmin)/(np.max(A)-vmin)

def expbg(x, y, new_x):
    y = normalize(y[np.argsort(x)])
    x = x[np.argsort(x)]
    return np.interp(new_x, x, y[-1]*normalize(np.cumsum(y)))

def importexp(fn, norm=True, nmin=True, delimiter=' ', col=1):
    t = np.loadtxt(fn, delimiter=delimiter, usecols=(0,col))
    E = t[:, 0]
    I = t[:, 1:]

    I = I[np.argsort(E)]
    E = E[np.argsort(E)]
    if norm:
        for i in np.arange(np.size(I, axis=1)):
            I[:, i] = normalize(I[:, i], nmin=nmin)
    return E, I

def calc_xasmag(clcf, xr, shft=0.0, norm=True):
    # Data import
    E, Gz, Gr, Gl = np.loadtxt(clcf, unpack=True, skiprows=5, usecols=(0,2,4,6))

    # Spectral weight
    Gz = -(1/np.pi)*Gz
    Gr = -(1/np.pi)*Gr
    Gl = -(1/np.pi)*Gl
    Gi = Gz + Gl + Gr

    # G(w) and dichroism
    G = np.interp(xr, E - E[argrelmax(Gi)[0][0]] + shft, Gi)
    D = np.interp(xr, E - E[argrelmax(Gi)[0][0]] + shft, Gr - Gl)

    if norm:
        G = normalize(G)

    return G, D

def calc_xasz(clcf, xr, shft=0.0, norm=True):
    # Data import
    E, Gz = np.loadtxt(clcf, unpack=True, skiprows=5, usecols=(0,2))

    # Spectral weight
    Gz = -(1/np.pi)*Gz

    # G(w)
    G = np.interp(xr, E - E[argrelmax(Gz)[0][0]] + shft, Gz)

    if norm:
        G = normalize(G)

    return G


def calc_xas_Ld(clc, xr, clcb=0, shft=0.0, shftb=0.0, norm=True):
    nc = 2
    if clcb:
        nc = 1

    G_R_a, G_A_a, _ = spectralweight(clc, xr, shft=shft, nc=nc)
    XAS_Ld = G_A_a.filter(regex='L', axis=1).sum(axis=1)
    XAS_d1 = G_A_a.filter(regex='1', axis=1).sum(axis=1)

    if clcb:
        G_R_b, G_A_b, _ = spectralweight(clcb, xr, shft=shftb, nc=nc)
        XAS_Ld += G_A_b.filter(regex='L', axis=1).sum(axis=1)
        XAS_d2 = G_A_b.filter(regex='1', axis=1).sum(axis=1)
    else:
        XAS_d2 = G_A_a.filter(regex='2', axis=1).sum(axis=1)

    if norm:
        XAS_Ld /= np.amax(XAS_Ld)
        tmp = np.amax(XAS_d1 + XAS_d2)
        XAS_d1 /= tmp
        XAS_d2 /= tmp

    return XAS_Ld, XAS_d1, XAS_d2


def sum_norm(t, norm=1.0):
    out = np.zeros_like(t[0])
    cht = np.zeros((len(t[0]), len(t)))

    for i in np.arange(len(t)):
        cht[:, i] = t[i]
        out += cht[:, i]

    for i in np.arange(len(t)):
        cht[:, i] *= norm/np.amax(out)

    out *= norm/np.amax(out)

    return tuple(np.column_stack((cht, out)).T)

def calc_xps_corelevel(clcf, xr, SO, doSO=True, shft=0.0, norm=True):
    rE, rG = np.loadtxt(clcf, unpack=True, skiprows=5, usecols=(0,2))
    rG = -(1/np.pi)*rG
    rE = rE - rE[np.argmax(rG)] + shft

    G_I_32 = np.interp(xr, rE, rG)
    G_I_12 = np.interp(xr, rE+SO, rG/2)

    if doSO:
        if norm:
            G_I_32, G_I_12, _ = sum_norm((G_I_32, G_I_12))
        return G_I_32, G_I_12
    else:
        if norm:
            G_I_32 /= np.amax(G_I_32)
        return G_I_32, 0

def calc_xps_vb(clc, xr, csec, G_R_nb, nc=1, clcb=False, shft=0.0, shftb=0.0, norm=True, sym=True, shiftEF=True):
    G_R, G_A, _ = spectralweight(clc, xr, shft=shft, invert=False, nc=nc, sym=sym, shiftEF=shiftEF)
    G_R_d1 = csec[0]*G_R.filter(regex='1', axis=1).sum(axis=1)
    G_R_Ld = csec[1]*G_R.filter(regex='L', axis=1).sum(axis=1)

    if clcb:
        G_R_b, G_A_b, _ = spectralweight(clcb, xr, shft=shftb, invert=False, nc=nc, sym=sym, shiftEF=shiftEF)
        G_R_Ld += csec[1]*G_R_b.filter(regex='L', axis=1).sum(axis=1)
        G_R_d2 = csec[2]*G_R_b.filter(regex='1', axis=1).sum(axis=1)
    elif nc == 2:
        G_R_d2 = csec[2]*G_R.filter(regex='2', axis=1).sum(axis=1)
    else:
        G_R_d2 = 0

    G_R_nb = csec[1]*G_R_nb

    if norm:
        G_R_d1, G_R_d2, G_R_Ld, G_R_nb, _ = sum_norm((G_R_d1, G_R_d2, G_R_Ld, G_R_nb))

    return G_R_d1, G_R_d2, G_R_Ld, G_R_nb


def spectralweight(specfile, xr, shft=0.0, wL=0.0, nc=2, fermi=False, Eg=0.0, invert=True, sym=True, shiftEF=True):
    t  = np.loadtxt(specfile, skiprows=5)
    E  = t[:, 0]
    G  = -(1/np.pi)*t[:, 2::2]

    if sym:
        ncols = 4
    else:
        ncols = 2

    if nc == 2:
        names = np.array([
                 'E',
                 '1eu', '1tu', '1ed', '1td',
                 'Leu', 'Ltu', 'Led', 'Ltd',
                 '2eu', '2tu', '2ed', '2td'
                ])
    else:
        if sym:
            names = np.array([
                     'E',
                     '1eu', '1tu', '1ed', '1td',
                     'Leu', 'Ltu', 'Led', 'Ltd',
                    ])
        else:
            names = np.array([ 'E', '1u', '1d', 'Lu', 'Ld'])

    if invert:
        G_R = pd.DataFrame(np.flipud(np.column_stack((-E, G[:, :ncols*(1+nc)]))), columns=names)
    else:
        G_R = pd.DataFrame(np.column_stack((E, G[:, :ncols*(1+nc)])), columns=names)

    G_A = pd.DataFrame(np.column_stack((E, G[:, ncols*(1+nc):])), columns=names)

    if sym:
        G_R.loc[:, ['Leu', 'Ltu', 'Led', 'Ltd']] *= nc
        G_A.loc[:, ['Leu', 'Ltu', 'Led', 'Ltd']] *= nc
    else:
        G_R.loc[:, ['Lu', 'Ld']] *= nc
        G_A.loc[:, ['Lu', 'Ld']] *= nc

    if invert:
        EF  = G_R['E'][argrelmax(G_R.loc[:, names[1:]].sum(axis=1).values)[0][-1]]
    else:
        EF  = G_R['E'][argrelmax(G_R.loc[:, names[1:]].sum(axis=1).values)[0][0]]

    if not shiftEF:
        EF = 0

    fR = []
    fA = []
    for n in names:
        if n != 'E':
            if wL > 0:
                fR.append(gaussian(np.array(G_R['E'].values) - EF + shft, np.array(G_R[n].values), xr, wL))
                fA.append(gaussian(np.array(G_A['E'].values) - EF + shft, np.array(G_A[n].values), xr, wL))
            else:
                fR.append(np.interp(xr, G_R['E'].values - EF + shft, G_R[n].values))
                fA.append(np.interp(xr, G_A['E'].values - EF + shft, G_A[n].values))

    G_R = np.array(fR).T
    G_A = np.array(fA).T


    G_R = pd.DataFrame(G_R, columns=names[1:])
    G_A = pd.DataFrame(G_A, columns=names[1:])

    if fermi:
        for n in names:
            if n != 'E':
                G_R[n] = G_R[n]/(np.exp( (xr+Eg/2)/fermi) + 1)
                G_A[n] = G_A[n]/(np.exp(-(xr-Eg/2)/fermi) + 1)

    return G_R, G_A, EF


def fermienergy(specfile, nc=2):
    t  = np.loadtxt(specfile, skiprows=5)
    E  = t[:, 0]
    G  = -(1/np.pi)*t[:, 2::2]

    G_R = pd.DataFrame(np.column_stack((E, G[:, :4*(1+nc)])))
    EF  = G_R[0][argrelmax(G_R.loc[:, 1:].sum(axis=1).values)[0][0]]

    return EF
