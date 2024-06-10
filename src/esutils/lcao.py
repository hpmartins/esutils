import numpy as np

def nst(p, pos):
    n = sst(p, pos)
    return np.arange(n, n+(3 if pos[p,0] == 1 else 5))

def sst(p, pos):
    return (3*(pos[:p,0]==1).sum() + 5*(pos[:p,0]==2).sum())

def slater(r, t1, t2, a, b, pps, ppp, pds, pdp):
    lmn     = r/np.linalg.norm(r)
    l, m, n = lmn

    if (t1 == 2) & (t2 == 2):
        return 0

    if (t1 == 1) & (t2 == 1):
        if a == b:
            return (lmn[a]**2)*pps + (1 - lmn[a]**2)*ppp
        else:
            return lmn[a]*lmn[b]*(pps - ppp)

    if (t1 == 2) & (t2 == 1):
        tmp = int(b)
        b = int(a)
        a = int(tmp)

    a += 1
    b += 1

    if (a == 1):
        if (b == 1):
            return (l*(n**2-(1/2)*(l**2+m**2))*pds-np.sqrt(3)*l*(n**2)*pdp)
        elif (b == 2):
            return ((np.sqrt(3)/2)*l*(l**2-m**2)*pds+l*(1-l**2+m**2)*pdp)
        elif (b == 3):
            return (np.sqrt(3)*(l**2)*m*pds+m*(1-2*l**2)*pdp)
        elif (b == 4):
            return (np.sqrt(3)*l*m*n*pds-2*l*m*n*pdp)
        elif (b == 5):
            return (np.sqrt(3)*(l**2)*n*pds+n*(1-2*l**2)*pdp)
    elif (a == 2):
        if (b == 1):
            return (m*(n**2-(1/2)*(l**2+m**2))*pds-np.sqrt(3)*m*(n**2)*pdp)
        elif (b == 2):
            return ((np.sqrt(3)/2)*m*(l**2-m**2)*pds-m*(1+l**2-m**2)*pdp)
        elif (b == 3):
            return (np.sqrt(3)*(m**2)*l*pds+l*(1-2*m**2)*pdp)
        elif (b == 4):
            return (np.sqrt(3)*(m**2)*n*pds+n*(1-2*m**2)*pdp)
        elif (b == 5):
            return (np.sqrt(3)*l*m*n*pds-2*l*m*n*pdp)
    elif (a == 3):
        if (b == 1):
            return (n*(n**2-(1/2)*(l**2+m**2))*pds+np.sqrt(3)*n*(l**2+m**2)*pdp)
        elif (b == 2):
            return ((np.sqrt(3)/2)*n*(l**2-m**2)*pds-n*(l**2-m**2)*pdp)
        elif (b == 3):
            return (np.sqrt(3)*l*m*n*pds-2*l*m*n*pdp)
        elif (b == 4):
            return (np.sqrt(3)*(n**2)*m*pds+m*(1-2*n**2)*pdp)
        elif (b == 5):
            return (np.sqrt(3)*(n**2)*l*pds+l*(1-2*n**2)*pdp)

def lcao_oct(pps, pds, D, Ep, Dq, pos=0, typ=''):
    if typ == 'single':
        pos = np.array([
                [1,  0,  0,  0,  0],
                [1,  0, -2,  0,  0],
                [1,  0, -1,  1,  0],
                [1,  0, -1, -1,  0],
                [1,  0, -1,  0,  1],
                [1,  0, -1,  0, -1],
                [2,  0, -1,  0,  0]
        ])
        pds = np.array([pds])
        D   = np.array([D])
        Dq  = np.array([Dq])
    elif typ == 'double':
        pos = np.array([
                [1, 0,  0,  0,  0],
                [1, 0, -2,  0,  0],
                [1, 0, -1,  1,  0],
                [1, 0, -1, -1,  0],
                [1, 0, -1,  0,  1],
                [1, 0, -1,  0, -1],

                [1, 1,  0,  0,  0],
                [1, 1,  2,  0,  0],
                [1, 1,  1,  1,  0],
                [1, 1,  1, -1,  0],
                [1, 1,  1,  0,  1],
                [1, 1,  1,  0, -1],

                [2, 0, -1,  0,  0],
                [2, 1,  1,  0,  0]
        ])
        pds = np.array(pds)
        D   = np.array(D)
        Dq  = np.array(Dq)

    ppp = -0.30*pps
    pdp = -0.45*pds
    Ed  = D + Ep

    sz = 3*(pos[:,0]==1).sum() + 5*(pos[:,0]==2).sum()
    H = np.zeros((sz,sz))

    char_l = []
    char_c = []
    for p in np.arange(len(pos)):
        char_l = char_l + (3 if pos[p,0]==1 else 5)*([1] if pos[p,0]==1 else [2])
        char_c = char_c + (3 if pos[p,0]==1 else 5)*[int(pos[p,1])]
        for q in np.arange(len(pos)):
            t1 = int(pos[p,0])
            t2 = int(pos[q,0])

            r = pos[q,2:] - pos[p,2:]
            n = np.linalg.norm(r)

            i = int(pos[p,1])

            if n == 0:
                for m in nst(p, pos):
                    if t1 == 1:
                        H[m, m] = Ep

                    if t1 == 2:
                        H[m, m] = Ed[i]
                        mm = m - sst(p, pos)
                        if mm <= 1:
                            H[m, m] += 6*Dq[i]
                        else:
                            H[m, m] -= 4*Dq[i]
            elif n <= np.sqrt(2):
                for ii in (np.arange(3) if (t1 == 1) else np.arange(5)):
                    for jj in (np.arange(3) if (t2 == 1) else np.arange(5)):
                        T = slater(r, t1, t2, ii, jj, pps, ppp, pds[i], pdp[i])
                        H[ii+sst(p,pos),jj+sst(q,pos)] = T
                        H[jj+sst(q,pos),ii+sst(p,pos)] = T

    eva, evc = np.linalg.eig(H)

    char_l = np.array(char_l)
    char_c = np.array(char_c)
    char_a = np.zeros((sz, 4))
    for q in np.arange(sz):
        vec = np.real(evc[:, q])

        for cn in np.unique(char_c):
            char_a[q, 2*cn  ] = ((vec[(char_l == 1) & (char_c == cn)])**2).sum()
            char_a[q, 2*cn+1] = ((vec[(char_l == 2) & (char_c == cn)])**2).sum()

    out = np.real(np.column_stack((eva, char_a)))
    oc = np.real(np.round(out, 3))
    
    E, countE = np.unique(oc[:, 0], return_counts=True)
    out = np.zeros((np.size(E), np.size(oc, axis=1) - 1))
    
    for i in np.arange(np.size(E)):
       out[i, :] = oc[oc[:, 0] == E[i], 1:].sum(axis=0)

    # out = out[out[:,0].argsort()]
    # E = out[:, 0]
    # out = out[:, 1:]
    return E, out


def get_oct_p(pps, pds, D, Ep, Dq, typ):
    E, listOct = lcao_oct(pps, pds, D, Ep, Dq, typ=typ)
    tmp = listOct[:, [0,2]].sum(axis=1)/listOct.sum(axis=1)
    listOct = listOct[tmp > 0.7, :]
    E       = -E[tmp > 0.7]
    listOct = listOct[:, [0,2]]
    return E, listOct
