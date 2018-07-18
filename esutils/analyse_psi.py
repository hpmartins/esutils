import numpy as np

psif = ''
e1 = 'Fe'
e2 = 'Mn'
nt = 30
n1 = 5
n2 = 4
ns = int(nt/3)

def fmt(n_d1, n_Ld, n_d2, pct, w):
    n_d1 = int(n_d1)
    n_Ld = 10 - int(n_Ld)
    n_d2 = int(n_d2)
    pct  = 100*pct
    s_Ld = str(n_Ld)
    if n_Ld <= 1:
        s_Ld = ''
    if n_Ld < 1:
        n_Ld = ''
    else:
        if w:
            if n_Ld > 1:
                n_Ld = '\\uL^' + s_Ld + ''
            else:
                n_Ld = '\\uL'
        else:
            n_Ld = 'L' + s_Ld
    if w:
        return '$d^{} {:5} d^{}$ & {:4.1f}\% \\\\'.format(n_d1, n_Ld, n_d2, pct)
    else:
        return 'd{} {:2} d{} {:4.1f}%'.format(n_d1, n_Ld, n_d2, pct)

psi = np.loadtxt(psif, skiprows=8, usecols=(1,2), converters={2: lambda s: int(s, 16)})
psi = psi[psi[:, 0] != 0, :]
psi[:, 0] = psi[:, 0]**2
psil = np.column_stack((psi[:, 0], [list(map(int, bin(int(i))[2:(2+nt)].zfill(nt))) for i in psi[:, 1]]))
nlist = np.column_stack((psil[:, 0], psil[:, 1:(1+ns)].sum(axis=1), psil[:, (1+ns):(1+2*ns)].sum(axis=1), psil[:, (1+2*ns):].sum(axis=1)))

d1d2_ctb = {}
for pct,occd1,occLd,occd2 in nlist:
    if occd1 not in d1d2_ctb:
        d1d2_ctb[occd1] = {}

    if occLd not in d1d2_ctb[occd1]:
        d1d2_ctb[occd1][occLd] = {}

    if occd2 not in d1d2_ctb[occd1][occLd]:
        d1d2_ctb[occd1][occLd][occd2] = 0

for pct,occd1,occLd,occd2 in nlist:
    d1d2_ctb[occd1][occLd][occd2] += pct

occ = [-1, 0, 0, 0]
for occd1 in d1d2_ctb.keys():
    for occLd in d1d2_ctb[occd1].keys():
        for occd2 in d1d2_ctb[occd1][occLd].keys():
            occ = np.vstack((occ, [occd1, occLd, occd2, d1d2_ctb[occd1][occLd][occd2]]))

occ = occ[(occ[:, 3] > 0.01*np.amax(occ[:, 3])) | ((occ[:, 1] == 10.0) & (occ[:, 3] > 0.001*np.amax(occ[:, 3]))), :]
occ = occ[np.argsort(occ[:, 3]), :]
occ = occ[::-1]
occ[:, 3] /= np.sum(occ[:, 3])

print('')
print('<{}> = {:4.1f}'.format(e1, (occ[:, 0]*occ[:, 3]).sum()))
print('<{}> = {:4.1f}'.format(e2, (occ[:, 2]*occ[:, 3]).sum()))
print('<Ld> = {:4.1f}'.format(2*(occ[:, 1]*occ[:, 3]).sum()))

print('\n{}: '.format(e1), end='')
for n in np.unique(occ[:, 0]):
    if n <= n1:
        nL = ''
    else:
        nL = 'L' + str(int(n - n1))
    print('{:.1f}% d{:d}{}, '.format(100*occ[occ[:, 0] == n, 3].sum(axis=0), int(n), nL), end='')

print('\n{}: '.format(e2), end='')
for n in np.unique(occ[:, 2]):
    if n <= n2:
        nL = ''
    else:
        nL = 'L' + str(int(n - n2))
    print('{:.1f}% d{:d}{}, '.format(100*occ[occ[:, 2] == n, 3].sum(axis=0), int(n), nL), end='')

print('\n')
for n_d1, n_Ld, n_d2, pct in occ:
    print(fmt(n_d1, n_Ld, n_d2, pct, 1))
