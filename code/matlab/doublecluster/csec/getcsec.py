from numpy import array
from re import search, findall
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#from sys import argv

nshell = {1: 'K', 2: 'L', 3: 'M', 4: 'N', 5: 'O'}

#Z  = argv[0]
#n  = argv[1]
#l  = argv[2]
#hn = argv[3]

# n = 2 e l = 1 -> 2p
# n = 3 e l = 2 -> 3d, etc
Z = 22
n = 3
l = 2
hn = 1840

Z = str(Z).zfill(2)
b = nshell[n] + str(l+1)


with open(Z + '.txt', 'r') as f:
	data = f.read()

data = data.replace('*', '').strip().replace('  ', ' ').replace('\n ', '\n').replace(' \n', '\n')
blocks = data.split('END OF DATA')[3:-1]

for block in blocks:
	bname = search('BLOCK:(.+)', block).groups()[0]
	if b == bname:
		table = array(findall('([0-9\.E\-\+]+) ([0-9\.E\-\+]+)\n', block), dtype=float)
		table[:, 0] = 1000*table[:, 0]
		f = interp1d(table[:, 0], table[:, 1])
		print(f(hn))
		break

print(table)
g = plt.figure()
asdf = g.add_subplot(1, 1, 1)
asdf.semilogy(table[:, 0], table[:, 1])
asdf.scatter(hn, f(hn))
asdf.set_xlim((1840, 2500))
asdf.set_ylim((3e2, 1e4))
g.show()