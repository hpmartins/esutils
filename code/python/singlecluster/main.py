import numpy as np
import scipy.linalg as sp
import yaml
import pandas as pd
import scipy.spatial as spatial
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix, lil_matrix
from scipy.special import wofz
from operator import itemgetter
from code import interact
from sys import exc_info
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=999999)
pd.set_option('display.max_rows', 1000)
pd.set_option('display.expand_frame_repr', True)
pd.set_option('display.max_colwidth', 300)
pd.set_option('display.width', 300)
pd.set_option('display.colheader_justify', 'left')
#np.set_printoptions(precision=4, formatter={'float': '{: 0.3f}'.format})

def fmt(x, y):
    return 'x: {x:0.2f}\ny: {y:0.2f}'.format(x=x, y=y)

class FollowDotCursor(object):
    """Display the x,y location of the nearest data point."""
    def __init__(self, ax, x, y, tolerance=5, formatter=fmt, offsets=(-20, 20)):
        try:
            x = np.asarray(x, dtype='float')
        except (TypeError, ValueError):
            x = np.asarray(mdates.date2num(x), dtype='float')
        y = np.asarray(y, dtype='float')
        self._points = np.column_stack((x, y))
        self.offsets = offsets
        self.scale = x.ptp()
        self.scale = y.ptp() / self.scale if self.scale else 1
        self.tree = spatial.cKDTree(self.scaled(self._points))
        self.formatter = formatter
        self.tolerance = tolerance
        self.ax = ax
        self.fig = ax.figure
        self.ax.xaxis.set_label_position('top')
        self.dot = ax.scatter(
            [x.min()], [y.min()], s=130, color='green', alpha=0.7)
        self.annotation = self.setup_annotation()
        plt.connect('motion_notify_event', self)

    def scaled(self, points):
        points = np.asarray(points)
        return points * (self.scale, 1)

    def __call__(self, event):
        ax = self.ax
        # event.inaxes is always the current axis. If you use twinx, ax could be
        # a different axis.
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
        elif event.inaxes is None:
            return
        else:
            inv = ax.transData.inverted()
            x, y = inv.transform([(event.x, event.y)]).ravel()
        annotation = self.annotation
        x, y = self.snap(x, y)
        annotation.xy = x, y
        annotation.set_text(self.formatter(x, y))
        self.dot.set_offsets((x, y))
        bbox = ax.viewLim
        event.canvas.draw()

    def setup_annotation(self):
        """Draw and hide the annotation box."""
        annotation = self.ax.annotate(
            '', xy=(0, 0), ha = 'right',
            xytext = self.offsets, textcoords = 'offset points', va = 'bottom',
            bbox = dict(
                boxstyle='round,pad=0.5', fc='yellow', alpha=0.75),
            arrowprops = dict(
                arrowstyle='->', connectionstyle='arc3,rad=0'))
        return annotation

    def snap(self, x, y):
        """Return the value in self.tree closest to x, y."""
        dist, idx = self.tree.query(self.scaled((x, y)), k=1, p=1)
        try:
            return self._points[idx]
        except IndexError:
            # IndexError: index out of bounds
            return self._points[0]

def keyboard(banner=None):
    try:
        raise None
    except:
        frame = exc_info()[2].tb_frame.f_back
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        interact(banner=banner, local=namespace)
    except SystemExit:
        return 

def obj_dic(d):
	top = type('new', (object,), d)
	seqs = tuple, set, frozenset
	for i, j in d.items():
		if isinstance(j, dict):
			setattr(top, i, obj_dic(j))
		elif isinstance(j, seqs):
			setattr(top, i, 
				type(j)(obj_dic(sj) if isinstance(sj, dict) else sj for sj in j))
		elif isinstance(j, list):
			setattr(top, i, np.array(j))
		else:
			setattr(top, i, j)
	return top

def wspin(i):
	return 1 if i<5 else -1

def normalize(A):
	return np.array([(x-min(A))/(max(A)-min(A)) for x in A])

def normalize2(A):
	return np.array([x/max(A) for x in A])

def export_pd(name, pd):
	pd.to_csv(name, index=False, float_format='%.6f')

def export_np(name, ar, h):
	np.savetxt(name, ar, fmt='%.6f', header=h, delimiter=',', comments='')

def import_exp(fl):
	try:
		x, y = np.loadtxt(fl, unpack=True, usecols=((0,1)))
	except FileNotFoundError:
		return False, False
		
	y = y[np.argsort(x)]
	x = x[np.argsort(x)]
	return x, y

def import_exp2(fl):
	try:
		x, y = np.loadtxt(fl, unpack=True, usecols=((0,1)))
		y = y[np.argsort(x)]
		x = x[np.argsort(x)]
	except FileNotFoundError:
		return []
		
	return np.vstack((x, y)).T
	
def calc_exp_bg(x, y, new_x):
	return np.interp(new_x, x[::-1], [y[0]*a for a in normalize(np.cumsum(y[::-1]))])

def voigt(lst, new_x, HWHM_G, HWHM_L, lor_a):
	if len(lst) <= 0:
		return np.zeros(new_x.shape)

	x = lst[:, 0]
	y = lst[:, 1]
		
	wL = [HWHM_L*(1+lor_a*(a-x[0])) for a in x]
	SIGMA = HWHM_G/np.sqrt(2*np.log(2));
	out = np.zeros(len(new_x))
	for p in range(len(new_x)):
		tmp = 0
		for q in range(len(y)):
			z = ((new_x[p]-x[q]) + 1j*wL[q])/SIGMA/np.sqrt(2)
			w = wofz(z)
			V = y[q]*(w.real)/SIGMA/np.sqrt(2*np.pi);
			tmp = tmp + V
		out[p] = tmp
	out = [x-min(out) for x in out]
	return np.array(out)

def cfgname(d, L, C, c):
	p1='d'+str(d)
	p2=('c' if c>0 else '')
	p3=('L'+(str(L) if L > 1 else '') if L>0 else '')
	p4=('C'+(str(C) if C > 1 else '') if C>0 else '')
	return p1+p2+p3+p4

def occupation(s, n, v):
	table = {}
	base = s.b if (n < 0) else s.b[n]
	vect = s.evc if (n < 0) else s.evc[n][:,v]
	for i in base[:,-1]:
		nm = cfgname(base[i,0], cnt_orb(base[i,2]), cnt_orb(base[i,3]), cnt_orb(base[i,4]))
		table[nm] = (table[nm]+vect[i]*vect[i]) if (nm in table) else vect[i]*vect[i]

	s_table = sorted(table.items(), key=itemgetter(1), reverse=True)
	s = ""
	for cfg, occ in s_table:
		s = s + '{:>7}: {:>5.2f}\n'.format(cfg, 100*occ)
	s = s + '{:>6} = {:>5.2f}\n'.format('<n>', np.dot(vect*base[:, 0], vect))
	return s

def lcao_oct(pm, pos):
	sz = 3*(pos[:,0]==1).sum() + 5*(pos[:,0]==2).sum()
	H = np.zeros((sz,sz))

	char_l = []
	for p in np.arange(len(pos)):
		char_l = char_l + (3 if pos[p,0]==1 else 5)*([1] if pos[p,0]==1 else [2])
		for q in np.arange(len(pos)):
			t1 = pos[p,0]
			t2 = pos[q,0]
			
			i = pos[p,1]
			j = pos[q,1]
			
			if i == j:
				for m in nst(p, pos):
					H[m, m] = pm.Ep if (t1 == 1) else (pm.Ed + (6*pm.Dq if m<=1 else -4*pm.Dq))
			else:
				r = pos[q,2:] - pos[p,2:]
				n = np.linalg.norm(r)
				if n < 1.5:
					for ii in (np.arange(3) if (t1 == 1) else np.arange(5)):
						for jj in (np.arange(3) if (t2 == 1) else np.arange(5)):
							T = slater(pm, r, t1, t2, ii, jj)
							H[ii+sst(p,pos),jj+sst(q,pos)] = T
							H[jj+sst(q,pos),ii+sst(p,pos)] = T

	eva, evc = sp.eigh(H)
	
	char_d = np.zeros(sz)
	char_p = np.zeros(sz)
	char_l = np.array(char_l)
	
	for q in np.arange(sz):
		vec = evc[:,q]
		char_p[q] = ((vec[char_l == 1])**2).sum()
		char_d[q] = ((vec[char_l == 2])**2).sum()
	
	
	if pm.rs_no_elec:
		is_pure_p = (np.around(char_p, 5) == 1.0)
		eva    = eva[is_pure_p]
		char_p = char_p[is_pure_p]
		char_d = char_d[is_pure_p]
		out = pd.DataFrame(np.vstack((eva, 2*char_p, 2*char_d, 2*(char_d+char_p))).T, columns=['E', 'I', 'd', 't'])
	else:
		tmp = np.vstack([eva, 2*char_p, 2*char_d, 2*(char_d+char_p), char_d/char_p]).T
		tmp = tmp[tmp[:,4].argsort(),]
		tmp = tmp[:18,:4]
		out = pd.DataFrame(tmp, columns=['E', 'I', 'd', 't'])

	out.loc[:, 'E'] = np.around(out.E.values, decimals=3)
	out.loc[:, 'I'] = np.around(out.I.values, decimals=3)
	out.loc[:, 'd'] = np.around(out.d.values, decimals=3)
	out.loc[:, 't'] = np.around(out.t.values, decimals=3)
	out = out.groupby('E', as_index=False).sum().sort_values('I').loc[:, ['E', 'I', 'd', 't']]
	out = out.loc[:, ['E', 'I']]
	return out

def nst(p, pos):
	n = sst(p, pos)
	return np.arange(n, n+(3 if pos[p,0] == 1 else 5))

def sst(p, pos):
	return (3*(pos[:p,0]==1).sum() + 5*(pos[:p,0]==2).sum())

def slater(pm, r, t1, t2, a, b):
	lmn     = r/np.linalg.norm(r)
	l, m, n = lmn
	
	if (t1 == 2) & (t2 == 2):
		return 0
	
	if (t1 == 1) & (t2 == 1):
		if a == b:
			return (lmn[a]**2)*pm.pps + (1 - lmn[a]**2)*pm.ppp
		else:
			return lmn[a]*lmn[b]*(pm.pps - pm.ppp)
	
	if (t1 == 2) & (t2 == 1):
		tmp = int(b)
		b = int(a)
		a = int(tmp)
	
	a += 1
	b += 1
	
	if (a == 1):
		if (b == 1):
			return (l*(n**2-(1/2)*(l**2+m**2))*pm.pds-np.sqrt(3)*l*(n**2)*pm.pdp)
		elif (b == 2):
			return ((np.sqrt(3)/2)*l*(l**2-m**2)*pm.pds+l*(1-l**2+m**2)*pm.pdp)
		elif (b == 3):
			return (np.sqrt(3)*(l**2)*m*pm.pds+m*(1-2*l**2)*pm.pdp)
		elif (b == 4):
			return (np.sqrt(3)*l*m*n*pm.pds-2*l*m*n*pm.pdp)
		elif (b == 5):
			return (np.sqrt(3)*(l**2)*n*pm.pds+n*(1-2*l**2)*pm.pdp)
	elif (a == 2):
		if (b == 1):
			return (m*(n**2-(1/2)*(l**2+m**2))*pm.pds-np.sqrt(3)*m*(n**2)*pm.pdp)
		elif (b == 2):
			return ((np.sqrt(3)/2)*m*(l**2-m**2)*pm.pds-m*(1+l**2-m**2)*pm.pdp)
		elif (b == 3):
			return (np.sqrt(3)*(m**2)*l*pm.pds+l*(1-2*m**2)*pm.pdp)
		elif (b == 4):
			return (np.sqrt(3)*(m**2)*n*pm.pds+n*(1-2*m**2)*pm.pdp)
		elif (b == 5):
			return (np.sqrt(3)*l*m*n*pm.pds-2*l*m*n*pm.pdp)
	elif (a == 3):
		if (b == 1):
			return (n*(n**2-(1/2)*(l**2+m**2))*pm.pds+np.sqrt(3)*n*(l**2+m**2)*pm.pdp)
		elif (b == 2):
			return ((np.sqrt(3)/2)*n*(l**2-m**2)*pm.pds-n*(l**2-m**2)*pm.pdp)
		elif (b == 3):
			return (np.sqrt(3)*l*m*n*pm.pds-2*l*m*n*pm.pdp)
		elif (b == 4):
			return (np.sqrt(3)*(n**2)*m*pm.pds+m*(1-2*n**2)*pm.pdp)
		elif (b == 5):
			return (np.sqrt(3)*(n**2)*l*pm.pds+l*(1-2*n**2)*pm.pdp)

def set_orb(v, idx, x):
	mask = 1 << (9 - idx)
	v &= ~mask
	if x:
		v |= mask
	return v

def tog_orb(v, idx):
	return v^(1<<(9-idx))

def has_orb(v, idx):
	return ((v & (1 << (9 - idx))) != 0)
	
def cnt_orb(v):
	return bin(v).count('1')

class Cluster:	
	def obj_dic(self, d):
		top = type('clustercfg', (object,), d)
		seqs = tuple, set, frozenset
		for i, j in d.items():
			if isinstance(j, dict):
				setattr(top, i, self.obj_dic(j))
			elif isinstance(j, seqs):
				setattr(top, i, 
					type(j)(self.obj_dic(sj) if isinstance(sj, dict) else sj for sj in j))
			elif isinstance(j, list):
				setattr(top, i, np.array(j))
			else:
				setattr(top, i, j)
		return top
	
	def kanamori(self):
		F2 = self.pm.F_red*self.pm.F2/49
		F4 = self.pm.F_red*self.pm.F4/441
		
		C = 35*F4
		B = F2 - 5*F4
		A = self.pm.U + (14/9)*B - (7/9)*C
		
		self.pm.u = A + 4*B + 3*C
		self.pm.ul = A - B + C
		self.pm.j = (5/2)*B + C
	
	def __init__(self, f):
		pm = self.obj_dic(yaml.safe_load(open(f)))
		pm.ppp = -0.3*pm.pps
		pm.pds = pm.Ts/np.sqrt(3)
		pm.pdp = -0.45*pm.pds
		pm.Tp = -0.5248*pm.Ts
		pm.Ed = pm.D + pm.Ep - np.sum(pm.start)*pm.U
		pm.Ec = pm.Ed - pm.Dc + np.sum(pm.start)*pm.U
		tmp = np.array([pm.Tp, pm.Tp, pm.Tp, pm.Ts, pm.Ts, pm.Tp, pm.Tp, pm.Tp, pm.Ts, pm.Ts])
		pm.T_list = np.multiply(pm.Tm, tmp)
		self.pm = pm
		self.kanamori()
		if (self.pm.J < 0):
			self.pm.J = self.pm.j


class ClusterState:
	def __init__(self, pm, wc):
		self.pm = pm
		self.wc = wc

	def aux_gen_basis(self, start):
		# Base normal
		nst = np.sum(start)
		bs = int(''.join(str(x) for x in start), 2)
		allcfgs = np.arange(1024)
		states = np.array([a for a in allcfgs[(allcfgs & bs) == bs] if cnt_orb(a) <= nst+self.pm.maxholes])
		base_def = np.zeros((states.shape[0], 5), dtype=int)
		base_def[:, 1] = states
		base_def[:, 2] = states^bs
		
		# Base coerente
		if self.pm.coherent.nonzero():
			base_coherent = np.zeros([2*base_def.shape[0], 5], dtype=int)
			nn = 0
			for i in self.pm.coherent.nonzero()[0]:
				if has_orb(bs, i):
					base_coherent[nn, :] = [nst-1, tog_orb(int(bs), i), 0, 0, set_orb(0, i, 1)]
					nn = nn + 1
				
				if ~has_orb(bs, i):
					cl1 = np.array([a for a in base_def[:, 1] if has_orb(a, i)])
					cl2 = np.array([set_orb(a, i, 0) for a in base_def[:, 2] if has_orb(a, i)])
					cl3 = np.array([set_orb(0, i, 1) for a in base_def[:, 2] if has_orb(a, i)])
					try:
						base_coherent[nn:(nn+cl1.shape[0]), :] = np.array(np.vstack((cl1, cl1, cl2, cl3, np.zeros((cl1.shape[0])))).T, dtype=int)
					except:
						keyboard()
					nn = nn + cl1.shape[0]

			base_coherent = base_coherent[0:nn, :]
		
		out = np.vstack((base_def, base_coherent))
		out[:, 0] = np.array([cnt_orb(a) for a in out[:, 1]])
		return np.column_stack((out[out[:, 0].argsort(), :], np.arange(out.shape[0])))


	def gen_basis(self):
		if (self.wc == 0 or self.wc == 2):
			self.b = self.aux_gen_basis(self.pm.start)
		elif self.wc == 1:
			self.b = {}
			for r in np.where(self.pm.start == 0)[0]:
				new_start = list(self.pm.start)
				new_start[r] = 1
				self.b[r] = self.aux_gen_basis(new_start)
		elif self.wc == -1:
			self.b = {}
			for r in np.arange(10):
				if self.pm.start[r]:
					new_start = list(self.pm.start)
					new_start[r] = 0
					self.b[r] = self.aux_gen_basis(new_start)
				else:
					if (self.pm.rs_no_elec or not(self.pm.start.any())):
						tmp = self.aux_gen_basis(self.pm.start)
						self.b[r] = tmp[~has_orb(tmp[:, 1], r), :]	
						self.b[r][:, 2] = set_orb(self.b[r][:, 2], r, 1)
						self.b[r][:, 5] = np.arange(self.b[r].shape[0])


	def fill_diag(self, i):
		nEd = i[0]
		nU  = nEd*(nEd-1)/2
		nDq = 6*cnt_orb(i[1] & 0b0001100011) - 4*cnt_orb(i[1] & 0b1110011100)
		nup = cnt_orb(i[1] & 0b1111100000)
		ndn = cnt_orb(i[1] & 0b0000011111)
		nJ  = (nup*(nup-1))/2 + (ndn*(ndn-1))/2
		nEp = cnt_orb(i[2])
		nEc = cnt_orb(i[3]) - cnt_orb(i[4])
		
		out = nEd*self.pm.Ed - nEp*self.pm.Ep - nEc*self.pm.Ec + nDq*self.pm.Dq - nJ*self.pm.J
		
		if self.wc == 2:
			out -= nEd*self.pm.U/0.83
		
		if self.pm.pp:
			npp = cnt_orb(i[2] & 0b1110011100) - cnt_orb(i[2] & 0b0001100011)
			out += npp*(self.pm.pps - self.pm.ppp)
		
		if self.pm.kanamori:
			nu = cnt_orb((i[1] & 0b1111100000) & ((i[1] & 0b0000011111) << 5))
			nul = nU - nu
			out += nu*self.pm.u + nul*self.pm.ul
		else:
			out += nU*self.pm.U
		
		return out
		
		
	def aux_gen_hamiltonian(self, b):
		ijH = np.zeros((120*b.shape[0], 3))
		nn = 0
		nl = np.unique(b[:, 0])
		
		comb = {}
		for p in nl:
			comb[p] = b[b[:, 0] == p, :]

		for p in nl:
			for i in comb[p]:
				ijH[nn, :] = [i[-1], i[-1], self.fill_diag(i)]
				nn = nn + 1
				if p < nl[-1]:
					for j in comb[p+1]:
						d = int(i[1]^j[1])
						if (bin(d).count('1') == 1):
							k = j - i
							if (k[2] == d) and (k[3] == 0) and (k[4] == 0):
								T = self.pm.T_list[9 - (d.bit_length() - 1)]
								ijH[nn, :] = [i[-1], j[-1], T]
								ijH[nn+1, :] = [j[-1], i[-1], T]
								nn += 2
							if (k[2] == 0) and (k[3] == d) and (k[4] == 0):
								T = self.pm.Tc
								ijH[nn, :] = [i[-1], j[-1], T]
								ijH[nn+1, :] = [j[-1], i[-1], T]
								nn += 2
							if (k[2] == 0) and (k[3] == 0) and (k[4] == -d):
								T = self.pm.Tc
								ijH[nn, :] = [i[-1], j[-1], T]
								ijH[nn+1, :] = [j[-1], i[-1], T]
								nn += 2

		ijH = ijH[0:nn, :]
		return coo_matrix((ijH[:, 2], (ijH[:, 0], ijH[:, 1])), shape=(b.shape[0], b.shape[0])).tocsr()
		
	def gen_hamiltonian(self):
		if ((self.wc == 0) or (self.wc == 2)):
			self.h = self.aux_gen_hamiltonian(self.b)
		elif self.wc == 1:
			self.h = {}
			for r in np.where(self.pm.start == 0)[0]:
				self.h[r] = self.aux_gen_hamiltonian(self.b[r])
		elif self.wc == -1:
			self.h = {}
			for r in np.arange(10):
				if self.pm.start[r]:
					self.h[r] = self.aux_gen_hamiltonian(self.b[r])
				else:
					if self.pm.rs_no_elec:
						self.h[r] = self.aux_gen_hamiltonian(self.b[r])
						for i in np.arange(self.h[r].shape[0]):
							self.h[r][i,i] = self.h[r][i,i] - (self.pm.D - self.pm.U) ### Pq???????
							# Pq fizemos isso? Supomos que o cluster proximo ao simples
							# seria um simples, entao seria o mesmo shift?
						

	def diagonalize(self):
		if ((self.wc == 0) or (self.wc == 2)):
			eva, evc = sp.eigh(self.h.toarray())
			if self.wc == 0:
				self.eva = eva[0]
				self.evc = evc[:, 0]
			else:
				self.eva = eva
				self.evc = evc
		else:
			self.eva = {}
			self.evc = {}
			for r, m in self.h.items():
				self.eva[r], self.evc[r] = sp.eigh(m.toarray())


class ClusterOperator:
	def __init__(self, si, sf):
		bi = si.b
		bf = sf.b
		self.d = {}
		self.p = {}
		self.c = {}
		
		for n in bf.keys():
			self.d[n] = lil_matrix((bi.shape[0], bf[n].shape[0]))
			self.p[n] = lil_matrix((bi.shape[0], bf[n].shape[0]))
			self.c[n] = lil_matrix((bi.shape[0], bf[n].shape[0]))
			for j in bf[n]:
				for i in bi:
					k = j - i
					#k[0]: numero de eletrons
					#k[1]: configuracao
					#k[2]: posicoes dos L
					#k[3]: posicoes dos C
					#k[4]: posicoes dos c
					if (k[1] == 0) and (bin(k[2]).count('1') == 1) and (k[3] == 0) and (k[4] == 0):
						self.p[n][i[-1],j[-1]] = 1
					if (bin(k[1]).count('1') == 1) and (k[0]/sf.wc == 1) and (k[2] == 0):
						if (k[3] == 0) and (k[4] == 0):
							self.d[n][i[-1],j[-1]] = 1
						if (bin(k[3]).count('1') == 1) and (k[4] == 0):
							self.c[n][i[-1],j[-1]] = 1
						if (k[3] == 0) and (bin(k[4]).count('1') == 1):
							self.c[n][i[-1],j[-1]] = 1
			self.d[n].tocsr()
			self.p[n].tocsr()
			self.c[n].tocsr()
	
			
class ClusterTransition:
	def __init__(self, i, o, f):
		if f.wc == 2:
			self.E = f.eva - min(f.eva)
			self.I = np.square(np.dot(f.evc.T, i.evc))
		else:
			self.p1 = 1e4
			self.d = type('', (), {})()
			self.p = type('', (), {})()
			self.c = type('', (), {})()
			self.d.E = {}
			self.p.E = {}
			self.c.E = {}
			self.d.I = {}
			self.p.I = {}
			self.c.I = {}
			for n in f.eva.keys():
				self.d.E[n] = f.wc*(f.eva[n] - i.eva)
				self.p.E[n] = f.wc*(f.eva[n] - i.eva)
				self.c.E[n] = f.wc*(f.eva[n] - i.eva)
				self.p1 = min(self.p1, np.min(f.eva[n] - i.eva))
				self.d.I[n] = np.asarray(np.square(np.dot(np.dot(i.evc.T, o.d[n].todense()), f.evc[n]))).squeeze()
				self.p.I[n] = np.asarray(np.square(np.dot(np.dot(i.evc.T, o.p[n].todense()), f.evc[n]))).squeeze()
				self.c.I[n] = np.asarray(np.square(np.dot(np.dot(i.evc.T, o.c[n].todense()), f.evc[n]))).squeeze()
				
				self.d.E[n] = self.d.E[n][self.d.I[n]/self.d.I[n].max() >= 1e-3]
				self.p.E[n] = self.p.E[n][self.p.I[n]/self.p.I[n].max() >= 1e-3]
				self.c.E[n] = self.c.E[n][self.c.I[n]/self.c.I[n].max() >= 1e-3]
				self.d.I[n] = self.d.I[n][self.d.I[n]/self.d.I[n].max() >= 1e-3]
				self.p.I[n] = self.p.I[n][self.p.I[n]/self.p.I[n].max() >= 1e-3]
				self.c.I[n] = self.c.I[n][self.c.I[n]/self.c.I[n].max() >= 1e-3]

def filter(E, I):
	tmp = np.where(I/np.max(I) >= 1e-3)
	return E[tmp], I[tmp], tmp[0]

def cones(a):
	return np.ones(a.shape[0], dtype=int)

def occupation_str(s, n, v):
	table = {}
	
	if ((n < 0) and (v < 0)):
		base = s.b
		vect = s.evc
	
	if ((n < 0) and (v >= 0)):
		base = s.b
		vect = s.evc[:, v]
	
	if ((n >= 0) and (v >= 0)):
		base = s.b[n]
		vect = s.evc[n][:, v]
		
	for i in base[:,-1]:
		nm = cfgname(base[i,0], cnt_orb(base[i,2]), cnt_orb(base[i,3]), cnt_orb(base[i,4]))
		table[nm] = (table[nm]+vect[i]*vect[i]) if (nm in table) else vect[i]*vect[i]

	s_table = sorted(table.items(), key=itemgetter(1), reverse=True)
	out = ''
	for cfg, occ in s_table:
		if (100*occ > 1e-2):
			if (len(out) > 0):
				out = out + '; '
			out = out + '{}: {:.2f}'.format(cfg, 100*occ)
	#return '{:.2f} // {}'.format(np.dot(vect*base[:, 0], vect), out)
	return out
	
def transitions(pm, gnd, cor, o_emi, rem, o_abs, add):
	# Core
	E, I, idx = filter(cor.eva - min(cor.eva) + pm.coreopt.p, np.square(np.dot(cor.evc.T, gnd.evc)))
	E = np.around(E, decimals=3)
	cor_comp = np.array([occupation_str(cor, -1, v) for v in idx])
	A_core = pd.DataFrame({'spec': 'C',
	                       'E': np.concatenate([E, E + pm.coreopt.so]),
	                       'I': np.concatenate([I, I/2]),
					       'orb': np.array(np.concatenate([3*cones(E), cones(E)]), dtype=int),
						   'char': np.nan,
						   'comp': np.concatenate([cor_comp, cor_comp])})
		
	# Removal
	E = []
	I = []
	o = []
	c = []
	cmp = []
	for n in rem.eva.keys():
		#d
		tE, tI, tidx = filter(rem.wc*(rem.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_emi.d[n].todense()), rem.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['d']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(rem, n, v) for v in tidx])])
		
		#p
		tE, tI, tidx = filter(rem.wc*(rem.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_emi.p[n].todense()), rem.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['p']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(rem, n, v) for v in tidx])])
		
		#c
		tE, tI, tidx = filter(rem.wc*(rem.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_emi.c[n].todense()), rem.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['c']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(rem, n, v) for v in tidx])])

	A_rem = pd.DataFrame({'spec': 'R', 'E': E, 'I': I, 'orb': np.array(o, dtype=int), 'char': c, 'comp': cmp})

	# addition
	E = []
	I = []
	o = []
	c = []
	cmp = []
	for n in add.eva.keys():
		#d
		tE, tI, tidx = filter(add.wc*(add.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_abs.d[n].todense()), add.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['d']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(add, n, v) for v in tidx])])
		
		#p
		tE, tI, tidx = filter(add.wc*(add.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_abs.p[n].todense()), add.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['p']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(add, n, v) for v in tidx])])
		
		#c
		tE, tI, tidx = filter(add.wc*(add.eva[n] - gnd.eva), np.asarray(np.square(np.dot(np.dot(gnd.evc.T, o_abs.c[n].todense()), add.evc[n]))).squeeze())
		E = np.concatenate([E, tE])
		I = np.concatenate([I, tI])
		o = np.concatenate([o, (n+1)*cones(tE)])
		c = c + ['c']*tE.shape[0]
		cmp = np.concatenate([cmp, np.array([occupation_str(add, n, v) for v in tidx])])
	
	A_add = pd.DataFrame({'spec': 'A', 'E': E, 'I': I, 'orb': np.array(o, dtype=int), 'char': c, 'comp': cmp})
	
	shift = max(A_rem.E.values) + pm.vbopt.p
	A_rem.loc[:, 'E'] = A_rem.E.values - shift
	A_add.loc[:, 'E'] = A_add.E.values - shift
	return pd.concat([A_core, A_rem, A_add])

def write_output(C):
	s = """SINGLECLUSTER CALCULATION: {}

# Ground state properties
{}
Gap: {}

# Core transitions
{}

# Removal transitions
{}

# Addition transitions
{}

### Operators ###
{}
"""
	COR = C.A[C.A.spec == 'C'].groupby(['E', 'comp'], as_index=False).sum().loc[:, ['E', 'I', 'comp']]
	REM = C.A[C.A.spec == 'R'].groupby(['E', 'comp'], as_index=False).sum().loc[:, ['E', 'I', 'comp']]
	XAS = C.A[(C.A.spec == 'A') & (C.A.char == 'p')].groupby(['E', 'comp'], as_index=False).sum().loc[:, ['E', 'I', 'comp']]
	XAS.loc[:, 'E'] = XAS.E.values - XAS.E.min() + C.pm.xasopt.p
	gap = C.A[C.A.spec == 'A'].E.values.min() - C.A[C.A.spec == 'R'].E.values.max()
	
	s_op = ''
	for n in C.rem.eva.keys():
		s_op = s_op + "\n# ORBITAL {}\n".format(n)
		
		s_op = s_op + "# Removal d\n"
		idxlst = np.nonzero(C.emi.d[n])
		for i in np.arange(len(idxlst[0])):
			s_op = s_op + "{} -> {}\n".format(bin(C.gnd.b[idxlst[0][i], 1])[2:].zfill(10), bin(C.rem.b[n][idxlst[1][i], 1])[2:].zfill(10))
			
		s_op = s_op + "# Removal c\n"
		idxlst = np.nonzero(C.emi.c[n])
		for i in np.arange(len(idxlst[0])):
			s_op = s_op + "{} -> {}\n".format(bin(C.gnd.b[idxlst[0][i], 1])[2:].zfill(10), bin(C.rem.b[n][idxlst[1][i], 1])[2:].zfill(10))
			
		s_op = s_op + "# Removal p\n"
		idxlst = np.nonzero(C.emi.p[n])
		for i in np.arange(len(idxlst[0])):
			s_op = s_op + "{} -> {}\n".format(bin(C.gnd.b[idxlst[0][i], 1])[2:].zfill(10), bin(C.rem.b[n][idxlst[1][i], 1])[2:].zfill(10))

	with open(C.cname + '.out', 'w') as f:
		f.write(s.format(C.cname, occupation(C.gnd, -1, -1), gap, COR, REM, XAS, s_op))


def run_cluster(cname):
	# Parameters
	vnop = 300  # number of energy points - TODO: personalizavel no input pra cada caso
	swbw = 0.3  # bar width para sw, vb, cs, xas
	vbbw = 0.3
	csbw = 0.5
	xsbw = 0.3
	nrnd = 3    # numero de decimais ao arredondar a energia
	itol = 1e-5 # transicoes menores que itol*max sao cortadas
	
	# Cria cluster
	C = Cluster(cname + '.in')
	C.cname = cname
	if C.pm.plot:
		f = plt.figure(figsize=(10,8))

	# Cria estados
	C.gnd = ClusterState(C.pm,  0)
	C.cor = ClusterState(C.pm,  2)
	C.rem = ClusterState(C.pm, -1)
	C.add = ClusterState(C.pm,  1)

	# Gera bases
	C.gnd.gen_basis()
	C.cor.gen_basis()
	C.rem.gen_basis()
	C.add.gen_basis()

	# Gera hamiltonianos
	C.gnd.gen_hamiltonian()
	C.cor.gen_hamiltonian()
	C.rem.gen_hamiltonian()
	C.add.gen_hamiltonian()

	# Diagonaliza
	C.gnd.diagonalize()
	print(occupation(C.gnd, -1, -1))
	# occupation(S, N, V)
	# S: C.gnd, C.cor. C.rem, C.add
	# N: qual orbital, de 0 a 9 (pro caso de removal e addition)
	# V: qual vetor, o numero dele, estilo matlab mesmo
	# occupation(C.rem, 0, 1) -> ocupacao do primeiro autovetor de remocao do primeiro t2g down.
	# Fica legal se ligar a flag dev, dai o terminal fica aberto no fim estilo keyboard e eh
	# so usar estilo matlab, occupation(...)
	C.cor.diagonalize()
	C.rem.diagonalize()
	C.add.diagonalize()

	# Operadores
	C.abs = ClusterOperator(C.gnd, C.add)
	C.emi = ClusterOperator(C.gnd, C.rem)

	# Transicoes
	C.tabs = ClusterTransition(C.gnd, C.abs, C.add)
	C.temi = ClusterTransition(C.gnd, C.emi, C.rem)
	C.A = transitions(C.pm, C.gnd, C.cor, C.emi, C.rem, C.abs, C.add)
	C.A.loc[:, 'E'] = np.around(C.A.E.values, decimals=nrnd)
	
	# Octaedro
	pos = [[1, 1,  1,  0,  0],
		   [1, 2,  0,  1,  0],
		   [1, 3, -1,  0,  0],
		   [1, 4,  0, -1,  0],
		   [1, 5,  0,  0,  1],
		   [1, 6,  0,  0, -1],
		   [2, 7,  0,  0,  0]]
	OCT = lcao_oct(C.pm, np.array(pos))

	###### VALENCE BAND ######
	# Experimental data
	VB_EXP = import_exp2(C.pm.vbopt.exp)
	
	# Init
	VB = C.A[(C.A.spec == 'R') & (C.A.char != 'c')]
	VB_D = VB.loc[VB.char == 'd', ['E', 'I']].groupby('E', as_index=False).sum()
	VB_P = VB.loc[VB.char == 'p', ['E', 'I']].groupby('E', as_index=False).sum()
	VB_OCT = OCT.copy()
	
	# Energy signal
	VB_D.loc[:, 'E'] = -VB_D.E.values
	VB_P.loc[:, 'E'] = -VB_P.E.values
	VB_OCT.loc[:, 'E'] = -VB_OCT.E.values
	
	# Cross-sections
	csecs = normalize2([C.pm.vbopt.csecd, C.pm.vbopt.csecp])
	VB_D.loc[:, 'I'] = csecs[0]*VB_D.I.values
	VB_P.loc[:, 'I'] = csecs[1]*VB_P.I.values
	VB_OCT.loc[:, 'I'] = csecs[1]*C.pm.vbopt.octmult*VB_OCT.I.values

	# Energy range
	if len(VB_EXP) > 0:
		VB_E = np.linspace(min(VB_EXP[:,0]), max(VB_EXP[:,0]), vnop)
		if C.pm.vbopt.bg:
			VB_BG = VB_EXP[-1, 1]*normalize(np.interp(VB_E, VB_EXP[:,0], np.cumsum(VB_EXP[:,1])))
		else:
			VB_BG = np.zeros(VB_E.shape) 
	else:
		VB_E = np.linspace(-2, 15, vnop)
		VB_BG = np.zeros(VB_E.shape) # TODO: background do calculo
	
	# Voigts
	VOIGT_VB_D = voigt(VB_D.sort_values('E').values, VB_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	VOIGT_VB_P = voigt(VB_P.sort_values('E').values, VB_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	VOIGT_OCT  = voigt(VB_OCT.sort_values('E').values, VB_E, C.pm.vbopt.wG_oct, C.pm.vbopt.wL, C.pm.vbopt.aL)

	# Normalization
	max_voigts = max(VOIGT_VB_D + VOIGT_VB_P + VOIGT_OCT)
	VOIGT_VB_D = max(VOIGT_VB_D)*normalize(VOIGT_VB_D)/max_voigts
	VOIGT_VB_P = max(VOIGT_VB_P)*normalize(VOIGT_VB_P)/max_voigts
	VOIGT_OCT  = max(VOIGT_OCT)*normalize(VOIGT_OCT)/max_voigts
	VOIGT_T = normalize(VOIGT_VB_D + VOIGT_VB_P + VOIGT_OCT)
	max_bars = 2*max(VB_D.I.max(), VB_P.I.max(), VB_OCT.I.max())
	VB_D.loc[:, 'I'] = VB_D.I.values/max_bars
	VB_P.loc[:, 'I'] = VB_P.I.values/max_bars
	VB_OCT.loc[:, 'I'] = VB_OCT.I.values/max_bars
	vb_mult_bg = 1/max(VOIGT_T + VB_BG)
	VOIGT_T = normalize(VOIGT_T + VB_BG)
	VOIGT_VB_D = vb_mult_bg*VOIGT_VB_D
	VOIGT_VB_P = vb_mult_bg*VOIGT_VB_P
	VOIGT_OCT  = vb_mult_bg*VOIGT_OCT
	VB_BG = VOIGT_T[-1]*normalize(VB_BG)

	# Export
	if C.pm.export:
		export_pd('out/'+cname+'_VB_BARS.OUT', pd.concat([pd.DataFrame(VB_D.values, columns=['Ed', 'Id']), pd.DataFrame(VB_P.values, columns=['Ep', 'Ip']), pd.DataFrame(VB_OCT.values, columns=['Eoct', 'Ioct'])], axis=1))
		export_np('out/'+cname+'_VB_VOIGT.OUT', np.vstack((VB_E, VOIGT_VB_D, VOIGT_VB_P, VOIGT_OCT, VOIGT_T, VB_BG)).T, 'E,Id,Ip,Ioct,Itot,Ibg')
	
	# PLOT
	if C.pm.plot:
		# Figure initialization
		vb_plt = f.add_subplot(2, 2, 3)
		vb_plt.set_title('Valence Band')
		vb_plt.set_xlabel('Binding Energy (eV)')
		vb_plt.set_ylabel('Norm. Intensity')
		vb_plt.yaxis.set_major_locator(plt.NullLocator())
		vb_plt.set_xlim((min(VB_E), max(VB_E)))
		vb_plt.invert_xaxis()
		
		# Experimental data
		if len(VB_EXP) > 0:
			vb_plt.scatter(VB_EXP[:,0], normalize(VB_EXP[:,1]), color='black')
		
		# Transitions
		vb_plt.bar(VB_P.values[:,0], VB_P.values[:,1], width=vbbw, align='center', color='blue')
		vb_plt.bar(VB_D.values[:,0], VB_D.values[:,1], width=vbbw, align='center', color='red')
		vb_plt.bar(VB_OCT.values[:,0], VB_OCT.values[:,1], width=vbbw, align='center', color='cyan')
		
		# Voigts
		vb_plt.plot(VB_E, VOIGT_VB_D, aa=True, color='red', lw=1)
		vb_plt.plot(VB_E, VOIGT_VB_P, aa=True, color='blue', lw=1)
		vb_plt.plot(VB_E, VOIGT_OCT, aa=True, color='cyan', lw=1)
		vb_plt.plot(VB_E, VOIGT_T, aa=True, color='black', lw=2)
		vb_plt.plot(VB_E, VB_BG, aa=True, color='black', lw=0.5, linestyle='dashed')

	###### XAS ######
	# Experimental data
	XAS_EXP = import_exp2(C.pm.xasopt.exp)
	
	# Init
	XAS = C.A[(C.A.spec == 'A') & (C.A.char == 'p')].loc[:, ['E', 'I']].groupby('E', as_index=False).sum()
	XAS.loc[:, 'E'] = XAS.E.values - XAS.E.min() + C.pm.xasopt.p
	
	# Energy range
	if len(XAS_EXP) > 0:
		XAS_E = np.linspace(min(XAS_EXP[:,0]), max(XAS_EXP[:,0]), vnop)
	else:
		XAS_E = np.linspace(XAS.E.min() - 10, XAS.E.max() + 10, vnop)
	
	# Voigts
	VOIGT_XAS = voigt(XAS.values, XAS_E, C.pm.xasopt.wG, C.pm.xasopt.wL, C.pm.xasopt.aL)
	
	# Export
	if C.pm.export:
		export_pd('out/'+cname+'_XAS_BARS.OUT', XAS)
		export_np('out/'+cname+'_XAS_VOIGT.OUT', np.vstack((XAS_E, VOIGT_XAS)).T, 'E,I')
	
	# PLOT
	if C.pm.plot:
		# Figure initialization
		xas_plt = f.add_subplot(2, 2, 4)
		xas_plt.set_title('XAS')
		xas_plt.set_xlabel('Photon Energy (eV)')
		xas_plt.yaxis.set_major_locator(plt.NullLocator())
		xas_plt.set_xlim((min(XAS_E), max(XAS_E)))
		xas_plt.set_ylim((-0.1, 1.1))
	
		# Experimental data
		if len(XAS_EXP) > 0:
			xas_plt.scatter(XAS_EXP[:,0], normalize(XAS_EXP[:,1]), color='black')
		
		# Transitions
		xas_plt.bar(XAS.values[:,0], XAS.values[:,1], width=xsbw, align='center', color='blue')
		
		# Voigt
		xas_plt.plot(XAS_E, 0.70*normalize(VOIGT_XAS), aa=True, color='blue', lw=1)
	
	###### CORE ######
	# Experimental data
	CORE_EXP = import_exp2(C.pm.coreopt.exp)

	# Init
	CORE = C.A[C.A.spec == 'C']
	CORE_3 = CORE[CORE.orb.isin([3])].loc[:, ['E', 'I']]
	CORE_1 = CORE[CORE.orb.isin([1])].loc[:, ['E', 'I']]

	# Energy range and background
	if len(CORE_EXP) > 0:
		CORE_E = np.linspace(min(CORE_EXP[:,0]), max(CORE_EXP[:,0]), vnop)
		if C.pm.coreopt.bg:
			CORE_BG = CORE_EXP[-1, 1]*normalize(np.interp(CORE_E, CORE_EXP[:,0], np.cumsum(CORE_EXP[:,1])))
		else:
			CORE_BG = np.zeros(CORE_E.shape)
	else:
		CORE_E = np.linspace(CORE.E.min() - 10, CORE.E.max() + 10, vnop)
		CORE_BG = np.zeros(CORE_E.shape) # TODO: background do calculo
	
	# Voigts
	VOIGT_CORE_3 = voigt(CORE_3.values, CORE_E, C.pm.coreopt.wG, C.pm.coreopt.wL, C.pm.coreopt.aL)
	VOIGT_CORE_1 = voigt(CORE_1.values, CORE_E, C.pm.coreopt.wG, C.pm.coreopt.wL, C.pm.coreopt.aL)
	VOIGT_CORE_T = VOIGT_CORE_3 + VOIGT_CORE_1
	VOIGT_CORE_3 = VOIGT_CORE_3/max(VOIGT_CORE_T)
	VOIGT_CORE_1 = VOIGT_CORE_1/max(VOIGT_CORE_T)
	VOIGT_CORE_T = normalize(VOIGT_CORE_T) + CORE_BG # TODO: Renormalizar

	# Export
	if C.pm.export:
		export_pd('out/'+cname+'_CORE_BARS.OUT', pd.concat([pd.DataFrame(CORE_3.values, columns=['E32', 'I32']), pd.DataFrame(CORE_1.values, columns=['E12', 'I12'])], axis=1))
		export_np('out/'+cname+'_CORE_VOIGT.OUT', np.vstack((CORE_E, VOIGT_CORE_3, VOIGT_CORE_1, CORE_BG, VOIGT_CORE_T)).T, 'E,I32,I12,Ibg,Itot')

	# PLOT
	if C.pm.plot:
		# Figure initialization
		cs_plt = f.add_subplot(2, 2, 2)
		cs_plt.set_title('Core Level')
		cs_plt.set_xlabel('Binding Energy (eV)')
		cs_plt.set_ylabel('Norm. Intensity')
		cs_plt.yaxis.set_major_locator(plt.NullLocator())
		cs_plt.set_xlim((min(CORE_E), max(CORE_E)))
		cs_plt.set_ylim((-0.1, 1.3))
		cs_plt.invert_xaxis()
		
		# Experimental data
		if len(CORE_EXP) > 0:
			cs_plt.scatter(CORE_EXP[:,0], normalize(CORE_EXP[:,1]), color='black')
		
		# Transitions
		cs_plt.bar(CORE_3.values[:,0], CORE_3.values[:,1], width=csbw, align='center', color='red')
		cs_plt.bar(CORE_1.values[:,0], CORE_1.values[:,1], width=csbw, align='center', color='orange')
		
		# Voigts
		cs_plt.plot(CORE_E, VOIGT_CORE_3, aa=True, color='red', lw=1)
		cs_plt.plot(CORE_E, VOIGT_CORE_1, aa=True, color='orange', lw=1)
		cs_plt.plot(CORE_E, VOIGT_CORE_T, aa=True, color='black', lw=2)
		cs_plt.plot(CORE_E, CORE_BG, aa=True, color='black', lw=0.5, linestyle='dashed')



	###### SPECTRAL WEIGHT ######
	# Init
	SW = C.A[(C.A.spec == 'R') | (C.A.spec == 'A')]
	SW_D_UP = SW[((SW.char == 'd') | (SW.char == 'c')) & (SW.orb <= 5)].loc[:, ['E', 'I']].groupby('E', as_index=False).sum()
	SW_D_DN = SW.loc[((SW.char == 'd') | (SW.char == 'c')) & (SW.orb  > 5), ['E', 'I']].groupby('E', as_index=False).sum()
	SW_P_UP = SW[(SW.char == 'p') & (SW.orb <= 5)].loc[:, ['E', 'I']].groupby('E', as_index=False).sum()
	SW_P_DN = SW[(SW.char == 'p') & (SW.orb  > 5)].loc[:, ['E', 'I']].groupby('E', as_index=False).sum()
	SW_OCT_UP = OCT.copy()
	SW_OCT_DN = OCT.copy()
	
	# Signal
	SW_OCT_UP.loc[:, 'I'] = SW_OCT_UP.I.values/2
	SW_OCT_DN.loc[:, 'I'] = SW_OCT_DN.I.values/2

	# Energy range
	SW_E = np.linspace(-12, 12, vnop)
	
	# Voigts
	V_SW_D_UP = voigt(SW_D_UP.values, SW_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_D_DN = voigt(SW_D_DN.values, SW_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_P_UP = voigt(SW_P_UP.values, SW_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_P_DN = voigt(SW_P_DN.values, SW_E, C.pm.vbopt.wG, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_OCT_UP  = voigt(SW_OCT_UP.values, SW_E, C.pm.vbopt.wG_oct, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_OCT_DN  = voigt(SW_OCT_DN.values, SW_E, C.pm.vbopt.wG_oct, C.pm.vbopt.wL, C.pm.vbopt.aL)
	V_SW_T_UP = V_SW_D_UP + V_SW_P_UP + V_SW_OCT_UP
	V_SW_T_DN = V_SW_D_DN + V_SW_P_DN + V_SW_OCT_DN

	# Signal again
	SW_D_DN.loc[:, 'I'] = -SW_D_DN.I.values
	SW_P_DN.loc[:, 'I'] = -SW_P_DN.I.values
	SW_OCT_DN.loc[:, 'I'] = -SW_OCT_DN.I.values
	
	# Export
	if C.pm.export:
		export_pd('out/'+cname+'_SW_BARS.OUT', pd.concat([pd.DataFrame(SW_D_UP.values, columns=['Edup', 'Idup']),
		                                       pd.DataFrame(SW_D_DN.values, columns=['Eddn', 'Iddn']),
											   pd.DataFrame(SW_P_UP.values, columns=['Epup', 'Ipup']),
											   pd.DataFrame(SW_P_DN.values, columns=['Epdn', 'Ipdn']),
											   pd.DataFrame(SW_OCT_UP.values, columns=['Eoctup', 'Ioctup']),
											   pd.DataFrame(SW_OCT_DN.values, columns=['Eoctdn', 'Ioctdn'])], axis=1))
		export_np('out/'+cname+'_SW_VOIGT.OUT', np.vstack((SW_E, V_SW_D_UP, -V_SW_D_DN, V_SW_P_UP, -V_SW_P_DN, V_SW_OCT_UP, -V_SW_OCT_DN, V_SW_T_UP, -V_SW_T_DN)).T, 'E,dup,ddn,pup,pdn,octup,octdn,totup,totdn')

	# PLOT
	if C.pm.plot:
		# Figure initialization
		sw_plt = f.add_subplot(2, 2, 1)
		sw_plt.set_title('Spectral Weight')
		sw_plt.set_xlabel('Energy (eV)')
		sw_plt.set_ylabel('Intensity')
		sw_plt.yaxis.set_major_locator(plt.NullLocator())
		sw_plt.set_xlim((-12.0, 12.0))
		sw_plt.set_ylim((-1.1*max(V_SW_T_DN), 1.1*max(V_SW_T_UP)))
		
		# Transitions
		sw_plt.bar(SW_D_UP.values[:,0], SW_D_UP.values[:,1], width=swbw, align='center', color='red')
		sw_plt.bar(SW_D_DN.values[:,0], SW_D_DN.values[:,1], width=swbw, align='center', color='red')
		sw_plt.bar(SW_P_UP.values[:,0], SW_P_UP.values[:,1], width=swbw, align='center', color='blue')
		sw_plt.bar(SW_P_DN.values[:,0], SW_P_DN.values[:,1], width=swbw, align='center', color='blue')
		sw_plt.bar(SW_OCT_UP.values[:,0], C.pm.vbopt.octmult*SW_OCT_UP.values[:,1], width=swbw, align='center', color='cyan')
		sw_plt.bar(SW_OCT_DN.values[:,0], C.pm.vbopt.octmult*SW_OCT_DN.values[:,1], width=swbw, align='center', color='cyan')
		
		# Voigts
		sw_plt.plot(SW_E, V_SW_D_UP, aa=True, color='red', lw=1)
		sw_plt.plot(SW_E, -V_SW_D_DN, aa=True, color='red', lw=1)
		sw_plt.plot(SW_E, V_SW_P_UP, aa=True, color='blue', lw=1)
		sw_plt.plot(SW_E, -V_SW_P_DN, aa=True, color='blue', lw=1)
		sw_plt.plot(SW_E, V_SW_OCT_UP, aa=True, color='cyan', lw=1)
		sw_plt.plot(SW_E, -V_SW_OCT_DN, aa=True, color='cyan', lw=1)
		sw_plt.plot(SW_E, V_SW_T_UP, aa=True, color='black', lw=2)
		sw_plt.plot(SW_E, -V_SW_T_DN, aa=True, color='black', lw=2)
		
		# Guidelines
		sw_plt.plot([0, 0], [-100, 100], color='black', aa=True)
		sw_plt.plot([-100, 100], [0, 0], color='black', aa=True)
		
		#cursor = FollowDotCursor(sw_plt, OCT.values[:,0], OCT.values[:,1]/2)
	
	write_output(C)
	
	if C.pm.plot:
		f.tight_layout()
		f.show()

	if C.pm.dev:
		keyboard()
		
	return C


