import numpy as np
from math import pi, sqrt, erfc, exp
from scipy.special import wofz
import matplotlib.pyplot as plt
import yaml
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=999999)

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


inp1 = 'SrRuO3'
inp2 = 'SrTiO3'

C1_WEIGHT = 0.5
C2_WEIGHT = 0.5

pm1 = obj_dic(yaml.safe_load(open(inp1+'.in')))
pm2 = obj_dic(yaml.safe_load(open(inp2+'.in')))

#### CORE ####

CS_E_MIN = 445
CS_E_MAX = 500
CS_E_NOP =  10*(CS_E_MAX - CS_E_MIN) + 1# number of points

CS_1_csec = 1.0
CS_2_csec = 0.36

CS_1_BARS = np.loadtxt('out/'+inp1+'_CORE_BARS.OUT', unpack=True, skiprows=1, delimiter=',').T
CS_2_BARS = np.loadtxt('out/'+inp2+'_CORE_BARS.OUT', unpack=True, skiprows=1, delimiter=',').T  

CS_1_BARS[:, 1] = CS_1_csec*C1_WEIGHT*CS_1_BARS[:, 1]
CS_1_BARS[:, 3] = CS_1_csec*C1_WEIGHT*CS_1_BARS[:, 3]
CS_2_BARS[:, 1] = CS_2_csec*C2_WEIGHT*CS_2_BARS[:, 1]
CS_2_BARS[:, 3] = CS_2_csec*C2_WEIGHT*CS_2_BARS[:, 3]

CS_E = np.linspace(CS_E_MIN, CS_E_MAX, CS_E_NOP)
CS_1_32_VOIGT = voigt(CS_1_BARS[:, 0:2], CS_E, pm1.coreopt.wG, pm1.coreopt.wL, pm1.coreopt.aL)
CS_1_12_VOIGT = voigt(CS_1_BARS[:, 2:4], CS_E, pm1.coreopt.wG, pm1.coreopt.wL, pm1.coreopt.aL)
CS_2_32_VOIGT = voigt(CS_2_BARS[:, 0:2], CS_E, pm2.coreopt.wG, pm2.coreopt.wL, pm2.coreopt.aL)
CS_2_12_VOIGT = voigt(CS_2_BARS[:, 2:4], CS_E, pm2.coreopt.wG, pm2.coreopt.wL, pm2.coreopt.aL)
CS_T_VOIGT = CS_1_32_VOIGT + CS_1_12_VOIGT + CS_2_32_VOIGT + CS_2_12_VOIGT
CS_T_VOIGT_MAX = max(CS_T_VOIGT)
CS_T_VOIGT = CS_T_VOIGT/CS_T_VOIGT_MAX
np.savetxt('out/' + inp1+inp2 + '_CORE_VOIGT.OUT', np.vstack((CS_E, CS_T_VOIGT)).T, delimiter=' ', fmt='%.6f')

exp = np.loadtxt('exp/STROcore.dat', unpack=True)
exp = exp.T
plt.scatter(exp[:,0], exp[:,1]/max(exp[:,1]))
plt.plot(CS_E, CS_T_VOIGT, color='black')
plt.plot(CS_E, (CS_1_32_VOIGT+CS_1_12_VOIGT)/CS_T_VOIGT_MAX, color='red')
plt.plot(CS_E, (CS_2_32_VOIGT+CS_2_12_VOIGT)/CS_T_VOIGT_MAX, color='green')
plt.gca().invert_xaxis()
plt.show()

#### VB ####

VB_E_MIN = -2
VB_E_MAX = 14
VB_E_NOP =  10*(CS_E_MAX - CS_E_MIN) + 1# number of points #max precisa ser maior que min..

#to confuso com essas csec aqui...
VB_1_csec_d = 1.0 
VB_1_csec_p = 1
VB_1_csec_oct = 1
VB_2_csec_d = 0.04 
VB_2_csec_p = 0.09
VB_2_csec_oct = 0.09

VB_1_BARS = np.loadtxt('out/'+inp1+'_VB_BARS.OUT', unpack=True, skiprows=1, delimiter=' ', usecols=((0,1,))).T
#VB_2_BARS = np.loadtxt('out/'+inp2+'_VB_BARS.OUT', unpack=True, skiprows=1, delimiter=' ').T
#0   1 2   3  4    5
#Ed,Id,Ep,Ip,Eoct,Ioct

VB_1_BARS[:, 1] = VB_1_csec_d*C1_WEIGHT*VB_1_BARS[:, 1]
VB_1_BARS[:, 3] = VB_1_csec_p*C1_WEIGHT*VB_1_BARS[:, 3]
VB_1_BARS[:, 5] = VB_1_csec_oct*C1_WEIGHT*VB_1_BARS[:, 5]
VB_2_BARS[:, 1] = VB_2_csec_d*C2_WEIGHT*VB_2_BARS[:, 1]
VB_2_BARS[:, 3] = VB_2_csec_p*C2_WEIGHT*VB_2_BARS[:, 3]
VB_2_BARS[:, 5] = VB_2_csec_oct*C2_WEIGHT*VB_2_BARS[:, 5]

VB_E = np.linspace(VB_E_MIN, VB_E_MAX, VB_E_NOP)
VB_1_d_VOIGT = voigt(VB_1_BARS[:, 0:2], VB_E, pm1.vbopt.wG, pm1.vbopt.wL, pm1.vbopt.aL)
VB_1_p_VOIGT = voigt(VB_1_BARS[:, 2:4], VB_E, pm1.vbopt.wG, pm1.vbopt.wL, pm1.vbopt.aL)
VB_1_oct_VOIGT = voigt(VB_1_BARS[:, 4:6], VB_E, pm1.vbopt.wG, pm1.vbopt.wL, pm1.vbopt.aL)
VB_2_d_VOIGT = voigt(VB_2_BARS[:, 0:2], VB_E, pm2.vbopt.wG, pm2.vbopt.wL, pm2.vbopt.aL)
VB_2_p_VOIGT = voigt(VB_2_BARS[:, 2:4], VB_E, pm2.vbopt.wG, pm2.vbopt.wL, pm2.vbopt.aL)
VB_2_oct_VOIGT = voigt(VB_2_BARS[:, 4:6], VB_E, pm2.vbopt.wG, pm2.vbopt.wL, pm2.vbopt.aL)

VB_T_VOIGT = VB_1_d_VOIGT + VB_1_p_VOIGT + VB_1_oct_VOIGT + VB_2_d_VOIGT + VB_2_p_VOIGT + VB_2_oct_VOIGT
VB_T_VOIGT_MAX = max(VB_T_VOIGT)
VB_T_VOIGT = VB_T_VOIGT/VB_T_VOIGT_MAX
np.savetxt('out/' + inp1+inp2 + '_VB_VOIGT.OUT', np.vstack((VB_E, VB_T_VOIGT)).T, delimiter=' ', fmt='%.6f')

exp = np.loadtxt('exp/STROVB.dat', unpack=True)
exp = exp.T
plt.scatter(exp[:,0], exp[:,1]/max(exp[:,1]))
plt.plot(VB_E, VB_T_VOIGT, color='black')
plt.plot(VB_E, (VB_1_d_VOIGT+VB_1_p_VOIGT)/VB_T_VOIGT_MAX, color='red')
plt.plot(VB_E, (VB_2_d_VOIGT+VB_2_p_VOIGT)/VB_T_VOIGT_MAX, color='green')
plt.plot(VB_E, (VB_1_oct_VOIGT+VB_2_oct_VOIGT)/VB_T_VOIGT_MAX, color='blue')
plt.gca().invert_xaxis()
plt.show()

#### XAS ####

XAS_E_MIN = 525
XAS_E_MAX = 555
XAS_E_NOP =  10*(CS_E_MAX - CS_E_MIN) + 1# number of points

XAS_1_BARS = np.loadtxt('out/'+inp1+'_XAS_BARS.OUT', unpack=True, skiprows=1, delimiter=',').T
XAS_2_BARS = np.loadtxt('out/'+inp2+'_XAS_BARS.OUT', unpack=True, skiprows=1, delimiter=',').T

XAS_1_BARS[:, 1] = C1_WEIGHT*XAS_1_BARS[:, 1]
XAS_2_BARS[:, 1] = C2_WEIGHT*XAS_2_BARS[:, 1]


XAS_E = np.linspace(XAS_E_MIN, XAS_E_MAX, XAS_E_NOP)
XAS_1_VOIGT = voigt(XAS_1_BARS[:, 0:2], XAS_E, pm1.xasopt.wG, pm1.xasopt.wL, pm1.xasopt.aL)
XAS_2_VOIGT = voigt(XAS_2_BARS[:, 0:2], XAS_E, pm2.xasopt.wG, pm2.xasopt.wL, pm2.xasopt.aL)
XAS_T_VOIGT = XAS_1_VOIGT + XAS_2_VOIGT
XAS_T_VOIGT_MAX = max(XAS_T_VOIGT)
XAS_T_VOIGT = XAS_T_VOIGT/XAS_T_VOIGT_MAX
np.savetxt('out/' + inp1+inp2 + '_XAS_VOIGT.OUT', np.vstack((XAS_E, XAS_T_VOIGT)).T, delimiter=' ', fmt='%.6f')

exp = np.loadtxt('exp/STROXAS.dat', unpack=True)
exp = exp.T
plt.scatter(exp[:,0], exp[:,1]/max(exp[:,1]))
plt.plot(XAS_E, XAS_T_VOIGT, color='black')
plt.plot(XAS_E, XAS_1_VOIGT/CS_T_VOIGT_MAX, color='red')
plt.plot(XAS_E, XAS_2_VOIGT/CS_T_VOIGT_MAX, color='green')
plt.gca().invert_xaxis()
plt.show()
