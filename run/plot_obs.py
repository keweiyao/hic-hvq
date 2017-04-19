#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import glob
import fortranformat as ff
from scipy.interpolate import interp1d
import HepPlot
from matplotlib.offsetbox import AnchoredText

# pp->c reference
ip = np.linspace(0.0,100,10000)
dip = ip[1]-ip[0]
x0, y0 = np.loadtxt(sys.argv[1]).T[:,2:]
pp_dXc_dpT = interp1d(x0, y0, fill_value='extrapolate')
pp_Xc = dip*np.sum(ip*pp_dXc_dpT(ip))
# PbPb->c reference
x1, y1 = np.loadtxt(sys.argv[2]).T[:,2:]
AA_dXc_dpT = interp1d(x1, y1, fill_value='extrapolate')
AA_Xc = dip*np.sum(ip*AA_dXc_dpT(ip))
# pp->D reference
x2, y2 = np.loadtxt(sys.argv[3]).T
pp_dXD_dpT = interp1d(x2, y2, fill_value='extrapolate')

label = [r'Duke $c$', r'Duke, $D$', r'Marlene, $c$', r'Marlene, $D$']

def unpack_file(filename):
	with h5py.File(filename, 'r') as f:
		Qpv = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'quark' and k.split('-')[1] == 'pv'], axis=0)
		QpT0 = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'quark' and k.split('-')[1] == 'pT0'], axis=0)
		Dpv = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'meson' and k.split('-')[1] == 'pv'], axis=0)
		DpT0 = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'meson' and k.split('-')[1] == 'pT0'], axis=0)
		
		Qw = QpT0*AA_dXc_dpT(QpT0)
		Dw = DpT0*AA_dXc_dpT(DpT0)
		return Qpv.T, QpT0, QpT0*AA_dXc_dpT(QpT0), \
			   Dpv.T, DpT0, DpT0*AA_dXc_dpT(DpT0), \
			   np.sum(Dw)/np.sum(Qw)

def get_Raa(pmu, weight, AAtot, ref):
	pT = np.sqrt(pmu[0]**2 + pmu[1]**2)
	H, be = np.histogram(pT, bins=100, range=[0.0, 50], normed=True, weights=weight)
	x = (be[1:] + be[:-1])*0.5
	return x, H/x/ref(x)*AAtot

def get_v2(pmu, weight):
	N = 11
	pTbins = np.exp(np.linspace(0., np.log(40.), N+1))-1.
	pT = np.sqrt(pmu[0]**2 + pmu[1]**2)
	q2 = (pmu[0]**2 - pmu[1]**2)/(pmu[0]**2 + pmu[1]**2)
	v2 = np.zeros(N)
	for i in range(N):
		index = (pT > pTbins[i]) & (pT < pTbins[i+1])
		v2[i] = np.average(q2[index], weights=weight[index])
	x = (pTbins[1:] + pTbins[:-1])*0.5
	return x, v2


mx, myD, myc = np.loadtxt(sys.argv[4]).T
for i, filename in enumerate(sys.argv[5:]):
	Qp, QpT0, Qw, Dp, DpT0, Dw, ratio_DoverC = unpack_file(filename)
	plt.subplot(2,1,1)
	x, y = get_Raa(Qp, Qw, AA_Xc, pp_dXc_dpT)
	plt.plot(x, y, 'b-', label = label[0])
	plt.plot(mx, myc, 'r-', label = label[2])

	plt.subplot(2,1,2)
	print (ratio_DoverC)
	x, y = get_Raa(Dp, Dw, AA_Xc*ratio_DoverC, pp_dXD_dpT)
	plt.plot(x, y, 'b-', label = label[1])
	plt.plot(mx, myD, 'r-', label = label[3])
"""
	plt.subplot(3,1,3)
	x, y = get_v2(Qp, Qw)
	plt.plot(x, y, 'y-' , label = label[0])
	x, y = get_v2(Dp, Dw)
	plt.plot(x, y, 'go--', label = label[1])
"""


plt.subplot(2,1,1)
plt.axis([0,40,0,1.2])
plt.xlabel(r'$p_T$ [GeV]')
plt.ylabel(r'$R_{AA}$')
plt.xticks([])

plt.legend(loc='best', framealpha=0.)

ax = plt.subplot(2,1,2)
plt.axis([0,40,0,1.1])

plt.xlabel(r'$p_T$ [GeV]')
plt.ylabel(r'$R_{AA}$')
color = ['g','b','y']
alpha = [0.2, 0.3, 0.4, 0.5]

pi =  HepPlot.plot_info('/LHC-2760-PbPb-D-Raa-central/Table17.yaml')
for i, p in enumerate(pi.plist):
	ax.errorbar(p.x, p.y,[[-p.stat[0]], [p.stat[1]]], fmt='gD', markersize=5,
				label='ALICE 0-20%' if i==0 else '')
	for it, c, a in zip(p.syslist, color, alpha):
		ax.fill_between([p.xl, p.xh], [p.y+it[0], p.y+it[0]],  [p.y+it[1], p.y+it[1]], alpha=a, color=c)
#ax.set_xlabel(pi.xlabel)
#ax.set_ylabel(pi.ylabel)
#box = AnchoredText(pi.description, loc=1, frameon=False)
#ax.add_artist(box)
#pi.recommand_scale(ax)

plt.legend(loc='best', framealpha=0.)
"""
plt.subplot(3,1,3)
plt.plot([0,40],[0,0])
plt.axis([0,40,-0.05, 0.25])
plt.xlabel(r'$p_T$ [GeV]')
plt.ylabel(r'$v_2$')
"""

plt.suptitle("LHC 2760 GeV, Pb+Pb-->D+X\n"+r"Elastic $K=1$, self consistent $m_D$")
plt.subplots_adjust(hspace=0.)
	
plt.show()
