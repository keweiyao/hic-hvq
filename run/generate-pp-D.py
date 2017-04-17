#!/usr/bin/env python3

import numpy as np
import fortranformat as ff
import sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import glob
import h5py

# pp->c reference
ip = np.linspace(0.0,100,10000)
dip = ip[1]-ip[0]
x0, y0 = np.loadtxt(sys.argv[1]).T[:,2:]
ppref = interp1d(x0, y0, fill_value='extrapolate')
pptot = dip*np.sum(ip*ppref(ip))

def unpack_file(filename):
	with h5py.File(filename, 'r') as f:
		Qpv = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'quark' and k.split('-')[1] == 'pv'], axis=0)
		QpT0 = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'quark' and k.split('-')[1] == 'pT0'], axis=0)
		Dpv = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'meson' and k.split('-')[1] == 'pv'], axis=0)
		DpT0 = np.concatenate([f[k].value for k in f.keys() if k.split('-')[0] == 'meson' and k.split('-')[1] == 'pT0'], axis=0)

		Qw = QpT0*ppref(QpT0)
		Dw = DpT0*ppref(DpT0)
		print (len(Qw)/len(Dw))
		return Qpv.T, QpT0, QpT0*ppref(QpT0), \
			   Dpv.T, DpT0, DpT0*ppref(DpT0), \
			   np.sum(Dw)/np.sum(Qw)

plt.plot(x0, y0, 'r-')

Qp, QpT0, Qw, Dp, DpT0, Dw, ratio = unpack_file(sys.argv[2])

DpT = np.sqrt(Dp[0]**2+ Dp[1]**2)
y, x = np.histogram(DpT, bins=1000, range=[0.0, 70], weights=Dw, normed=True)
x = 0.5*(x[1:] + x[:-1])
yD = y/x*pptot*ratio
plt.plot(x, yD, 'b-')
with open('pp-D-spectra.dat','w') as f:
	for xx, yy in zip(x,yD):
		f.write("{} {}\n".format(xx, yy))

plt.semilogy()
plt.show()



