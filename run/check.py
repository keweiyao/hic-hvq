#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt


def FONLL_spectra(filename, nuA='Pb', nuB='Pb', sqrts=2760, mass=1.3):
	f = h5py.File(filename, 'r')
	ds = f['{}-{}-{}'.format(sqrts, nuA, nuB)]
	b1 = ds['b1'].value
	b2 = ds['b2'].value
	normy = np.linspace(-1,1,13)
	pT = ds['pT'].value
	dsigma = ds['dsigma_bdb_dpT2_dy'].value*pT
	f.close()

	dnormy = normy[1]-normy[0]
	pTmin, pTmax, pTmid = pT[0], pT[-1], pT[29]
	dpT1, dpT2 = pT[1]-pT[0], pT[-1] - pT[-2]
	bmin, bmax, db = b1[0], b1[-1], b1[1]-b1[0]
	def interp(b1i, b2i, yi, pTi):
		mT = np.sqrt(pTi**2 + mass**2)
		normyi = yi/np.arccosh(sqrts/2./mT)
		index_b1 = int((b1i-bmin)/db)
		index_b2 = int((b2i-bmin)/db)
		index_normy = int((normyi+1.)/dnormy)
		A1 = np.max([np.min([15, index_b1]), 0])
		A2 = np.max([np.min([15, index_b2]), 0])
		A3 = np.max([np.min([11, index_normy]), 0])
		rx1 = (b1i-b1[A1])/db
		rx2 = (b2i-b2[A2])/db
		rx3 = (normyi-normy[A3])/dnormy
		if pTi < pTmid:
			index_pT = int((pTi - pTmin)/dpT1)
			A4 = np.max([np.min([28, index_pT]), 0])
			rx4 = (pTi-pT[A4])/dpT1
		else:
			index_pT = int((pTi - pTmid)/dpT2) + 29
			A4 = np.max([np.min([58, index_pT]), 29])
			rx4 = (pTi-pT[A4])/dpT2

		resultA = 0.0
		ra1 = [1.-rx1, rx1]
		ra2 = [1.-rx2, rx2]
		ra3 = [1.-rx3, rx3]
		ra4 = [1.-rx4, rx4]
		for i1 in range(2):
			for i2 in range(2):
				for i3 in range(2):
					for i4 in range(2):
						resultA += dsigma[A1+i1, A2+i2, A3+i3, A4+i4] \
									*ra1[i1]*ra2[i2]*ra3[i3]*ra4[i4]
		return resultA
	interp = np.vectorize(interp)
	return interp

def p4_to_pT_and_y(px, py, pz, E):
	return np.sqrt(px**2+py**2), 0.5*np.log((E+pz)/(E-pz))
def p4_to_pT_phi_and_y(px, py, pz, E):
	return np.sqrt(px**2+py**2), np.arctan2(py,px), 	\
		0.5*np.log((E+pz)/(E-pz))
def read_text_file(filename):
	"""
	Read a text file into a nested list of bytes objects,
	skipping comment lines (#).

	"""
	with open(filename, 'rb') as f:
		return [l.split() for l in f if not l.startswith(b'#')]

mass = 1.3
run_id = 1
FONLL = FONLL_spectra("./hic-hvq/share/hvq/lhc/lhc-PbPb-b1-b2-y-pT.hdf5",
							nuA='Pb', nuB='Pb', sqrts=2760, mass=mass)
fmeson = h5py.File('./HeavyFlavorResult.hdf5', 'r')
pT0, y0 = p4_to_pT_and_y(*fmeson['meson-p0-{}'.format(run_id)].value.T)
pT, y = p4_to_pT_and_y(
						*fmeson['meson-p-{}'.format(run_id)].value.T)
s1, s2 = fmeson['meson-s1-s2-{}'.format(run_id)].value.T
pid = fmeson['meson-pid-{}'.format(run_id)].value
cut = (pT > 3.0) & (pT < 4.0)
weight = FONLL(s1,s2,y0, pT0)
plt.hist(y[cut], weights=weight[cut], color='r', bins=7, range=[-3,3], histtype='step')
#plt.hist(pT0, color='r', bins=101, range=[0,50], histtype='step')

ID, charge, fmass, px, py, pz, y, eta, pT0, y0, s1, s2 = (
			np.array(col, dtype=dtype) for (col, dtype) in
			zip(
				zip(*read_text_file('particles_out_{}.dat'.format(run_id))),
				(2*[int] + 10*[float])
			)
		)

pT = np.sqrt(px**2 + py**2)
weight = FONLL(s1,s2, y0, pT0)
cut = np.array([iid in [411, 421, 10411] for iid in ID]) & (pT > 3.0) & (pT < 4.0)
print(np.count_nonzero(cut))
plt.hist(y[cut], weights=weight[cut], color='b', bins=7, range=[-3,3], histtype='step')
#plt.hist(pT0[cut], color='b', bins=101, range=[0,50], histtype='step')
plt.show()

