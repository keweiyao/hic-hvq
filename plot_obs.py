#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from scipy.interpolate import interp1d

x0, y0 = np.loadtxt(sys.argv[1]).T
ppref = interp1d(x0, y0)
pptot = (y0*x0).sum()*(x0[1]-x0[0])
x0, y1 = np.loadtxt(sys.argv[2]).T
AAtot = (y1*x0).sum()*(x0[1]-x0[0])
AAinit = interp1d(x0, y1)
ratio = AAtot/pptot
print(ratio)


flist = [h5py.File(fname, 'r') for fname in sys.argv[3:] ]
label = [r'LBT $C_{22}$', r'LBT $C_{22, 23}$', r'LBT $C_{22, 23, 32}$', 'LGV w/ Ein'\
		, 'LGV w/o Ein']

def get_v2_light(pmu):
	N = 15
	pTbins = np.linspace(0, 4, N+1)
	pT = np.sqrt(pmu[1]**2 + pmu[2]**2)
	q2 = (pmu[1]**2 - pmu[2]**2)/(pmu[1]**2 + pmu[2]**2)
	v2 = np.zeros(N)
	for i in range(N):
		index = (pT > pTbins[i]) & (pT < pTbins[i+1])
		v2[i] = np.mean(q2[index])
	x = (pTbins[1:] + pTbins[:-1])*0.5
	return x, v2
	

def get_Raa(pmu, pT0, AAtot):
	pT = np.sqrt(pmu[1]**2 + pmu[2]**2)
	weight = pT0*AAinit(pT0)
	H, be = np.histogram(pT, bins=25, range=[0.1, 50], normed=True, weights=weight)
	x = (be[1:] + be[:-1])*0.5
	ppy = ppref(x)
	return x, H/x/ppy*AAtot

def get_v2(pmu, pT0):
	N = 15
	#pTbins = np.concatenate([np.linspace(0, 10, n1+1)[:-1], np.linspace(10, 30., n2+1)])
	pTbins = np.exp(np.linspace(0., np.log(20.), N+1))-1.
	pT = np.sqrt(pmu[1]**2 + pmu[2]**2)
	q2 = (pmu[1]**2 - pmu[2]**2)/(pmu[1]**2 + pmu[2]**2)
	weight = pT0*AAinit(pT0)
	v2 = np.zeros(N)
	for i in range(N):
		index = (pT > pTbins[i]) & (pT < pTbins[i+1])
		v2[i] = np.average(q2[index], weights=weight[index])
	x = (pTbins[1:] + pTbins[:-1])*0.5
	return x, v2


#pmu = np.loadtxt('lighthadrons.dat').T
#x0, y0 = get_v2_light(pmu)
for i in range(371): #129
	if i%10 != 0:
		continue
	plt.clf()
	
	#plt.subplot(2,1,2)
	#plt.plot(x0, y0, '-', linewidth=3., color='grey', label = 'Light hadrons (hydro)')
	
	for j, f in enumerate(flist):
		pmu = f['p-%d'%i].value.T
		pT0 = f['init_pT'].value
		
		plt.subplot(2,1,1)
		x, y = get_Raa(pmu, pT0, AAtot)
		plt.plot(x, y, label = label[j])
		
		plt.subplot(2,1,2)
		x, y = get_v2(pmu, pT0)
		plt.plot(x, y, 'o-', label = label[j])
		
	plt.subplot(2,1,1)
	plt.plot([0,50],[1,1])
	plt.axis([0,50,0,1.4])
	plt.legend()
	plt.xlabel(r'$p_T$ [GeV]')
	plt.ylabel(r'$R_{AA}$')
	
	plt.subplot(2,1,2)
	plt.plot([0,50],[0,0])
	plt.axis([0, 20, -0.1, 0.1])
	plt.legend()
	plt.xlabel(r'$p_T$ [GeV]')
	plt.ylabel(r'$(p_x^2 - p_y^2)/(p_x^2 + p_y^2)$')
	plt.suptitle("Charm quark: LHC 2760 GeV, 30-50%-averaged TRENTo-IC + Vishnew hydro")
	plt.pause(0.1)
	
	
plt.show()
