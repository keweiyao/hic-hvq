import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt

f = h5py.File(sys.argv[1], 'r')

sd = np.zeros([261, 261])
taa = np.zeros([261, 261])

for k in f.keys():
	if k[0] == 'e':
		sd = sd + f[k].value
	if k[0] == 'T':
		taa = taa + f[k].value
f.close()

f = h5py.File('initial.hdf' ,'w')
sd /= 1000
taa /= 1000

f.create_dataset('TAB_0', data=taa)
ds = f.create_dataset('event_0', data=sd)
ds.attrs.create('dx', data=0.1)
f.close()
plt.imshow(sd.T)
plt.show()
plt.imshow(taa.T)
plt.show()
