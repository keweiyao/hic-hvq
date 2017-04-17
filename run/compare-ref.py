#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

scale = [np.pi*1e-9, 1.]
flist = sys.argv[1:]
for i, f in enumerate(flist):
	x, y  = np.loadtxt(f).T
	plt.plot(x, y*scale[i], label=os.path.basename(f))
	plt.semilogy()
plt.legend()
plt.show()

