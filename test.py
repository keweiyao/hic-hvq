#!/usr/bin/env python3

import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
# This is "main" script of LBT-LGV code
# First import the heavy quark event module
import event
# To declare an event module, we need a bunch of options
# They are split into two categories
# Event options: 
# 1. Medium options: options are organized in a dictionary by { 'keyward' : value }

# For example, a medium-option dictionary config a simulation with a static medium
static_config = {	'type'	  : 'static',	
					'static_dt' : 1.  }
# The following medium information is to be used later
box_info = {'Temp'  : 0.3, 
			'Vx'	: 0.0, 
			'Vy'	: 0.0, 
			'Vz'	: 0.0   }

# Another example, a medium-option dictionary config a dynamical medium
# with the path to the hydro-history file 
hydro_history_filepath = sys.argv[1]
dynamic_config = {  'type'	  : 'dynamic', 
					'hydrofile' : hydro_history_filepath	}

# Physics option
# 1. An example for linear-Boltzmann evolution
LBT_config = {  'physics'   : 'LBT',
                                '2->2'    : True,
                                '2->3'    : False,
                                '3->2'    : False,
                                'Nf'        : 3,
                                'mass'    : 1.3 }  
                                
# 1. An example for Langevin evolution
LGV_config = {  'physics'   : 'LGV',
                                'dt_lrf'        : 0.02,
                                'elastic'   : True,
                                'Einstein'  : False,
                                'Nf'        : 3,
                                'mass'    : 1.3 } 

# Initialization option
# Initizlize HQ in a static box [-L, L]^3
box_init = {	'type'  : 'box',
				'L'	 : 10.,
				'pmax'  : 10.   }

TAA = np.loadtxt(sys.argv[2]).T**2
realistic_init =  { 'type'		  : 'A+B',
					'sample power'  : 1.,
					'pTmin'		 : 0.1,
					'pTmax'		 : 70.,
					'ymin'		  : -1.,
					'ymax'		  : 1.,
					'TAB'		   : TAA,
					'dxy'		   : 0.1   }
		
e1 = event.event(   medium_flags=dynamic_config , 
					physics_flags=LGV_config   )

e1.initialize_HQ(   NQ=200000,
					init_flags=realistic_init   )

# Run Model
f = h5py.File("lgv-no-ein.hdf5", 'w')
Init_pT = e1.Init_pT()
f.create_dataset('init_pT', data=Init_pT)

for i in range(500):
	print("t = %1.2f [fm/c]"%e1.sys_time() )
	status = e1.perform_hydro_step(StaticPropertyDictionary=box_info)
	
	if i%10 == 0:
		dsp, dsx = e1.HQ_hist()
		f.create_dataset('p-%d'%i, data=dsp)
	
	if not status:
		break
f.close()
