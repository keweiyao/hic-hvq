duke heavy quark project
==========================

An ebe package for local tests.

Requirements:

  cmake2.8+
  
  FORTRAN

  c++11, boost lib and hdf5 lib

  python3, h5py lib and Cython installed

First TRENTo generate the xy-IC for both hard (T_AA) and soft (T_R) event.
The pT distribution is reweighted as FONLL+EPS09nlo.

First run medium evolution:

   + TRENTo + freestream + Vishnew + frzout + UrQMD

With the trento IC and hydro evolution history, run heavy quark evolution:

   + TRENTo/FONLL + linear boltzmann / Langevin + frag&Recomb

To build packages into run/

.. code::

  ./makepkg run

Then go to run/ and run event

.. code::
  
  cd run
  ./job-wrapper inputfile

+ Initial Condition File: initial.hdf
+ Hydro history and hypersurface: JetData.h5 and surface.dat
+ Final state heavy quark data: hvq-final.dat (OSCAR1997A)
+ Final state heavy meson data: h-meson-final.dat (OSCAR1997A)
+ Final state light hadron data: particles_out.dat


