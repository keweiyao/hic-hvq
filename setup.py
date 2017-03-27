from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob
import numpy


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


#-------------HqEvo Module------------------
# LBT
filelist_LBT = list(set(['./HQ-Evo/cython/HqEvo.pyx'] + glob('./HQ-Evo/src/*.cpp'))-set(["./HQ-Evo/src/main.cpp"]))
extensions = [
        Extension('HqEvo', filelist_LBT, language="c++", extra_compile_args=["-std=c++11", '-march=native'],libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"])]
setup(ext_modules=cythonize(extensions))

# LGV
filelist_LGV = list(set(['./HQ-Evo/cython/HqLGV.pyx'] + glob('./HQ-Evo/src/*.cpp'))-set(["./HQ-Evo/src/main.cpp"]))
extensions = [
        Extension('HqLGV', filelist_LGV, language="c++", extra_compile_args=["-std=c++11", '-march=native'],libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"])
]
setup(ext_modules=cythonize(extensions))

#-----------Medium Reader------------------
setup(
    ext_modules = cythonize([
    Extension("medium", ["./Medium-Reader/medium.pyx"], language="c++", 
              libraries=["m"])
    ]),
     include_dirs=[numpy.get_include()]
)

#-------------Event Module compiled at last------------------
filelist1 = ["./Event/cython/event.pyx", "./Event/src/utility.cpp"]
extensions = [
	Extension(
		'event',
		filelist1,
		language="c++",
		extra_compile_args=["-std=c++11"],
		libraries=["m"])
]

setup(
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()]
        )


"""
    ext_modules=cythonize(extensions,
			compiler_directives={ 'c_string_type':str, 'c_string_encoding':ascii})
"""
