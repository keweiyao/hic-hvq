from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob
import numpy
import os

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


#-------------HqEvo Module------------------
#includes=['']
libs=[path for path in os.environ['LD_LIBRARY_PATH'].split(':') if path]
#-------------HqEvo Module------------------
fileLBT = [     'cython/HqEvo.pyx',
                        'src/matrix_elements.cpp',
                        'src/utility.cpp',
                        'src/Xsection.cpp',
                        'src/sample_methods.cpp',
                        'src/rates.cpp']
fileLBT = ['HQ-Evo/'+it for it in fileLBT]
fileLGV = [     'cython/HqLGV.pyx',
                        'src/Langevin.cpp',
                        'src/matrix_elements.cpp',
                        'src/qhat_matrix_elements.cpp',
                        'src/utility.cpp',
                        'src/qhat_Xsection.cpp',
                        'src/qhat.cpp',
                        'src/TLorentz.cpp']
fileLGV = ['HQ-Evo/'+it for it in fileLGV]
modules = [
        Extension('HqEvo',
                         sources=fileLBT,
                         language="c++",
#                        include_dirs=includes,
                         library_dirs=libs,
                         extra_compile_args=["-std=c++11", '-fPIC'],
                         libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"]),
                Extension('HqLGV',
                         sources=fileLGV,
                         language="c++",
#                        include_dirs=includes,
                         library_dirs=libs,
                         extra_compile_args=["-std=c++11", '-fPIC'],
                         libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"]),
]

data = ('/share/hvq/tables/', 
	[fn for fn in glob('./HQ-Evo/tables/*.hdf5')] )

setup(
        ext_modules=cythonize(modules),
	data_files=[data]
)


###-----------Medium Reader------------------
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
