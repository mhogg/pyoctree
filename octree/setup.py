
# python setup.py build_ext --inplace --compiler=msvc

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("octree",
               sources=["octree.pyx","cOctree.cpp"],
               include_dirs=[numpy.get_include()],
               extra_compile_args=['/openmp'],
               #extra_link_args=['/openmp'],
               language="c++"),]

setup(
    name = "Octree structure",
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

