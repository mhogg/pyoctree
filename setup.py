# -*- coding: utf-8 -*-

# Copyright (C) 2017 Michael Hogg

# This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

from setuptools import setup, find_packages
from setuptools.extension import Extension
from codecs import open
from os import path
import numpy

# get current path
here = path.abspath(path.dirname(__file__))

# function to open the readme file
def readme():
    with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
        return f.read()
        
# find the version
exec(open(path.join('pyoctree','version.py')).read())

try:
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext
    use_cython = False
else:
    use_cython = True

# Supply correct openmp compiler arguments depending on os. Based on SO question:
# http://stackoverflow.com/questions/30985862/how-to-identify-compiler-before-defining-cython-extensions
# However, need to add extra link args in addition to build args, as these are required for gcc.
# From http://stackoverflow.com/questions/16737260/how-to-tell-distutils-to-use-gcc, the distutils --compiler 
# option expects "unix", "msvc", "cygwin", "mingw32", "bcpp", or "emx". Therefore, should add support for 
# all these, but will currently support only mscv, mingw32, and unix. 
# The current code has been tested on Windows 10 and CentOS 6. The CentOS compiler type is "unix", which uses gcc by default. 
BUILD_ARGS = {}
BUILD_ARGS['msvc']    = ['/openmp', '/EHsc', ]
BUILD_ARGS['mingw32'] = ['-fopenmp', ]
BUILD_ARGS['unix']    = ['-fopenmp', ]  # On CentOS, "compiler" variable equals 'unix' for gcc compiler
LINK_ARGS = {}
LINK_ARGS['msvc']     = []
LINK_ARGS['mingw32']  = []
LINK_ARGS['unix']     = ['-lgomp', ]    # gcc requires -lgomp for linking. Otherwise get unrecognised symbol error

# Custom class to add compiler depended build and link arguments
class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        compiler = self.compiler.compiler_type
        for ext in self.extensions:
            if compiler in BUILD_ARGS:
                ext.extra_compile_args = BUILD_ARGS[compiler]
            if compiler in LINK_ARGS:
                ext.extra_link_args    = LINK_ARGS[compiler]
        build_ext.build_extensions(self)

cmdclass    = {}
ext_modules = []
if use_cython:  
    ext_modules += [ Extension("pyoctree.pyoctree", sources=["pyoctree/pyoctree.pyx","pyoctree/cOctree.cpp"],include_dirs=[numpy.get_include()],language="c++")]
    cmdclass.update({ 'build_ext': build_ext_compiler_check })
else:
    ext_modules += [ Extension("pyoctree.pyoctree", sources=["pyoctree/pyoctree.cpp","pyoctree/cOctree.cpp"],include_dirs=[numpy.get_include()],language="c++")]
    
setup(
    name = 'pyoctree',
    version = __version__,
    description = 'Octree structure containing 3D triangular mesh model',
    long_description = readme(),
    license = 'MIT license',
    keywords = ["octree","triangle","mesh","python","cython","ray","tracing"],    
    author = 'Michael Hogg',
    author_email = 'michael.christopher.hogg@gmail.com',
    url = "https://github.com/mhogg/pyoctree",
    download_url = "https://github.com/mhogg/pyoctree/releases", 
    packages = find_packages(),
    include_package_data = True,
    package_data = {'': ['README.rst','LICENSE.txt'],
                    'pyoctree': ['Examples\*']},				
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",                           
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        ],
    ext_modules = ext_modules,
    cmdclass = cmdclass,
)
