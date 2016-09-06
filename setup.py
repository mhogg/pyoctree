# -*- coding: utf-8 -*-

# Copyright (C) 2016 Michael Hogg

# This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

from setuptools import setup
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
    use_cython = False
else:
    use_cython = True

cmdclass    = {}
ext_modules = []
if use_cython:  
    ext_modules += [ Extension("pyoctree.pyoctree", sources=["pyoctree/pyoctree.pyx","pyoctree/cOctree.cpp"],include_dirs=[numpy.get_include()],extra_compile_args=['/openmp'],language="c++")]
    cmdclass.update({ 'build_ext':build_ext })
else:
    ext_modules += [ Extension("pyoctree.pyoctree", sources=["pyoctree/pyoctree.cpp","pyoctree/cOctree.cpp"],include_dirs=[numpy.get_include()],extra_compile_args=['/openmp'],language="c++")]
    
setup(
    name = 'pyoctree',
    version = __version__,
    description = 'Octree structure containing 3D triangular mesh model',
    long_description = readme(),
    license = 'MIT license',
    keywords = ["octree","triangle","mesh","python","cython"],    
    author = 'Michael Hogg',
    author_email = 'michael.christopher.hogg@gmail.com',
    url = "https://github.com/mhogg/pyoctree",
    download_url = "https://github.com/mhogg/pyoctree/releases", 
    packages = ['','pyoctree'],
    package_data = {'':['LICENSE.txt','README.md','setup.py','Examples/*']},
    classifiers = [
        "Programming Language :: Python",                                  
        "Programming Language :: Cython",         
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",                                                   
        "Development Status :: 4 - Beta",                                  
        "Environment :: Other Environment", 
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent",     
        ],
    ext_modules = ext_modules,
    cmdclass = cmdclass,
)
