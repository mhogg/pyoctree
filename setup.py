# -*- coding: utf-8 -*-

# Copyright (C) 2015 Michael Hogg

# This file is part of octree - See LICENSE.txt for information on usage and redistribution

import octree

from distutils.core import setup
from distutils.extension import Extension
import numpy
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass    = {}
ext_modules = []
if use_cython:  
    ext_modules += [ Extension("octree.octree", sources=["octree/octree.pyx","octree/cOctree.cpp"],include_dirs=[numpy.get_include()],extra_compile_args=['/openmp'],language="c++")]
    cmdclass.update({ 'build_ext':build_ext })
else:
    ext_modules += [ Extension("octree.octree", sources=["octree/octree.cpp","octree/cOctree.cpp"],include_dirs=[numpy.get_include()],extra_compile_args=['/openmp'],language="c++")]
    
setup(
    name = 'octree',
    version = octree.__version__,
    description = 'Octree structure to represent 3D model containing tri faces',
    license = 'MIT license',
    keywords = ["octree","triangle","mesh","python","cython"],    
    author = 'Michael Hogg',
    author_email = 'michael.christopher.hogg@gmail.com',
    url = "https://github.com/mhogg/octree",
    download_url = "https://github.com/mhogg/octree/releases", 
    packages = ['','octree'],
    package_data = {'':['LICENSE.txt','README.md','setup.py','Examples/*']},
    classifiers = [
        "Programming Language :: Python",                                  
        "Programming Language :: Cython",         
        "Programming Language :: Python :: 2",             
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",                                                   
        "Development Status :: 4 - Beta",                                  
        "Environment :: Other Environment", 
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent",     
        ],
    ext_modules = ext_modules,
    cmdclass = cmdclass,
    long_description = """Octree structure to represent 3D model containing tri faces""",
)
