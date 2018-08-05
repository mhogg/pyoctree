pyoctree
========

Octree structure containing a 3D triangular mesh model. To be used for
ray tracing / shadow casting.

Written in C++ for speed, but exposed to Python using Cython.

.. image:: https://img.shields.io/pypi/v/pyoctree.svg
   :target: https://pypi.python.org/pypi/pyoctree/
   :alt: Latest PyPI version
   
Details
-------

Pyoctree uses an adaptive structure, so it will automatically divide
branches to ensure that there are no more than 200 objects per leaf.

Requirements
------------

-  Python 2.7 or Python >= 3.5

Optional
--------

-  vtk >= v6.2.0 or >= v7.0 (for outputting a vtk file for viewing octree structure in Paraview)
-  A C++ compiler for building the extension module from the provided cpp file (already cythonized). Suggested compilers are:
   -  *The Microsoft C++ Compiler for Python 2.7* if using Python 2 on Windows
   -  *Microsoft Visual C++ 2015 (14.0)* if using Python 3 on Windows
   -  *gcc* on Linux
   -  *Mingw32* on Windows or Linux 
- Cython >= v0.20 and a C++ compiler to build from source     

Note that a compiler is not required if installing using the provided Python wheel.
   
Installation
------------

Intersection testing uses parallel processing via OpenMP. To use more than a 
single processor, use the provided Python wheel or compile from source using a 
compiler that supports OpenMP. Then set value of environment variable
OMP\_NUM\_THREADS to the number of desired processors.

Note that the compilers provided by the Anaconda Python distribution *do not* support OpenMP.

1. Building from source
~~~~~~~~~~~~~~~~~~~~~~~

To compile *without* OpenMP, open a command prompt, browse to the base directory containing the setup.py file and type:

.. code::

   python setup.py install
   
To compile *with* OpenMP, open a command prompt, browse to the base directory containing the setup.py file and type:

.. code::

    python setup.py install --openmp
   
2. Installation using Python wheel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the python wheel from `releases <https://github.com/mhogg/pyoctree/releases>`_ i.e. 
pyoctree-0.2.7-cp36-cp36m-win_amd64.whl for Python 3.6 on Windows 64-bit. Then, open a command 
prompt, browse to the download directory and type:

.. code::

   pip install pyoctree-0.2.7-cp36-cp36m-win_amd64.whl
   
Note that Python wheels have been built *with* OpenMP.

Usage
-----

1. Creating the octree representation of a 3D triangular mesh model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from pyoctree import pyoctree as ot
    tree = ot.PyOctree(pointCoords,connectivity)

where:

-  pointCoords is a Nx3 numpy array of floats (dtype=float) representing
   the 3D coordinates of the mesh points

-  connectivity is a Nx3 numpy array of integers (dtype=np.int32)
   representing the point connectivity of each tri element in the mesh

2. Finding intersection between mesh object and ray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The octree can be used to quickly find intersections between the object
and a ray. For example:

.. code:: python

    import numpy as np
    startPoint = [0.0,0.0,0.0]
    endPoint   = [0.0,0.0,1.0]
    rayList    = np.array([[startPoint,endPoint]],dtype=np.float32)
    intersectionFound  = tree.rayIntersection(rayList)

Examples
--------

Some Jupyter notebooks are provided in the Examples directory on how to
use pyoctree.

Help
----

If help is required, please create an issue on Github.
