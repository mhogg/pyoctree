pyoctree
========

Octree structure containing a 3D triangular mesh model. To be used for
ray tracing / shadow casting.

Written in C++ for speed, but exposed to Python using Cython.

.. image:: https://img.shields.io/pypi/v/pyoctree.svg
   :target: https://pypi.python.org/pypi/pyoctree/
   :alt: Latest PyPI version
   
.. image:: https://img.shields.io/pypi/dm/pyoctree.svg
   :target: https://pypi.python.org/pypi/pyoctree/
   :alt: Number of PyPI downloads
   
Details
-------

Pyoctree uses an adaptive structure, so it will automatically divide
branches to ensure that there are no more than 200 objects per leaf.

Intersection testing uses parallel processing via OpenMP. To use more
than a single processor, set value of environment variable
OMP\_NUM\_THREADS to number of desired processors.

Requirements
------------

-  vtk >= v6.2.0 (optional, for outputting a vtk file for viewing octree
   structure in Paraview)
-  Cython >= v0.20 and a C++ compiler for building the extension module. The Microsoft C++ 
   Compiler for Python 2.7 or Mingw32 can both be used. Note that this is not required if
   installing using the Python wheel.
   
Installation
------------

1. Building from source
~~~~~~~~~~~~~~~~~~~~~~~

In a command prompt, browse to the base directory containing the setup.py file and type:

.. code::

   python setup.py install

2. Installation using Python wheel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the python wheel from `releases <https://github.com/mhogg/pyoctree/releases>`_ i.e. pyoctree-0.1.2-cp27-cp27m-win_amd64.whl. Then, open a command prompt, browse to the download directory and type:

.. code::

   pip install pyoctree-0.1.2-cp27-cp27m-win_amd64.whl

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
