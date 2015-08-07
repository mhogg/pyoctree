octree
======

Octree structure containing a triangular mesh. To be used for ray tracing / shadow casting.

Written in C++ for speed, but exposed to Python using Cython.


## Details

Octree structure is adaptive, so it will automatically divide branches to ensure that there are no more than 200 objects per leaf.

Intersection testing uses parallel processing via OpenMP. To use more than a single processor, set value of environment variable OMP_NUM_THREADS to number of desired processors.  


## Requirements

* vtk >= v6.2.0 (optional, for outputting a vtk file for viewing octree structure in Paraview)
* Cython >= v0.20 and a C++ compiler for compiling from source (optional)


## Usage

### 1. Creating the octree representation of a 3D triangular mesh 

```python
import octree
tree = octree.PyOctree(pointCoords,connectivity)
```

where:

* pointCoords is a Nx3 numpy array of floats (dtype=float) representing the 3D coordinates of the mesh points

* connectivity is a Nx3 numpy array of integers (dtype=np.int32) representing the point connectivity of each tri element in the mesh 


### 2. Finding intersection between mesh object and ray

The octree can be used to quickly find intersections between the object and a ray. For example:

```python
import numpy as np
startPoint = [0.0,0.0,0.0]
endPoint   = [0.0,0.0,1.0]
rayList    = np.array([[startPoint,endPoint]],dtype=np.float32)
intersectionFound  = tree.rayIntersection(rayList)
```


## Examples

Some IPython notebooks are provided in the Examples directory on how to use octree.


## Help

If help is required, please create an issue on Github.
