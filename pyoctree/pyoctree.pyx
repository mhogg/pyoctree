# -*- coding: utf-8 -*-

# Copyright (C) 2017 Michael Hogg

# This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

# cython: profile=False

from __future__ import print_function
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.string cimport string
from cython.operator cimport dereference as deref, preincrement as inc, predecrement as dec
ctypedef np.float64_t float64
ctypedef np.int32_t int32
 

cdef extern from "cOctree.h":

    cdef cppclass Intersection:
        double s
        vector[double] p
    
    cdef cppclass cLine:
        vector[double] p0,p1,dir
        cLine()
        cLine(vector[double] p0,vector[double] p1_dir,int isP1orDir)
        void getDir()
        void getP1()
            
    cdef cppclass cTri:
        int label
        vector[vector[double]] vertices
        vector[double] N
        double D
        cTri()
        cTri(int label, vector[vector[double]] vertices)
        void getN()
        void getD()
        
    cdef cppclass cOctNode:
        double size
        int level
        string nid
        vector[double] position
        vector[cOctNode] branches
        vector[int] data
        cOctNode()   
        int numPolys()
        bint isLeafNode()
        bint boxRayIntersect(cLine &ray)
        
    cdef cppclass cOctree:
        cOctree(vector[vector[double]] vertexCoords3D, vector[vector[int]] polyConnectivity)
        int numPolys()
        cOctNode root
        vector[Intersection] findRayIntersect(cLine ray)
        cOctNode* getNodeFromId(string nodeId)
        vector[cTri] polyList		
        vector[cOctNode*] getNodesFromLabel(int polyLabel)
        vector[bint] findRayIntersects(vector[cLine] &rayList)
        vector[bint] findRayIntersectsSorted(vector[cLine] &rayList)
        set[int] getListPolysToCheck(cLine &ray)
        vector[cOctNode*] getSortedNodesToCheck(cLine &ray)
    
    
cdef class PyOctree:

    cdef cOctree *thisptr
    cdef public PyOctnode root
    cdef public list polyList

    def __cinit__(self,double[:,::1] _vertexCoords3D, int[:,::1] _polyConnectivity):
    
        cdef int i, j
        cdef vector[double] coords
        cdef vector[vector[double]] vertexCoords3D
        
        coords.resize(3)
        for i in range(_vertexCoords3D.shape[0]):
            for j in range(3):
                coords[j] = _vertexCoords3D[i,j]
            vertexCoords3D.push_back(coords)
                
        cdef vector[int] connect
        cdef vector[vector[int]] polyConnectivity
        
        connect.resize(3)
        for i in range(_polyConnectivity.shape[0]):
            for j in range(3):
                connect[j] = _polyConnectivity[i,j]
            polyConnectivity.push_back(connect)
            
        # Create cOctree
        self.thisptr = new cOctree(vertexCoords3D,polyConnectivity)
            
        # Get root node
        cdef cOctNode *node = &self.thisptr.root
        self.root = PyOctnode_Init(node,self)
        node = NULL
        
        # Get polyList
        self.polyList = []
        cdef cTri *tri 
        for i in range(self.thisptr.polyList.size()):
            tri = &self.thisptr.polyList[i]
            self.polyList.append(PyTri_Init(tri,self))
            
    def __dealloc__(self):
        #print("Deallocating octree")
        del self.thisptr
        
    def getNodesFromLabel(self,int label):
        '''
        getNodesFromLabel(int label)
        
        Returns a list of PyOctnodes that the polygon with the given label
        is contained within 
        '''
        cdef cOctNode *node = NULL
        cdef vector[cOctNode*] nodes = self.thisptr.getNodesFromLabel(label)
        cdef int i
        cdef list nodeList = []
        for i in range(nodes.size()):
            node = nodes[i]
            nodeList.append(PyOctnode_Init(node,self))
        if len(nodeList)==1:
            return nodeList[0]
        else:
            return nodeList     
            
    def getNodesFromRay(self,np.ndarray[float,ndim=2] _rayPoints):
        '''
        getNodesFromRay(np.ndarray[float,ndim=2] _rayPoints)
        
        Returns a list of PyOctnodes that intersect with the given ray. The 
        nodes are sorted in order of distance from the ray origin
        '''
        cdef int i
        cdef vector[double] p0, p1
        p0.resize(3)
        p1.resize(3)
        for i in range(3):
            p0[i] = _rayPoints[0][i]
            p1[i] = _rayPoints[1][i]
        cdef cLine ray = cLine(p0,p1,0)
        cdef vector[cOctNode*] nodes = self.thisptr.getSortedNodesToCheck(ray)
        cdef cOctNode *node = NULL
        cdef list nodeList = []
        for i in range(nodes.size()):
            node = nodes[i]
            nodeList.append(PyOctnode_Init(node,self))
        if len(nodeList)==1:
            return nodeList[0]
        else:
            return nodeList        
                    
    def getNodeFromId(self,str nodeId):
        '''
        getNodeFromId(str nodeId)
        
        Returns a PyOctnode given the node string id i.e. '0' for root and 
        '0-0' for first branch
        '''
        # Encode str nodeId to byte string before passing to C++ code
        cdef cOctNode *node = self.thisptr.getNodeFromId(nodeId.encode())
        if node is NULL:
            return None
        else:
            return PyOctnode_Init(node,self)   

    def getListOfPolysToCheck(self,np.ndarray[float,ndim=2] _rayPoints):
        '''
        getListOfPolysToCheck(np.ndarray[float,ndim=2] _rayPoints)
        
        Returns a list of polygons that should be tested for intersections with
        the given ray
        '''
        cdef int i
        cdef vector[double] p0, p1
        p0.resize(3)
        p1.resize(3)
        for i in range(3):
            p0[i] = _rayPoints[0][i]
            p1[i] = _rayPoints[1][i]
        cdef cLine ray = cLine(p0,p1,0)
        cdef set[int] polySetToCheck = self.thisptr.getListPolysToCheck(ray)
        cdef int numPolys = polySetToCheck.size()
        cdef set[int].iterator it
        it = polySetToCheck.begin()
        s  = str(deref(it)); inc(it)
        while it!=polySetToCheck.end():
            s += ", " + str(deref(it))
            inc(it)
        return s
        
    def rayIntersection(self,np.ndarray[float,ndim=2] _rayPoints):
        '''
        rayIntersection(np.ndarray[float,ndim=2] _rayPoints)
        
        Finds and returns a list of all intersection points between the tree
        polys and the given ray 
        '''	
        cdef int i
        cdef vector[double] p0, p1
        p0.resize(3)
        p1.resize(3)
        for i in range(3):
            p0[i] = _rayPoints[0][i]
            p1[i] = _rayPoints[1][i]
        cdef cLine ray = cLine(p0,p1,0)
        cdef vector[Intersection] intersectList = self.thisptr.findRayIntersect(ray)
        numInts = intersectList.size()
        intList = []
        for i in range(numInts):
            intsect = Intersect()
            intsect.SetValues(intersectList[i].p,intersectList[i].s)
            intList.append(intsect)
        return intList

    def rayIntersections(self,np.ndarray[float,ndim=3] _rayList):
        '''
        rayIntersections(np.ndarray[float,ndim=3] _rayList
        
        For every ray in the list provided, returns a corresponding array of 
        integers indicating if a intersection occurred or not. This array can
        be used to create an shadow image of the tree mesh part
        '''
        cdef int i,j
        cdef vector[double] p0, p1
        cdef vector[cLine] rayList
        p0.resize(3)
        p1.resize(3)
        for i in range(_rayList.shape[0]):
            for j in range(3):
                p0[j] = _rayList[i][0][j]
                p1[j] = _rayList[i][1][j]
            rayList.push_back(cLine(p0,p1,0))
        cdef vector[bint] ints = self.thisptr.findRayIntersectsSorted(rayList)
        cdef np.ndarray[int32,ndim=1] foundInts = np.zeros(_rayList.shape[0],dtype=np.int32)
        for i in range(_rayList.shape[0]):
            foundInts[i] = ints[i]
        return foundInts
        
    cpdef int getNumberOfNodes(self):
        '''
        getNumberOfNodes()
        
        Returns the number of Octnodes in the Octree
        '''
        cdef int numNodes = 0
        self.countNodes(self.root,numNodes)
        return numNodes
        
    cdef void countNodes(self,PyOctnode node, int &numNodes):
        '''
        Utility function for getNumberOfNodes. Recursively counts number of Octnodes in Octree
        '''
        cdef PyOctnode branch
        (&numNodes)[0] = (&numNodes)[0] + 1  # Syntax here is fix for Cython bug
        if not node.isLeaf:
            for branch in node.branches:
                self.countNodes(branch,numNodes)
                
    def getOctreeRep(self,fileName='octree.vtu'):
        '''
        getOctreeRep(fileName='octree.vtu')
        
        Output a vtk representation of the Octree that can be viewed in Paraview. This
        requires the vtk python module
        '''
        
        try: import vtk
        except:
            print('Error: Cannot import required vtk module')
            return
        
        def getTree(node):
            if node.level==1:
                getNodeRep(node)             
            for branch in node.branches:
                getNodeRep(branch)
                getTree(branch)        
                
        def getNodeRep(node):
            offsets = {0:(-1,-1,-1),1:(+1,-1,-1),2:(+1,+1,-1),3:(-1,+1,-1),
                       4:(-1,-1,+1),5:(+1,-1,+1),6:(+1,+1,+1),7:(-1,+1,+1)} 
            connect = []
            for i in range(8):
                vi = node.position + 0.5*node.size*np.array(offsets[i])
                connect.append(len(vertexCoords))
                vertexCoords.append(vi)                       
            vertexConnect.append(connect)               
               
        # For every node in tree, get vertex coordinates and connectivity
        vertexCoords  = []     
        vertexConnect = []
        getTree(self.root)
        
        # Convert to vtk unstructured grid
        uGrid = vtk.vtkUnstructuredGrid()
        
        # 1. Convert vertices and add to unstructured grid
        coords = vtk.vtkFloatArray()
        coords.SetNumberOfComponents(3)
        for v in vertexCoords:
            coords.InsertNextTuple(tuple(v))
        vertices = vtk.vtkPoints()    
        vertices.SetData(coords)
        uGrid.SetPoints(vertices)
        
        # 2. Add element data to unstructured grid
        numElems = len(vertexConnect)
        for i in xrange(numElems):
            hexelem = vtk.vtkHexahedron()
            c = vertexConnect[i]
            for j in range(8):
                hexelem.GetPointIds().SetId(j,c[j])    
            uGrid.InsertNextCell(hexelem.GetCellType(), hexelem.GetPointIds()) 
            
        # Write unstructured grid to file
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(fileName)
        writer.SetInputData(uGrid)
        writer.SetDataModeToAscii()
        writer.Write()      
        
    property numPolys:
        def __get__(self):
            return self.thisptr.numPolys()
        
        
cdef class Intersect:
    cdef public double s
    cdef public np.ndarray p
    def __init__(self):
        self.s = 0.0
        self.p = np.zeros(3,dtype=float)
    cdef SetValues(self,vector[double] p, double s):
        cdef int i
        for i in range(3):
            self.p[i] = p[i]
        self.s = s
        
              
cdef class PyOctnode:
    
    cdef cOctNode *thisptr
    cdef public object parent

    def __cinit__(self,parent=None):
        self.thisptr = NULL
        self.parent  = parent
        # If parent is None, then create a new cOctNode instance. Otherwise, assume
        # that PyOctree instance is managed by the PyOctree
        if self.parent is None:
            self.thisptr = new cOctNode()
        
    def __dealloc__(self):
        # If parent is None, then cOctNodes are not managed by cOctree
        if self.parent is None:
            del self.thisptr
        
    def hasPolyLabel(self,label):
        '''
        hasPolyLabel(label)
        
        Checks if poly with given label is in the current node
        '''
        cdef int numPolys  = self.thisptr.numPolys()
        cdef int i
        for i in range(numPolys):
            if self.thisptr.data[i]==label:
                return True
        return False
        
    cdef printWarningMsg(self,s):
        print('PyOctnode is managed by PyOctree: %s is read-only' % s)
        
    def __str__(self):
        return "<%s, Id: %s, isLeaf: %r, numPolys: %d>" % ('PyOctnode', self.nid, self.isLeaf, self.numPolys)
        
    def __repr__(self):
        return "<%s %s>" % ('PyOctnode', self.nid)
        
    property isLeaf:
        '''Checks if node is a leaf (has no branches)'''
        def __get__(self):
            return self.thisptr.isLeafNode()
        def __set__(self,_isLeaf):
            pass

    property branches:
        '''
        Returns a list of all octNode branches. If node is a leaf, returns
        an empty list
        '''
        def __get__(self):
            branches = []
            cdef int numBranches = self.thisptr.branches.size()
            cdef int i
            cdef cOctNode *node = NULL
            for i in range(numBranches):
                node = &self.thisptr.branches[i]
                branches.append(PyOctnode_Init(node,self))
            node = NULL
            return branches
        def __set__(self,_branches):
            if self.parent is not None:
                self.printWarningMsg('PyOctnode.branches')
            else: pass

    property polyList:
        '''
        Returns a list of all the polygon indices (with respect to the cOctree
        polyList) within the given node
        '''
        def __get__(self):
            cdef list polyList = []
            cdef int numPolys  = self.thisptr.numPolys()
            cdef int i
            for i in range(numPolys):
                polyList.append(self.thisptr.data[i])
            return polyList
        def __set__(self,_polyList):
            pass
            
    property polyListAsString:
        '''
        Similar to polyList property, but returns a comma delimited string
        rather than a list
        '''
        def __get__(self):
            cdef int numPolys = self.thisptr.numPolys()
            cdef int i
            if self.numPolys == 0: return ""
            s = str(self.thisptr.data[0])
            for i in range(1, numPolys):
                s += ", " + str(self.thisptr.data[i])
            return s            
        def __set__(self,_polyListAsString):
            pass

    property level:
        '''octNode level'''
        def __get__(self):
            return self.thisptr.level
        def __set__(self,_level):
            if self.parent is None: self.thisptr.level = _level
            else: self.printWarningMsg('PyOctnode.level')

    property nid:
        '''octNode node id'''
        def __get__(self):
            # Encode to convert from byte code to string
            return self.thisptr.nid.encode()
        def __set__(self,_nid):
            if self.parent is None: self.thisptr.nid = _nid
            else: self.printWarningMsg('PyOctnode.nid')
            
    property numPolys:
        '''Number of polygons in given octNode'''
        def __get__(self):
            return self.thisptr.numPolys()
        def __set__(self,_numPolys):
            pass

    property size:
        '''Size of octNode bounding box'''
        def __get__(self):
            return self.thisptr.size
        def __set__(self,_size):
            if self.parent is None: self.thisptr.size = _size
            else: self.printWarningMsg('PyOctnode.size')

    property position:
        '''Coordinates of octNode centre'''
        def __get__(self):
            cdef int dims = self.thisptr.position.size()
            cdef int i
            cdef np.ndarray[float64,ndim=1] position = np.zeros(3,dtype=np.float64)
            for i in range(dims):
                position[i] = self.thisptr.position[i]
            return position
        def __set__(self,_position):
            cdef int i
            if self.parent is None:
                _position = np.array(_position)
                if not (_position.shape == (3,) or _position.shape == (1,3)):
                    print('Error: position must be a 1x3 array')
                    return
                for i in range(3):
                    self.thisptr.position[i] = _position[i]
            else: self.printWarningMsg('PyOctnode.position')


cdef class PyTri:
    cdef cTri *thisptr
    cdef public object parent
    def __cinit__(self,parent=None):
        self.thisptr = NULL
        self.parent  = parent
        # If parent is None, then create a new cTri instance. Otherwise, assume
        # that PyTri instance is managed by the PyOctree
        if self.parent is None:
            self.thisptr = new cTri()
    def __dealloc__(self):
        # If parent is None, then cTris are not managed by cOctree
        if self.parent is None:
            del self.thisptr
    def __str__(self):
        return "<%s %d>" % ('PyTri', self.label)
    def __repr__(self):
        return "<%s %d>" % ('PyTri', self.label)
    cdef printWarningMsg(self,s):
        print('PyTri is managed by PyOctree: %s is read-only' % s)
    property label:
        '''Tri label'''
        def __get__(self):
            return self.thisptr.label
        def __set__(self,_label):
            self.printWarningMsg('PyTri.label')
    property vertices:
        '''Array of tri vertices'''
        def __get__(self):
            cdef np.ndarray[float64,ndim=2] vertices = np.zeros((3,3))
            cdef int i, j
            for i in range(3):
                for j in range(3):
                    vertices[i,j] = self.thisptr.vertices[i][j]
            return vertices
        def __set__(self,_vertices):
            if self.parent is not None:
                self.printWarningMsg('PyTri.vertices')
                return                
            _vertices = np.array(_vertices,dtype=np.float64)
            if _vertices.shape != (3,3):
                print('Error: vertices must be a 3x3 array')
                return
            cdef int i,j
            for i in range(3):
                for j in range(3):
                    self.thisptr.vertices[i][j] = _vertices[i,j]  
    property N:
        '''Tri face normal'''
        def __get__(self):
            cdef np.ndarray[float64,ndim=1] N = np.zeros(3)
            cdef int i
            self.thisptr.getN()            
            for i in range(3):
                N[i] = self.thisptr.N[i]
            return N
        def __set__(self,_N):
            pass
    property D:
        '''Perp. distance from tri face to origin'''
        def __get__(self):
            self.thisptr.getD()
            return self.thisptr.D
        def __set__(self,_D):
            pass


# Need a global function to be able to point a cOctNode to a PyOctnode
cdef PyOctnode_Init(cOctNode *node, object parent):
    result = PyOctnode(parent)
    result.thisptr = node
    return result


# Need a global function to be able to point a cTri to a PyTri
cdef PyTri_Init(cTri *tri, object parent):
    result = PyTri(parent)
    result.thisptr = tri
    return result
