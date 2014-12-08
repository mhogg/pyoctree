
# cython: profile=False

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
            self.polyList.append(PyTri_Init(tri))
            
    def __dealloc__(self):
        #print "Deallocating octree"
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
                    
    def getNodeFromId(self,string nodeId):
        '''
        getNodeFromId(string nodeId)
        
        Returns a PyOctnode given the node string id i.e. '0' for root and 
        '0-0' for first branch
        '''
        cdef cOctNode *node = self.thisptr.getNodeFromId(nodeId)
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
        cdef vector[bint] ints = self.thisptr.findRayIntersects(rayList)
        cdef np.ndarray[int32,ndim=1] foundInts = np.zeros(_rayList.shape[0],dtype=np.int32)
        for i in range(_rayList.shape[0]):
            foundInts[i] = ints[i]
        return foundInts
    
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

    def __cinit__(self,parent):
        # self.thisptr will be set by global function PyOctnode_Init to point
        # to an existing cOctNode object
        self.thisptr = NULL
        self.parent  = parent
        
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
        
    def boxRayIntersect(self,np.ndarray[float,ndim=2] _rayPoints):
        '''
        boxRayIntersect(self,np.ndarray[float,ndim=2] _rayPoints)
        
        Determines if an intersection occurs between the given ray and the
        current node
        '''
        cdef int i
        cdef vector[double] p0, p1
        p0.resize(3)
        p1.resize(3)
        for i in range(3):
            p0[i] = _rayPoints[0][i]
            p1[i] = _rayPoints[1][i]
        cdef cLine ray = cLine(p0,p1,0)
        cdef bint result = self.thisptr.boxRayIntersect(ray)
        return result
        
    property isLeaf:
        def __get__(self):
            return self.thisptr.isLeafNode()  

    property branches:
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

    property polyList:
        def __get__(self):
            cdef list polyList = []
            cdef int numPolys  = self.thisptr.numPolys()
            cdef int i
            for i in range(numPolys):
                polyList.append(self.thisptr.data[i])
            return polyList
            
    property polyListAsString:
        def __get__(self):
            cdef int numPolys = self.thisptr.numPolys()
            cdef int i
            s = str(self.thisptr.data[0])
            for i in range(1, numPolys):
                s += ", " + str(self.thisptr.data[i])
            return s

    property level:
        def __get__(self):
            return self.thisptr.level  

    property nid:
        def __get__(self):
            return self.thisptr.nid  
                                 
    property numPolys:
        def __get__(self):
            return self.thisptr.numPolys()

    property size:
        def __get__(self):
            return self.thisptr.size

    property position:
        def __get__(self):
            cdef int dims = self.thisptr.position.size()
            cdef int i
            cdef np.ndarray[float64,ndim=1] position = np.zeros(3,dtype=np.float64)
            for i in range(dims):
                position[i] = self.thisptr.position[i]
            return position

    def __dealloc__(self):
        # No need to dealloc - cOctNodes are managed by cOctree
        pass


cdef class PyTri:
    cdef cTri *thisptr
    def __cinit__(self):
        self.thisptr = NULL     
    property label:
        def __get__(self):
            return self.thisptr.label
        #def __set__(self,label):
        #    self.thisptr.label = label 
    property vertices:
        def __get__(self):
            cdef np.ndarray[float64,ndim=2] vertices = np.zeros((3,3))
            cdef int i, j
            for i in range(3):
                for j in range(3):
                    vertices[i,j] = self.thisptr.vertices[i][j]
            return vertices
    property N:
        def __get__(self):
            cdef np.ndarray[float64,ndim=1] N = np.zeros(3)
            cdef int i
            for i in range(3):
                N[i] = self.thisptr.N[i]
            return N
    property D:
        def __get__(self):
            return self.thisptr.D	
    def __dealloc__(self):
        # No need to dealloc - cTris are managed by cOctree
        pass


# Need a global function to be able to point a cOctNode to a PyOctnode
cdef PyOctnode_Init(cOctNode *node, object parent):
    result = PyOctnode(parent)
    result.thisptr = node
    return result

    
# Need a global function to be able to point a cTri to a PyTri
cdef PyTri_Init(cTri *tri):
    result = PyTri()
    result.thisptr = tri
    return result
