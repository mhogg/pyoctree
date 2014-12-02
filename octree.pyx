
# cython: profile=True

#import string as pystring
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref, preincrement as inc, predecrement as dec
ctypedef np.float64_t float64
ctypedef np.int32_t int32


cdef extern from "math.h":
    float fabs(float val)
    float acos(float val)
    float asin(float val)
    float atan(float val)
    float fmin(float a, float b)
    float fmax(float a, float b)
 

cdef extern from "cOctree.h":
    
    cdef cppclass CLine:
        vector[double] p0,p1,dir
        CLine()
        CLine(vector[double] p0,vector[double] p1_dir,int isP1orDir)
        void getDir()
        void getP1()
            
    cdef cppclass cTri:
        int label
        vector[vector[double]] vertices
        vector[double] N
        double D
        #cTri()
        cTri(int label, vector[vector[double]] vertices)
        void getN()
        void getD()
        bint isInNode(vector[vector[double]] &lowUpp)
        
    cdef cppclass cOctNode:
        double size
        int level
        string nid
        vector[double] position
        vector[cOctNode] branches
        vector[int] data        
        int numPolys()
        bint isLeafNode()
        
    cdef cppclass cOctree:
        cOctree(vector[vector[double]] vertexCoords3D, vector[vector[int]] polyConnectivity)
        int numPolys()
        cOctNode root
        vector[cTri] polyList
        cOctNode* getNodeFromLabel(int polyLabel)
        cOctNode* getNodeFromId(string nodeId)
    
    #cdef double dotProduct( vector[double] v1, vector[double] v2 )
    #cdef vector[double] crossProduct(vector[double] v1, vector[double] v2)
    #cdef vector[double] vectSubtract( vector[double] a, vector[double] b )
    #cdef vector[double] vectAdd( vector[double] a, vector[double] b )
    #cdef vector[double] vectAdd( vector[double] a, vector[double] b, double sf )
    #cdef double distBetweenPoints( vector[double] a, vector[double] b )   
        
                    
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
        
    def getNodeFromLabel(self,int label):
        cdef cOctNode *node = self.thisptr.getNodeFromLabel(label)
        if node is NULL:
            return None
        else:
            return PyOctnode_Init(node,self)
            
    def getNodeFromId(self,string nodeId):
        cdef cOctNode *node = self.thisptr.getNodeFromId(nodeId)
        if node is NULL:
            return None
        else:
            return PyOctnode_Init(node,self)            
    
    property numPolys:
        def __get__(self):
            return self.thisptr.numPolys()
                
    def __dealloc__(self):
        print "Deallocating octree"
        del self.thisptr            

              
cdef class PyOctnode:
    
    cdef cOctNode *thisptr
    cdef public object parent

    def __cinit__(self,parent):
        # self.thisptr will be set by global function PyOctnode_Init to point
        # to an existing cOctNode object
        self.thisptr = NULL
        self.parent  = parent

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
        def __set__(self,label):
            self.thisptr.label = label 
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
