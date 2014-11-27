
import vtk
import numpy as np
import octree

# Read in stl file
reader = vtk.vtkSTLReader()
reader.SetFileName('knot.stl')
reader.MergingOn()
stl = reader.GetOutput()
stl.Update()

#writer = vtk.vtkSTLWriter()
#writer.SetFileName('knot_out.stl')
#writer.SetInput(stl)
#writer.Write()

# Get array of point coordinates
numPoints   = stl.GetNumberOfPoints()
pointCoords = np.zeros((numPoints,3),dtype=float)
for i in xrange(numPoints):
    pointCoords[i,:] = stl.GetPoint(i)
    
# Get polygon connectivity
numPolys     = stl.GetNumberOfCells()
connectivity = np.zeros((numPolys,3),dtype=int)
for i in xrange(numPolys):
    tri = stl.GetCell(i)
    ids = tri.GetPointIds()
    for j in range(3):
        connectivity[i,j] = ids.GetId(j)

# 
tree = octree.PyOctree(pointCoords,connectivity)

# Visualise octree
def getOctreeRep(octnode):
    global connect, numVerts, pc
    if not octnode.isLeaf:
        for branch in octnode.branches:
            pos  = branch.position
            size = branch.size
            halfSize = size/2.0
            n1 = pos + np.array([-halfSize,-halfSize,-halfSize])
            n2 = pos + np.array([+halfSize,-halfSize,-halfSize])  
            n3 = pos + np.array([+halfSize,+halfSize,-halfSize])  
            n4 = pos + np.array([-halfSize,+halfSize,-halfSize])  
            n5 = pos + np.array([-halfSize,-halfSize,+halfSize])
            n6 = pos + np.array([+halfSize,-halfSize,+halfSize])  
            n7 = pos + np.array([+halfSize,+halfSize,+halfSize])  
            n8 = pos + np.array([-halfSize,+halfSize,+halfSize])  
            verts = [n1,n2,n3,n4,n5,n6,n7,n8]
            connect.append(range(numVerts,numVerts+8))  
            numVerts += len(verts)
            for v in verts:
                pc.InsertNextTuple(tuple(v)) 
            getOctreeRep(branch)
    
connect = []
pc = vtk.vtkFloatArray()
pc.SetNumberOfComponents(3)
numVerts = 0    
getOctreeRep(tree.root)

pnts = vtk.vtkPoints()
pnts.SetData(pc)

uGrid = vtk.vtkUnstructuredGrid()
uGrid.SetPoints(pnts)
numElems = len(connect)
for i in range(numElems):
    hexelem = vtk.vtkHexahedron()
    c = connect[i]
    for j in range(8):
        hexelem.GetPointIds().SetId(j,c[j])    
    uGrid.InsertNextCell(hexelem.GetCellType(), hexelem.GetPointIds())  

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('octree.vtu')
writer.SetInput(uGrid)
writer.SetDataModeToAscii()
writer.Write()
