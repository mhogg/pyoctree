
#include "cOctree.h"

// ------------------------------------------------------

cLine::cLine() 
{
    // Default line is unit vector along x-axis
    p0.resize(3,0.0); p1.resize(3,0.0); dir.resize(3,0.0);
    p1[0]=1.0; dir[0]=1.0;
}

cLine::cLine(vector<double> &_p0, vector<double> &p1_dir, int isP1orDir)
{
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = _p0;
    if (isP1orDir==0) {
        p1 = p1_dir;
        getDir(); }
    else if (isP1orDir==1) {
        dir = p1_dir;
        getP1(); }
}

cLine::~cLine() {}

void cLine::getDir()
{
    vector<double> p0p1(3); double dmag=0.0;
    for (unsigned int i=0; i<3; i++) {
        p0p1[i] = p1[i]-p0[i];
        dmag   += pow(p0p1[i],2.0); }
    dmag = sqrt(dmag);
    dir  = p0p1;
    for (vector<double>::iterator it=dir.begin(); it!=dir.end(); ++it)
        *it /= dmag; 
}

void cLine::getP1()
{
    vector<double> p1(3);
    for (unsigned int i=0; i<3; i++)
        p1[i] = p0[i]+dir[i];
}

// ------------------------------------------------------

cTri::cTri()
{
    // Default tri
    label = 0;
    vertices.resize(3);
    for (vector<vector<double> >::iterator it=vertices.begin(); it!=vertices.end(); ++it)
        (*it).resize(3,0.0);
    vertices[1][0]=1.0; vertices[2][1]=1.0;
    getN();
    getD();
    getLowerVert();
    getUpperVert();    
}

cTri::cTri(int _label, vector<vector<double> > _vertices)
{
    label    = _label;
    vertices = _vertices;
    getN();
    getD();
    getLowerVert();
    getUpperVert();
}

cTri::~cTri() {
    //cout << "Destroying cTri" << endl;
}

void cTri::getN()
{
    vector<vector<double> > v = vertices;
    vector<double> v1(3),v2(3),v3;
    for (unsigned int i=0; i<3; i++) {
        v1[i] = v[0][i] - v[1][i];
        v2[i] = v[0][i] - v[2][i]; }
    v3 = crossProduct(v1,v2);
    double v3mag = sqrt(dotProduct(v3,v3));
    N = v3;
    for (vector<double>::iterator it=N.begin(); it!=N.end(); ++it)
        *it /= v3mag;
}

void cTri::getD()
{
    D = dotProduct(vertices[0],N);
}

void cTri::getLowerVert()
{
    lowVert.resize(3,1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            if (vertices[i][j] < lowVert[j]) 
            {
                lowVert[j] = vertices[i][j]; 
            }
        }
    } 
}

void cTri::getUpperVert()
{
    uppVert.resize(3,-1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            if (vertices[i][j] > uppVert[j]) 
            {
                uppVert[j] = vertices[i][j]; 
            }
        }
    } 
} 

bool cTri::isInNode(cOctNode &node)
{
    if (lowVert[0] > node.upp[0]) return false;
    if (lowVert[1] > node.upp[1]) return false;
    if (lowVert[2] > node.upp[2]) return false;
    if (uppVert[0] < node.low[0]) return false;
    if (uppVert[1] < node.low[1]) return false;
    if (uppVert[2] < node.low[2]) return false;
    return true;
}

bool cTri::isPointInTri(vector<double> &p)
{
    // Determines if point p is within a triangle by computing and
    // testing the barycentric coordinates (u, v, w) of p

    // Find Barycentric coordinates of point (u,v,w)
    vector<double> v0 = vectSubtract(vertices[1],vertices[0]);
    vector<double> v1 = vectSubtract(vertices[2],vertices[0]);
    vector<double> v2 = vectSubtract(p,vertices[0]);
    double d00 = dotProduct(v0, v0);
    double d01 = dotProduct(v0, v1);
    double d11 = dotProduct(v1, v1);
    double d20 = dotProduct(v2, v0);
    double d21 = dotProduct(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    // Use Barycentric coordinates to work out if point lies within the tri element
    double tol = 1.0e-06;	
    return ((v>=tol) && (w>=tol) && (u>=tol));
}

bool cTri::rayPlaneIntersectPoint(cLine &ray)
{
    double tol  = 1.0e-06;
    double sDen = dotProduct(ray.dir,N);
    if (fabs(sDen)> tol) // Normals cannot be perpendicular such that dot product equals 0
    {
        double sNum = D - dotProduct(ray.p0,N);
        double s = sNum / sDen;
        vector<double> p = vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    } 
    return false;
}

bool cTri::rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s)
{
    double tol  = 1.0e-06;
    double sDen = dotProduct(ray.dir,N);
    if (fabs(sDen)> tol) // Normals cannot be perpendicular such that dot product equals 0
    {
        double sNum = D - dotProduct(ray.p0,N);
        s = sNum / sDen;
        p = vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    } 
    return false;
}

// ------------------------------------------------------

cOctNode::cOctNode()
{  
    setupConstants();
    level = 0;
    nid   = "";
    size  = 1.0;
    position.resize(3,0.0);
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);
}

cOctNode::cOctNode(int _level, string _nid, vector<double> _position, double _size)
{
    setupConstants();
    level    = _level;
    nid      = _nid;
    position = _position;
    size     = _size;
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);    
}

cOctNode::~cOctNode() {
    //cout << "Calling destructor for cOctnode " << nid << endl;
}

void cOctNode::setupConstants() 
{
    NUM_BRANCHES_OCTNODE = 8; 
    MAX_OCTNODE_OBJECTS  = 100;   
}

bool cOctNode::isLeafNode() { return branches.size()==0; }

void cOctNode::getLowUppVerts() 
{
    low.resize(3);
    upp.resize(3);
    double halfSize = size/2.0;
    for (int i=0; i<3; i++) {
        low[i] = position[i] - halfSize;
        upp[i] = position[i] + halfSize;
    }
}

void cOctNode::addPoly(int _indx) { data.push_back(_indx); }

int cOctNode::numPolys() { return data.size(); }

void cOctNode::addNode(int _level, string _nid, vector<double> _position, double _size)
{
    branches.push_back(cOctNode(_level,_nid,_position,_size));
}

bool cOctNode::sphereRayIntersect(cLine &ray)
{
    // Quick test for determining if a ray is *likely* to intersect a given node

    // Radius of sphere that contains node
    double radius = distBetweenPoints(low,position);
    
    // Project centre of sphere (node.position) onto ray
    vector<double> oc = vectSubtract(position, ray.p0);
    double s = dotProduct(oc,ray.dir);
    vector<double> projpnt = vectAdd(ray.p0, ray.dir, s);
    double dist = distBetweenPoints(projpnt,position);
    
    // If distance between spherical centre and projected point is 
    // less than the radius of the sphere, then an intersection is
    // *possible*
    return (dist<=radius);
}

// ------------------------------------------------------


cOctree::cOctree(vector<vector<double> > _vertexCoords3D, vector<vector<int> > _polyConnectivity)
{
    MAX_OCTREE_LEVELS = 10;
    vertexCoords3D    = _vertexCoords3D;
    polyConnectivity  = _polyConnectivity;
    int _offsets[][3] = {{-1,-1,-1},{+1,-1,-1},{-1,+1,-1},{+1,+1,-1},
                         {-1,-1,+1},{+1,-1,+1},{-1,+1,+1},{+1,+1,+1}};
    
    for (int i=0; i<8; i++) {
        for (int j=0; j<3; j++) {
            branchOffsets[i][j] = _offsets[i][j];
        }
    }

    setupPolyList();    
    vector<double> position = getPositionRoot();
    double size = getSizeRoot();
    root = cOctNode(1,"0", position, size);
    insertPolys();
}

void cOctree::setupPolyList()
{
    int indx;
    vector<vector<double> > vertices(3,vector<double>(3,0.0));
    
    polyList.reserve(polyConnectivity.size());
    for (unsigned int i=0; i<polyConnectivity.size(); i++) {
        for (int j=0; j<3; j++) {
            indx = polyConnectivity[i][j];
            vertices[j] = vertexCoords3D[indx]; }
        polyList.push_back(cTri(i,vertices)); 
    }
}

int cOctree::numPolys() { return polyList.size(); }

void cOctree::insertPoly(cOctNode &node, cTri &poly)
{
    if (node.isLeafNode()) {
    
        if (poly.isInNode(node)) {
        
            if (node.numPolys() < node.MAX_OCTNODE_OBJECTS) {
                node.addPoly(poly.label);
            } else {
                node.addPoly(poly.label);
                if (node.level < MAX_OCTREE_LEVELS) {
                    splitNodeAndReallocate(node);
                }
            }
        }
        
    } else {
      
        for (unsigned int i=0; i<node.branches.size(); i++) {
            insertPoly(node.branches[i],poly);
        }
        
    }
}

void cOctree::insertPolys()
{
    for (int i=0; i<numPolys(); i++) {
        insertPoly(root,polyList[i]);
    }
}

vector<double> cOctree::getPositionRoot() {

    // Get low and upp
    vector<double> low, upp, position(3);
    low = vertexCoords3D[0];
    upp = vertexCoords3D[0];
    for (unsigned int i=1; i<vertexCoords3D.size(); i++) {
        for (int j=0; j<3; j++) {
            if (vertexCoords3D[i][j] < low[j]) { low[j] = vertexCoords3D[i][j]; }
            if (vertexCoords3D[i][j] > upp[j]) { upp[j] = vertexCoords3D[i][j]; }
        }
    }
    // Center of node is average of low and upp
    for (int i=0; i<3; i++) {
        position[i] = 0.5 * (low[i]+upp[i]);
    }
    return position;
}

double cOctree::getSizeRoot() {

    // Get low and upp
    vector<double> low, upp, range;
    low = vertexCoords3D[0];
    upp = vertexCoords3D[0];
    for (unsigned int i=1; i<vertexCoords3D.size(); i++) {
        for (int j=0; j<3; j++) {
            if (vertexCoords3D[i][j] < low[j]) { low[j] = vertexCoords3D[i][j]; }
            if (vertexCoords3D[i][j] > upp[j]) { upp[j] = vertexCoords3D[i][j]; }
        }
    }
    // Range is the size of the node in each coord direction
    range = vectSubtract(upp,low);
    double size = range[0];
    for (int i=1; i<3; i++) {
        if (range[i] > size) { size = range[i]; }
    }
    // Scale up size of node by 5%
    size *= 1.05;
    return size;
}

void cOctree::splitNodeAndReallocate(cOctNode &node)
{
    // Split node into 8 branches
    vector<double> position(3);
    for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
        for (int j=0; j<3; j++) {
            position[j] = node.position[j] + 0.25*node.size*branchOffsets[i][j]; }
        string nid = node.nid + "-" + NumberToString(i);
        node.addNode(node.level+1,nid,position,0.5*node.size); 
    }
    
    // Reallocate date from node to branches
    for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
        for (int j=0; j<node.numPolys(); j++) {
            int indx = node.data[j];
            if (polyList[indx].isInNode(node.branches[i])) {
                if (node.branches[i].numPolys() < node.MAX_OCTNODE_OBJECTS) {
                    node.branches[i].addPoly(indx);
                } else {
                    splitNodeAndReallocate(node.branches[i]);
                }
            }
        }
    }
    node.data.resize(0);
}

vector<cOctNode*> cOctree::getNodesFromLabel(int polyLabel)
{
    // Function for finding all the nodes that contains tri with given label 
    vector<cOctNode*> nodeList;
    findBranchesByLabel(polyLabel,root,nodeList);
    return nodeList;
}

void cOctree::findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList)
{
    // Recursive function used by getNodesFromLabel
    if (node.isLeafNode()) {
        vector<int>::iterator it;
        it = find(node.data.begin(),node.data.end(),polyLabel);
        if (it != node.data.end()) { nodeList.push_back(&node); }
    } else {
        for (unsigned int i=0; i<node.branches.size(); i++) {
            findBranchesByLabel(polyLabel, node.branches[i], nodeList);
        }
    }
}

cOctNode* cOctree::getNodeFromId(string nodeId)
{
    return findBranchById(nodeId,root);
}

cOctNode* cOctree::findBranchById(string nodeId, cOctNode &node)
{
    if (nodeId.compare(node.nid)==0) {
        return &node;
    } else {
        for (unsigned int i=0; i<node.branches.size(); i++) {
            cOctNode *branch = findBranchById(nodeId, node.branches[i]);
            if (branch != NULL) { return branch; }
        }
    }
    return NULL;
}

cOctree::~cOctree() 
{
    //cout << "Destroying the cOctree" << endl;
}

set<int> cOctree::getListPolysToCheck(cLine &ray)
{
    // Returns a list of all polygons that are within OctNodes hit by a given ray
    set<int> intTestPolys;
    getPolysToCheck(root,ray,intTestPolys);
    return intTestPolys;
}

void cOctree::getPolysToCheck(cOctNode &node, cLine &ray, set<int> &intTestPolys)
{
    // Utility function for getListPolysToCheck. Finds all OctNodes hit by a given ray
    // and returns a list of the objects contained within
    if (node.sphereRayIntersect(ray)) {
        if (node.isLeafNode()) {
            for (int i=0; i<node.numPolys(); i++) {
                intTestPolys.insert(node.data[i]); }
        } else {
            for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
                getPolysToCheck(node.branches[i],ray,intTestPolys);
            } 
        }
    }
}

vector<Intersection> cOctree::findRayIntersect(cLine &ray)
{   
    // Get polys to check
    set<int> polyListCheck = getListPolysToCheck(ray);
    
    // Loop through all polys in check list to find a possible intersection
    vector<Intersection> intersectList;
    set<int>::iterator it;
    vector<double> ip;
    double s;
    for (it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
        int polyLabel = *it;
        if (polyList[polyLabel].rayPlaneIntersectPoint(ray,ip,s)) {
            intersectList.push_back(Intersection(ip,s)); } 
    }
    
    // Sort list in terms of distance of the intersection from the ray origin
    sort(intersectList.begin(),intersectList.end());
    
    return intersectList;
}

vector<int> cOctree::findRayIntersects(vector<cLine> &rayList)
{
    int numRays = rayList.size();
    vector<int> foundIntsects(numRays,0);
    #pragma omp parallel for
    for (int i=0; i<numRays; i++) 
    {
        cLine *ray = &rayList[i]; 
        set<int> polyListCheck = getListPolysToCheck(*ray);
        for (set<int>::iterator it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
            int polyLabel = *it;
            if (polyList[polyLabel].rayPlaneIntersectPoint(*ray)) {
                foundIntsects[i] = 1; break; } 
        }
    }
    return foundIntsects;
}

// ------------------------------------------------------

double dotProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates dot product v1.v2
    double dp=0.0;
    for (unsigned int i=0; i<3; i++)
        dp += v1[i]*v2[i]; 
    return dp;
}

double distBetweenPoints(vector<double> &p1, vector<double> &p2)
{
    // Calculate the distance between points p1 and p2, |p1-p2|
    double sum=0.0;
    for (unsigned int i=0; i<3; i++)
        sum += pow((p1[i]-p2[i]),2.0);
    sum = sqrt(sum);
    return sum;
}

vector<double> crossProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates cross product v1xv2
    vector<double> cp(3);
    cp[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cp[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cp[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cp;
}

vector<double> vectAdd( vector<double> &a, vector<double> &b )
{
    // Vector addition, c=a+b
    return vectAdd(a, b, 1.0);
}

vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf )
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i] + sf*b[i];
    return c;
}
vector<double> vectSubtract( vector<double> &a, vector<double> &b )
{
    // Vector subtraction, c=a-b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i]-b[i];
    return c;
}

string NumberToString( int Number )
{
    // Converts integer to string
    ostringstream ss;
    ss << Number;
    return ss.str();
}

// ------------------------------------------------------
