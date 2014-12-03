
#include "cOctree.h"

// ------------------------------------------------------

CLine::CLine() 
{
    // Default line is unit vector along x-axis
    p0.resize(3,0.0); p1.resize(3,0.0); dir.resize(3,0.0);
    p1[0]=1.0; dir[0]=1.0;
}

CLine::CLine(vector<double> &lp0, vector<double> &p1_dir, int isP1orDir)
{
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = lp0;
    if (isP1orDir==0) {
        p1 = p1_dir;
        getDir(); }
    else if (isP1orDir==1) {
        dir = p1_dir;
        getP1(); }
}

CLine::~CLine() {}

void CLine::getDir()
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

void CLine::getP1()
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

bool cTri::isInNode(vector<vector<double> > &lowUppVerts)
{
    vector<double> nodelow, nodeupp;
    nodelow = lowUppVerts[0];
    nodeupp = lowUppVerts[1];
    if (lowVert[0] > nodeupp[0]) return false;
    if (lowVert[1] > nodeupp[1]) return false;
    if (lowVert[2] > nodeupp[2]) return false;
    if (uppVert[0] < nodelow[0]) return false;
    if (uppVert[1] < nodelow[1]) return false;
    if (uppVert[2] < nodelow[2]) return false;
    return true;
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
    MAX_OCTNODE_OBJECTS  = 500;   
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

// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------

/*
cOctNode* cOctree::getNodeFromLabel(int polyLabel)
{
    return findBranchByLabel(polyLabel,root);
}

cOctNode* cOctree::findBranchByLabel(int polyLabel, cOctNode &node)
{
    if (node.isLeafNode()) {
	    vector<int>::iterator it;
	    it = find(node.data.begin(),node.data.end(),polyLabel);
		if (it != node.data.end()) { 
		    return &node; }
	} else {
	    for (unsigned int i=0; i<node.branches.size(); i++) {
		    cOctNode *branch = findBranchByLabel(polyLabel, node.branches[i]);
			if (branch != NULL) { return branch; }
		}
	}
	return NULL;
}
*/

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
