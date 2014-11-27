
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>

using namespace std;

class cOctNode {
public:

    int MAX_OCTNODE_OBJECTS;
    int NUM_BRANCHES_OCTNODE;    
    double size;
    vector<double> position;
    vector<cOctNode> branches;
    int level;
    int nid;
    vector<int> data;
    vector<double> low, upp;
    
    cOctNode();
    cOctNode(int _level, int _nid, vector<double> _position, double _size);
    ~cOctNode();
    bool isLeafNode();
    void addPoly(int _indx);
    cOctNode addNode(int _level, int _nid, vector<double> _position, double size);
    int numPolys();
    void getLowUppVerts();
    void setupConstants();
};

class cTri {
public:

    int label;
    vector<vector<double> > vertices;
    vector<double> N;
    double D;
    vector<double> lowVert, uppVert;
    cTri();
    cTri(int label, vector<vector<double> > _vertices);
    ~cTri();
    void getN();
    void getD();
    bool isInNode(cOctNode &node);
    bool isInNode(vector<vector<double> > &lowUppVerts);
    void getLowerVert();
    void getUpperVert();
};

class cOctree {
public:

    int MAX_OCTREE_LEVELS;
    int branchOffsets[8][3];
    cOctNode root;
    vector<vector<double> > vertexCoords3D;
    vector<vector<int> > polyConnectivity;
    vector<cTri> polyList;
   
    cOctree(vector<vector<double> > _vertexCoords3D, vector<vector<int> > _polyConnectivity);
    ~cOctree();    
    int numPolys();
    void insertPoly(cOctNode &node, cTri &poly);
    void insertPolys();
    void setupPolyList();
    vector<double> getPositionRoot();
    double getSizeRoot();
    void splitNodeAndReallocate(cOctNode &node);
};

class CLine {
public:
    vector<double> p0,p1,dir;
    CLine();
    CLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~CLine();
    void getDir();
    void getP1();
};

// Function prototypes
double dotProduct( vector<double> &v1, vector<double> &v2 );
double distBetweenPoints(vector<double> &p1, vector<double> &p2);
vector<double> crossProduct( vector<double> &v1, vector<double> &v2 );
vector<double> vectAdd( vector<double> &a, vector<double> &b);
vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf);
vector<double> vectSubtract( vector<double> &a, vector<double> &b );
