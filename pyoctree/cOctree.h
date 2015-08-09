
// Copyright (C) 2015 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm> // find, sort
#include <utility>   // pair

// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

typedef struct intersection
{
    double s;
    vector<double> p;
    intersection() { p.resize(3,0.0); s=0.0; }
    intersection(vector<double> _p, double _s) { p=_p; s=_s;}
    // Overload operator < so we can use sort function of a std::vector
    bool operator < (const intersection& intersect) const { 
        return (s < intersect.s); }
} Intersection;

class cLine {
public:

    vector<double> p0,p1,dir;
    cLine();
    cLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~cLine();
    void getDir();
    void getP1();
};

class cOctNode {
public:

    static const int MAX_OCTNODE_OBJECTS  = 200;
    static const int NUM_BRANCHES_OCTNODE = 8;  
    double size;
    int level;
    string nid;
    vector<double> position;
    vector<cOctNode> branches;
    vector<int> data;
    vector<double> low, upp;
    cOctNode();
    cOctNode(int _level, string _nid, vector<double> _position, double _size);
    ~cOctNode();
    bool isLeafNode();
    int numPolys();
    void addPoly(int _indx);
    void addNode(int _level, string _nid, vector<double> _position, double _size);
    void getLowUppVerts();
    bool boxRayIntersect(cLine &ray);
    bool sphereRayIntersect(cLine &ray);
};

class cTri {
public:

    double D;
    int label;
    vector<vector<double> > vertices;
    vector<double> N;
    vector<double> lowVert, uppVert;
    cTri();
    cTri(int _label, vector<vector<double> > _vertices);
    ~cTri();
    bool isInNode(cOctNode &node);
    bool isPointInTri(vector<double> &p);
    bool rayPlaneIntersectPoint(cLine &ray, bool entryOnly);
    bool rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s);
    void getN();
    void getD();
    void getLowerVert();
    void getUpperVert();
};

class cOctree {
public:

    static const int MAX_OCTREE_LEVELS = 10;
    int branchOffsets[8][3];
    cOctNode root;
    vector<vector<double> > vertexCoords3D;
    vector<vector<int> > polyConnectivity;
    vector<cTri> polyList;
    cOctree(vector<vector<double> > _vertexCoords3D, vector<vector<int> > _polyConnectivity);
    ~cOctree();    
    double getSizeRoot();
    int numPolys();
    cOctNode* getNodeFromId(string nodeId);
    cOctNode* findBranchById(string nodeId, cOctNode &node);
    set<int> getListPolysToCheck(cLine &ray);    
    vector<double> getPositionRoot();	
    vector<Intersection> findRayIntersect(cLine &ray);    
    vector<int> findRayIntersects(vector<cLine> &rayList);
    vector<int> findRayIntersectsSorted(vector<cLine> &rayList);		
    vector<cOctNode*> getNodesFromLabel(int polyLabel);	
    vector<cOctNode*> getSortedNodesToCheck(cLine &ray);
    void insertPoly(cOctNode &node, cTri &poly);
    void insertPolys();
    void setupPolyList();
    void splitNodeAndReallocate(cOctNode &node);
    void findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList);
    void getPolysToCheck(cOctNode &node, cLine &ray, set<int> &intTestPolys);
    void getNodesToCheck(cOctNode &node, cLine &ray, vector<pair<cOctNode*,double> > &nodeList);
};

// Function prototypes
bool sortNodes(const pair<cOctNode*,double>&i, const pair<cOctNode*,double>&j);
double dotProduct( vector<double> &v1, vector<double> &v2 );
double distBetweenPoints(vector<double> &p1, vector<double> &p2);
string NumberToString( int Number );
vector<double> crossProduct( vector<double> &v1, vector<double> &v2 );
vector<double> vectAdd( vector<double> &a, vector<double> &b);
vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf);
vector<double> vectSubtract( vector<double> &a, vector<double> &b );
