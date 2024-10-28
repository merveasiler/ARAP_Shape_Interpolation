#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>

using namespace std;

struct Vertex
{
	int idx;
	float * coords = NULL;
	float * normal = NULL;

	// adjacencies
	vector< int > vertList;
	vector< int > triList;
	vector< int > edgeList;

	Vertex(int i, float* c) : idx(i), coords(c) { normal = new float[3]; };
	void setNormal(float x, float y, float z) {
		normal[0] = x;
		normal[1] = y;
		normal[2] = z;
	};
	~Vertex() {
		delete[] coords;
		delete[] normal;
		vertList.clear();
		triList.clear();
		edgeList.clear();
	};
};

struct Edge
{
	int idx;
	int * endVerts = NULL;
	float length;

	// adjacency
	vector< int > triList;

	Edge(int i, int* c) : idx(i), endVerts(c) {};
	~Edge() {
		delete[] endVerts;
		triList.clear();
	};
	void computeLength(Vertex* v1, Vertex* v2) {
		float* coords1 = v1->coords;
		float* coords2 = v2->coords;
		length = sqrt(pow(coords1[0] - coords2[0], 2) + pow(coords1[1] - coords2[1], 2) + pow(coords1[2] - coords2[2], 2));
	};
};

struct Triangle
{
	int idx;
	int * corners = NULL;
	float* center = NULL;	// center of gravity
	float* normal = NULL;
	float* areaVect = NULL; // unit area vector
	vector< int > edgeList;
	vector< int > triList; // neighbor tris
	vector< double > angleList;
	double angularSkewness;
	double squish;

	Triangle(int i, int* c) : idx(i), corners(c) {};
	~Triangle() {
		delete[] corners;
		delete[] center;
		delete[] normal;
		delete[] areaVect;
		edgeList.clear();
		triList.clear();
		angleList.clear();
	};

	float* crossProduct(float* vect1, float* vect2) {
		float* resultVect = new float[3];
		resultVect[0] = vect1[1] * vect2[2] - vect1[2] * vect2[1];
		resultVect[1] = -vect1[0] * vect2[2] + vect1[2] * vect2[0];
		resultVect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];
		return resultVect;
	}

	float computeLength(float* vect) {
		return sqrt(pow(vect[0], 2) + pow(vect[1], 2) + pow(vect[2], 2));
	}

	void normalize(float* vect) {
		float length = computeLength(vect);
		for (int i = 0; i < 3; i++)
			vect[i] = vect[i] / length;
	}

	void computeCenter(Vertex* v1, Vertex* v2, Vertex* v3) {
		center = new float[3];
		for (int i = 0; i < 3; i++)
			center[i] = (v1->coords[i] + v2->coords[i] + v3->coords[i]) / 3;
	}

	void computeNormal(Vertex* v1, Vertex* v2, Vertex* v3) {

		float edgeVect1[3], edgeVect2[3];
		for (int i = 0; i < 3; i++) {
			edgeVect1[i] = v2->coords[i] - v1->coords[i];
			edgeVect2[i] = v3->coords[i] - v2->coords[i];
		}

		normal = crossProduct(edgeVect1, edgeVect2);
		normalize(normal);
	}

	void computeAreaVector(Vertex* v1, Vertex* v2, Vertex* v3) {
		
		float vect1[3], vect2[3];
		for (int i = 0; i < 3; i++)
			vect1[i] = v3->coords[i] - v2->coords[i];
		for (int i = 0; i < 3; i++)
			vect2[i] = v1->coords[i] - v2->coords[i];

		// normalize
		float length1 = sqrt(pow(vect1[0], 2) + pow(vect1[1], 2) + pow(vect1[2], 2));
		float length2 = sqrt(pow(vect2[0], 2) + pow(vect2[1], 2) + pow(vect2[2], 2));
		for (int i = 0; i < 3; i++) {
			vect1[i] = vect1[i] / length1;
			vect2[i] = vect2[i] / length2;
		}

		// cross product
		areaVect = new float[3];
		areaVect[0] = vect1[1] * vect2[2] - vect1[2] * vect2[1];
		areaVect[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
		areaVect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

		// normalize
		float length = sqrt(pow(areaVect[0], 2) + pow(areaVect[1], 2) + pow(areaVect[2], 2));
		for (int i = 0; i < 3; i++)
			areaVect[i] = areaVect[i] / length;

	}

	void setAngSkewness(double angSkewness) {
		angularSkewness = angSkewness;
	}

	double getAngSkewness() {
		return angularSkewness;
	}

	void setSquish(double squish) {
		this->squish = squish;
	}

	double getSquish() {
		return squish;
	}

};


struct Tetrahedron
{
	int idx;			// also the id of base triangle of the tetrahedron
	Vertex* peak;		// peak point of the tetrahedron
	double transformationDiff;

	Tetrahedron(int i, Vertex* p) : idx(i), peak(p) {};
	~Tetrahedron() {
		delete peak;
	};

	void setTransformationDiff(double transformationDiff) {
		this->transformationDiff = transformationDiff;
	}

	double getTransformationDiff() {
		return transformationDiff;
	}

};

class Mesh
{
	friend class Transformation;

private:
	float volume;

	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;
	vector< Tetrahedron* > tets;

	// methods to construct mesh
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	int makeVertsNeighbor(int v1i, int v2i);

public:
	Mesh() {};
	~Mesh();
	void loadObj(const char* name);
	void loadOff(const char* name);

	// methods to reach mesh elements
	int getNumOfVerts() { return verts.size(); };
	int getNumOfTris() { return tris.size(); };
	int getNumOfEdges() { return edges.size(); };
	int getNumOfTets() { return tets.size(); };
	Vertex* getVertex(int i) { return verts[i]; };
	Triangle* getTriangle(int i) { return tris[i]; };
	Edge* getEdge(int i) { return edges[i]; };
	Tetrahedron* getTetrahedron(int i) { return tets[i]; };
	const vector<Vertex*> getAllVerts() { return verts; };
	const vector<Triangle*> getAllTris() { return tris; };
	const vector<Edge*> getAllEdges() { return edges; };
	const vector<Tetrahedron*> getAllTetrahedra() { return tets; };

	// methods to compute mesh features
	void computeTrisAngles();
	void computeVolume();
	float getVolume() { return volume;  };
	void tetrahedralizeSurface();
	void writeSurfaceTetrahedralizedVersion(const char* name);
	void loadOffFromMatrices();
};

#pragma once
