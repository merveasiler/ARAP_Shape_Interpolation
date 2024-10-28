#include "pch.h"
#include "Mesh.h"

// Load the mesh in the OBJ file format
void Mesh::loadObj(const char* name)
{
	int nVerts, nTris;
	int lineNumber = 0;
	float x, y, z;
	string line;
	string cornerInfo1, cornerInfo2, cornerInfo3;
	int cornerId1, cornerId2, cornerId3;

	ifstream inpFile(name);
	if (inpFile.is_open()) {

		// Read the comments, detect the lineNumber of the last comment
		while (true) {
			getline(inpFile, line);
			if (line[0] != '#')
				break;
			lineNumber++;
		}

		// Initialize reading, get the number of vertices and triangles
		inpFile.seekg(0, inpFile.beg);
		for (int i = 0; i<lineNumber; i++)
			getline(inpFile, line);
		sscanf_s(line.c_str(), "%*s %*s %*s %d %*s %d", &nVerts, &nTris);

		// Read the vertices
		for (int i = 0; i < nVerts; i++) {
			inpFile >> line >> x >> y >> z;
			addVertex(x, y, z);
		}

		// Read the normals of the vertices
		for (int i = 0; i < nVerts; i++) {
			inpFile >> line >> x >> y >> z;
			this->getVertex(i)->setNormal(x, y, z);
		}

		// Read the triangles
		for (int i = 0; i < nTris; i++) {
			inpFile >> line >> cornerInfo1 >> cornerInfo2 >> cornerInfo3;
			sscanf_s(cornerInfo1.c_str(), "%d//", &cornerId1);
			sscanf_s(cornerInfo2.c_str(), "%d//", &cornerId2);
			sscanf_s(cornerInfo3.c_str(), "%d//", &cornerId3);
			addTriangle(cornerId1 - 1, cornerId2 - 1, cornerId3 - 1);	// indices start from 1 in the file
		}

		inpFile.close();
	}

}

// Load the mesh in the OFF file format
void Mesh::loadOff(const char* name) {

	int nVerts, nTris, dummy;
	float x, y, z;
	string line;
	int cornerId1, cornerId2, cornerId3;

	ifstream inpFile(name);
	if (inpFile.is_open()) {

		// Read OFF
		getline(inpFile, line);

		// Get the number of vertices and triangles (dummy is the number of edges)
		inpFile >> nVerts >> nTris >> dummy;

		// Read the vertices
		for (int i = 0; i < nVerts; i++) {
			inpFile >> x >> y >> z;
			addVertex(x, y, z);
		}

		// Read the triangles (dummy is the number of vertices)
		for (int i = 0; i < nTris; i++) {
			inpFile >> dummy >> cornerId1 >> cornerId2 >> cornerId3;
			addTriangle(cornerId1, cornerId2, cornerId3);
		}

		inpFile.close();
	}

}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back(new Vertex(idx, c));
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	int* c = new int[3];
	c[0] = v1;
	c[1] = v2;
	c[2] = v3;

	tris.push_back(new Triangle(idx, c));
	tris[idx]->computeCenter(verts[v1], verts[v2], verts[v3]);		// to be used in squish metric 
	tris[idx]->computeAreaVector(verts[v1], verts[v2], verts[v3]);	// to be used in squish metric

	//set up structure
	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	for (int i = 0; i<3; i++) {
		int corner1 = c[i];
		int corner2 = c[(i + 1) % 3];
		int edge_id = makeVertsNeighbor(corner1, corner2);
		if (edge_id < 0) {	// if the edge was not defined before
			edge_id = edges.size();
			addEdge(corner1, corner2);
		}
		else {	// construct adjacent triangle relationship
			for (int j = 0; j < edges[edge_id]->triList.size(); j++) {
				Triangle* adjTriangle = tris[edges[edge_id]->triList[j]];
				adjTriangle->triList.push_back(idx);
				tris[idx]->triList.push_back(adjTriangle->idx);
			}
		}

		tris[idx]->edgeList.push_back(edge_id);
		edges[edge_id]->triList.push_back(idx);
	}

}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();
	int* c = new int[2];
	c[0] = v1;
	c[1] = v2;

	edges.push_back(new Edge(idx, c));
	edges[idx]->computeLength(verts[v1], verts[v2]);

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

int Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	for (unsigned int i = 0; i < verts[v1i]->edgeList.size(); i++) {
		Edge* e = edges[verts[v1i]->edgeList[i]];
		for (int j = 0; j < 2; j++)
			if (e->endVerts[j] == v2i)
				return e->idx;
	}

	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return -1;
}

Mesh::~Mesh()
{
	for (unsigned int i = 0; i<verts.size(); i++)
		delete verts[i];
	verts.clear();

	for (unsigned int i = 0; i<tris.size(); i++)
		delete tris[i];
	tris.clear();

	for (unsigned int i = 0; i<edges.size(); i++)
		delete edges[i];
	edges.clear();

	for (unsigned int i = 0; i<tets.size(); i++)
		delete tets[i];
	tets.clear();
}

