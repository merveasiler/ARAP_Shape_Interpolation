#include "pch.h"
#include "Mesh.h"

void Mesh::computeTrisAngles() {

	for (unsigned int i = 0; i < tris.size(); i++) {

		Triangle* triangle = tris[i];

		Edge* a = edges[triangle->edgeList[0]];
		Edge* b = edges[triangle->edgeList[1]];
		Edge* c = edges[triangle->edgeList[2]];

		// law of cosines
		double cos_a = (pow((b->length), 2) + pow((c->length), 2) - pow((a->length), 2)) / (2 * b->length * c->length);
		double cos_b = (pow((a->length), 2) + pow((c->length), 2) - pow((b->length), 2)) / (2 * a->length * c->length);
		double cos_c = (pow((a->length), 2) + pow((b->length), 2) - pow((c->length), 2)) / (2 * a->length * b->length);

		triangle->angleList.push_back(acos(cos_a));
		triangle->angleList.push_back(acos(cos_b));
		triangle->angleList.push_back(acos(cos_c));

	}

}

void Mesh::tetrahedralizeSurface() {

	// find each triangle's surface tetrahedron 
	for (int i = 0; i < tris.size(); i++) {
		
		//compute the peak point of the tetrahedron (the fourth point)
		Triangle* triangle = tris[i];
		float* peakCoords = new float[3];

		// get the pointers for the triangle corners
		Vertex* corners[3];
		for (int j = 0; j < 3; j++)
			corners[j] = verts[triangle->corners[j]];

		// find the center of the triangle
		triangle->computeCenter(corners[0], corners[1], corners[2]);

		// find the vector in the direction of the normal which will point to the peak point
		float edgeVector1[3], edgeVector2[3];
		for (int j = 0; j < 3; j++) {
			edgeVector1[j] = corners[1]->coords[j] - corners[0]->coords[j];
			edgeVector2[j] = corners[2]->coords[j] - corners[1]->coords[j];
		}
		float* dirVect = triangle->crossProduct(edgeVector1, edgeVector2);
		float scaleFactor = sqrt(triangle->computeLength(dirVect));
		for (int j = 0; j < 3; j++)
			dirVect[j] /= scaleFactor;

		// find the peak point
		for (int j = 0; j < 3; j++)
			peakCoords[j] = triangle->center[j] + dirVect[j];
		delete[] dirVect;

		Vertex* peak = new Vertex(i, peakCoords);
		Tetrahedron* tetrahedron = new Tetrahedron(i, peak);
		tets.push_back(tetrahedron);

	}

}

// Writes the tetrahedralized surface version of a mesh into a file 
void Mesh::writeSurfaceTetrahedralizedVersion(const char* name) {

	int nVerts = verts.size();
	int nTris = tris.size();

	ofstream offFile(name);
	if (offFile.is_open()) {

		// Initialize wrtiting, write "OFF", number of vertices and triangles to .off file
		offFile << "OFF" << endl;
		offFile << nVerts + nTris << " " << nTris * 4 << " " << "0" << endl;

		// Write the vertices into .off file
		for (int i = 0; i < nVerts; i++) {
			float* coords = verts[i]->coords;
			offFile << coords[0] << " " << coords[1] << " " << coords[2] << " " << endl;	// x, y, z coordinates of the i^{th} vertex
		}

		for (int i = 0; i < nTris; i++) {
			float* coords = tets[i]->peak->coords;
			offFile << coords[0] << " " << coords[1] << " " << coords[2] << " " << endl;
		}
		// Write the triangles into .off file
		for (int i = 0; i < nTris; i++) {
			int* corners = tris[i]->corners;
			offFile << "3" << " " << corners[0] << " " << corners[1] << " " << corners[2] << endl;
			offFile << "3" << " " << corners[0] << " " << corners[1] << " " << nVerts + i << endl;
			offFile << "3" << " " << corners[1] << " " << corners[2] << " " << nVerts + i << endl;
			offFile << "3" << " " << corners[2] << " " << corners[0] << " " << nVerts + i << endl;
		}

		offFile.close();
	}
}

void Mesh::computeVolume() {

	volume = 0;
	for (int i = 0; i < tris.size(); i++) {

		Triangle* triangle = tris[i];
		Vertex* corners[3];

		for (int j = 0; j < 3; j++)
			corners[j] = verts[triangle->corners[j]];
		float* a_cross_b = triangle->crossProduct(corners[0]->coords, corners[1]->coords);

		float a_cross_b_dot_c = 0;
		for (int j = 0; j < 3; j++)
			a_cross_b_dot_c += a_cross_b[j] * corners[2]->coords[j];

		if (a_cross_b_dot_c < 0)
			a_cross_b_dot_c *= -1;

		volume += a_cross_b_dot_c / 6;
	}

}

/*
void Mesh::loadOffFromMatrices() {

// Read the vertices
int nVerts = TV.size() / 3;
for (int i = 0; i < nVerts; i++)
addVertex(TV(i, 0), TV(i, 1), TV(i, 2));

// Read the triangles (dummy is the number of vertices)
int nTris = TF.size() / 3;
for (int i = 0; i < nTris; i++)
addTriangle(TF(i, 0), TF(i, 1), TF(i, 2));

}
*/