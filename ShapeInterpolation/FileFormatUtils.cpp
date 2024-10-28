#include "pch.h"

#include <string>
#include <fstream>
using namespace std;


/* Take the name of the obj file,
 * Create the off format of that file
 */
void convertObjToOff(string name) {
	int nVerts, nTris;
	int lineNumber = 0;
	int dummy = 0;
	float x, y, z;
	string line;
	string cornerInfo1, cornerInfo2, cornerInfo3;
	int cornerId1, cornerId2, cornerId3;

	ifstream objFile(name.c_str());
	ofstream offFile((name + ".off").c_str());

	if (objFile.is_open() && offFile.is_open()) {

		// Read the comments, detect the lineNumber of the last comment
		while (true) {
			getline(objFile, line);
			if (line[0] != '#')
				break;
			lineNumber++;
		}

		// Initialize reading, get the number of vertices and triangles from .obj file
		objFile.seekg(0, objFile.beg);
		for (int i = 0; i < lineNumber; i++)
			getline(objFile, line);
		sscanf_s(line.c_str(), "%*s %*s %*s %d %*s %d", &nVerts, &nTris);

		// Initialize wrtiting, write "OFF", number of vertices and triangles to .off file
		offFile << "OFF" << endl;
		offFile << nVerts << " " << nTris << " " << "0" << endl;

		// Read the vertices from .obj file and write into .off file
		for (int i = 0; i < nVerts; i++) {
			objFile >> line >> x >> y >> z;
			offFile << x << " " << y << " " << z << " " << endl;
		}

		// Read the normals of the vertices from .obj file
		for (int i = 0; i < nVerts; i++)
			objFile >> line >> x >> y >> z;

		// Read the triangles from .obj file, write into .off file
		for (int i = 0; i < nTris; i++) {
			objFile >> line >> cornerInfo1 >> cornerInfo2 >> cornerInfo3;
			sscanf_s(cornerInfo1.c_str(), "%d//", &cornerId1);
			sscanf_s(cornerInfo2.c_str(), "%d//", &cornerId2);
			sscanf_s(cornerInfo3.c_str(), "%d//", &cornerId3);
			// indices start from 1 in the .obj file:
			offFile << "3" << " " << cornerId1 - 1 << " " << cornerId2 - 1 << " " << cornerId3 - 1 << endl;
		}

		objFile.close();
		offFile.close();
	}
}
