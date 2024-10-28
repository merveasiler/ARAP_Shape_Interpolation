#include "pch.h"

/*
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/barycenter.h>
#include <Eigen/Core>
#include <string>

using namespace Eigen;
using namespace std;

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;


void tetrahedralizeMesh(string name) {

	// Load a surface mesh
	igl::readOFF(name, V, F);

	// Tetrahedralize the interior		pq1.414Y
	igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414", TV, TT, TF);

	// Write the surface mesh
	igl::writeOFF(name + "_tetrahedralized.off", TV, TF);
}
*/