#pragma once

#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Geometry>
#include <Eigen/SparseCholesky>
#include <ctime>

using namespace Eigen;
using namespace std;

class Transformation {

	Mesh* sourceMesh;
	Mesh* targetMesh;

	// Declare the vectors and matrices in the system of linear equations
	int N;					// Number of unknown variables
	SparseMatrix<double> A;	// A*X = b1, A*Y = b2, A *Z = b3
	MatrixXd B;				// b1, b2, b3 for x, y, z solutions resp.
	MatrixXd X;				// X = [x, y, z] component vectors

							// Declare the helper vector and matrices
	SparseMatrix<double> K;	// coefficients k1, k2, k3 from the rotation + scaling part
	SparseVector<double> C;

public:
	Transformation(Mesh* sMesh, Mesh* tMesh) : sourceMesh(sMesh), targetMesh(tMesh) {
		N = sourceMesh->getNumOfVerts() + sourceMesh->getNumOfTris();

		A.resize(N, N);
		//A.reserve(N * 4);
		A.setZero();

		B.resize(N, 3);
		B.setZero();

		X.resize(N, 3);

		K.resize(N, 3);
		//K.reserve(12);
		C.resize(N, 1);
		//C.reserve(4);

	};

	void computeMatrixOfVectorsFromPeakToCorners(Tetrahedron* tetrahedron, Mesh* mesh, Matrix3d & V);
	void transformTetToTet(Tetrahedron* sourceTet, Tetrahedron* targetTet, MatrixXd & M);
	void factorizeTransformationAsRandS(MatrixXd & M, Matrix3d & R, Matrix3d & S);
	void computeTranslation(Vertex* sourceVert, Vertex* targetVert, const Matrix3d & M, Vector3d & T);
	void interpolateLinearly(const Matrix3d & sourceMatrix, const Matrix3d & targetMatrix, float lambda, Matrix3d & resultMatrix);
	void interpolateByQuaternions(const Matrix3d & sourceMatrix, const Matrix3d & targetMatrix, float lambda, Matrix3d & resultMatrix);

	void fillVectorKs(int triangle_id, const Matrix3d & P_inverse);
	void fillVectorC(int triangle_id, const Matrix3d & P_inverse);

	void interpolateTetrahedron(Tetrahedron* sourceTet, Tetrahedron* targetTet, float lambda);
	void interpolateARAPByAlexa(float lambda, string resultMeshName);
	void writeInterpolatedMesh(string resultMeshName);

	/* ****************************************************** */
	/* MY OWN QUATERNION RELATED IMPLEMENTATIONS (not in use) */
	/* ****************************************************** */
	Matrix3d & interpolateRotation(double startPoint[3], double endPoint[3], float lambda);
	Quaterniond & convertMatrixToQuaternion(const Matrix3d & R);
	Matrix3d & convertQuaternionToMatrix(const Quaterniond & q);
	Quaterniond & slerp(const Quaterniond & q1, const Quaterniond & q2, float lambda);

};