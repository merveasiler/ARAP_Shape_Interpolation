#include "pch.h"
#include "Transformation.h"

/**
* @param tetrahedron;
* @param mesh; shape that the tetrahedron belongs to
* @param V; will include the resulting matrix where each column is a vector from the peak to one of the corners of tetrahedron
*/
void Transformation::computeMatrixOfVectorsFromPeakToCorners(Tetrahedron* tetrahedron, Mesh* mesh, Matrix3d & V) {

	// fetch the vertices of the tetrahedron
	Vertex * verts[4];

	for (int i = 0; i < 3; i++)
		verts[i] = mesh->getVertex(mesh->getTriangle(tetrahedron->idx)->corners[i]);
	verts[3] = tetrahedron->peak;

	// Compute matrix V = [v1-v4; v2-v4; v3-v4]
	for (int i = 0; i < 3; i++)	// p0, p1, p2 vertices of the tetrahedron
		for (int j = 0; j < 3; j++)	// x, y, z coordinates
			V(j, i) = verts[i]->coords[j] - verts[3]->coords[j];

}

/**
* @brief Computes the transformation matrix from a source tetrahedron to a target tetrahedron.

* @param sourceTet; source tetrahedron
* @param targetTet; target tetrahedron
* @param M; will include the resulting 3x3 transformation matrix consisting of rotation and scaling
*/
void Transformation::transformTetToTet(Tetrahedron* sourceTet, Tetrahedron* targetTet, MatrixXd & M) {

	// Find the affine transformation matrix M = P * Q^{-1}
	Matrix3d P, Q;

	// compute P and Q
	computeMatrixOfVectorsFromPeakToCorners(sourceTet, sourceMesh, P);
	computeMatrixOfVectorsFromPeakToCorners(targetTet, targetMesh, Q);

	/*	// REMOVE
	for (int i = 0; i < 3; i++) {
	Triangle* t = sourceMesh->getTriangle(0);
	Vertex* v0 = sourceMesh->getVertex(t->corners[0]);
	Vertex* v1 = sourceMesh->getVertex(t->corners[1]);
	Vertex* v2 = sourceMesh->getVertex(t->corners[2]);
	Vertex* v3 = sourceTet->peak;
	cout << v0->coords[i] << " " << v1->coords[i] << " " << v2->coords[i] << " " << v3->coords[i] << endl;
	}

	cout << "source peak to corners" << endl;
	cout << P << endl;
	cout << "target peak to corners" << endl;
	cout << Q << endl;
	*/

	// compute M
	M = Q * P.inverse();

}

/**
* @brief Computes the transformation matrix from a source tetrahedron to a target tetrahedron.

* @param M; the transformation matrix including combined rotation and scaling
* @param R; will include the resulting 3x3 rotation matrix of the transformation
* @param S; will include the resulting 3x3 scaling matrix of the transformation
*/
void Transformation::factorizeTransformationAsRandS(MatrixXd & M, Matrix3d & R, Matrix3d & S) {

	// apply SVD (singular value decomposition) to M = U * D * V^{t}
	JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);

	/* Notes:
	-> M = U * (V^{t} * V) * D * V^{t}
	-> R = U * V^{t} and S = V * D * V^{t}
	-> M = R * S
	*/

	// Compute R and S
	R = svd.matrixU() * svd.matrixV().transpose();
	S = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();

}

/**
* @brief Computes the translation vector from a source tetrahedron to a target tetrahedron.

* @param sourceVert; one of the vertices of the source tetrahedron
* @param targetVert; the corresponding vertex of the target tetrahedron
* @param M; the transformation matrix including combined rotation and scaling
* @param T; will include the resulting 3x1 translation vector of the transformation
*/
void Transformation::computeTranslation(Vertex* sourceVert, Vertex* targetVert, const Matrix3d & M, Vector3d & T) {

	Vector3d sourceVec, targetVec;

	for (int i = 0; i < 3; i++) {
		sourceVec(i) = sourceVert->coords[i];
		targetVec(i) = targetVert->coords[i];
	}

	T = targetVec - M * sourceVec;

}

/**
* @brief Applies linear interpolation from a source matrix to a target matrix in [0,1] and finds the matrix at time t.

* @param sourceMatrix
* @param targetMatrix
* @param lambda; some moment in [0,1]
* @param resultMatrix; will include the resulting 3x3 matrix at the time @lambda
*/
void Transformation::interpolateLinearly(const Matrix3d & sourceMatrix, const Matrix3d & targetMatrix, float lambda, Matrix3d & resultMatrix) {

	resultMatrix = sourceMatrix * (1 - lambda) + targetMatrix * (lambda);

}

/**
* @brief Applies quaternion interpolation from a source matrix to a target matrix in [0,1] and finds the matrix at time t.

* @param sourceMatrix
* @param targetMatrix
* @param lambda; some moment in [0,1]
* @param resultMatrix; will include the resulting 3x3 matrix at the time @lambda
*/
void Transformation::interpolateByQuaternions(const Matrix3d & sourceMatrix, const Matrix3d & targetMatrix, float lambda, Matrix3d & resultMatrix) {

	Quaterniond qSource, qTarget;
	qSource = sourceMatrix;										// convertMatrixToQuaternion(sourceMatrix);
	qTarget = targetMatrix;										// convertMatrixToQuaternion(targetMatrix);

	Quaterniond q_at_lambda = qSource.slerp(lambda, qTarget);	// slerp(qSource, qTarget, lambda);

	resultMatrix = q_at_lambda.toRotationMatrix();				// convertQuaternionToMatrix(q_at_lambda);

}

/**
*
*/
void Transformation::interpolateTetrahedron(Tetrahedron* sourceTet, Tetrahedron* targetTet, float lambda) {

	// Compute the transformation matrix from sourceTet to targetTet as a multiplication of rotating and scaling part
	MatrixXd M;
	transformTetToTet(sourceTet, targetTet, M);
	Matrix3d R, S;
	factorizeTransformationAsRandS(M, R, S);

	// Compute the translation by using the difference between the transformed and original meshes
	Vector3d T, T_at_lambda;
	computeTranslation(sourceTet->peak, targetTet->peak, M, T);
	T_at_lambda = T * lambda;												// Interpolate the translation part linearly

																			// Interpolate the rotation and scaling matrices
	Matrix3d M_at_lambda, R_at_lambda, S_at_lambda;
	interpolateLinearly(Matrix3d::Identity(), S, lambda, S_at_lambda);		// interpolate the scaling part linearly
	interpolateByQuaternions(Matrix3d::Identity(), R, lambda, R_at_lambda);	// interpolate the rotating part by quaternions
	M_at_lambda = R_at_lambda * S_at_lambda;

	// Construct the system of linear equations
	Matrix3d P, P_inverse;
	computeMatrixOfVectorsFromPeakToCorners(sourceTet, sourceMesh, P);
	P_inverse = P.inverse();	// [ U_alpha | U_beta | U_gama ]

								// Compute the coefficient matrices
	int tetrahedron_id = sourceTet->idx;
	fillVectorKs(tetrahedron_id, P_inverse);	// from transformation part(rotation + scaling)
	fillVectorC(tetrahedron_id, P_inverse);		// from translation part

												// Add the new coefficient to the system of linear equations
												// Prepare the vertex indices (the corresponding row/column ids in the coefficients)
	int indices[4];
	for (int i = 0; i < 3; i++)
		indices[i] = sourceMesh->getTriangle(tetrahedron_id)->corners[i];
	indices[3] = sourceMesh->getNumOfVerts() + tetrahedron_id;
	cout << tetrahedron_id << endl;
	// Smart sparse matrix addition for A += K * K.transpose() + C * C.transpose();
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			int ind1 = indices[i];
			int ind2 = indices[j];

			double value = 0;
			for (int k = 0; k < 3; k++)
				value += K.coeffRef(ind1, k) * K.coeffRef(ind2, k);

			A.coeffRef(ind1, ind2) += value + C.coeffRef(ind1) * C.coeffRef(ind2);
		}

	// Smart sparse vector addition for B.col(i) += K * M_at_lambda.row(i).transpose() + T_at_lambda(i) * C;
	for (int i = 0; i < 4; i++) {
		int ind1 = indices[i];
		for (int j = 0; j < 3; j++) {	// for each of x, y, z coordinates

			double value = 0;
			for (int k = 0; k < 3; k++)
				value += K.coeffRef(ind1, k) * M_at_lambda(j, k);

			B(ind1, j) += value + T_at_lambda(j) * C.coeffRef(ind1);
		}
	}

	/* ORGINAL VERSION
	A += K * K.transpose() + C * C.transpose();
	for (int i = 0; i < 3; i++)		// b1, b2, b3 for x, y, z solutions resp.
	B.col(i) += K * M_at_lambda.row(i).transpose() + T_at_lambda(i) * C;
	*/
}

/**
*
*/
void Transformation::interpolateARAPByAlexa(float lambda, string resultMeshName) {

	time_t begin_time, end_time, first_time;

	// Tetrahedralize the surfaces of the triangles
	time(&first_time);
	time(&begin_time);
	sourceMesh->tetrahedralizeSurface();
	targetMesh->tetrahedralizeSurface();
	time(&end_time);
	cout << "Elapsed time in tetrahedralization of surfaces: " << difftime(end_time, begin_time) << " seconds\n";

	// Construct the system of linear equations by computing the transformation over each tetrahedron
	time(&begin_time);
	for (int i = 0; i < sourceMesh->getNumOfTets(); i++)
		interpolateTetrahedron(sourceMesh->getTetrahedron(i), targetMesh->getTetrahedron(i), lambda);
	time(&end_time);
	cout << "Elapsed time in the interpolation computations: " << difftime(end_time, begin_time) << " seconds\n";

	time(&begin_time);
	// Solve the system of linear equations
	SimplicialLDLT < SparseMatrix<double> > solver(A);	//solver.compute(A);
	if (solver.info() != Success) {
		cout << "Could not decompose A.\n";
		return;
	}
	else
		cout << "A was decomposed successfully.\n";

	for (int i = 0; i < 3; i++) {
		X.col(i) = solver.solve(B.col(i));
		if (solver.info() != Success) {
			cout << "Could not solve the sparse linear system.\n";
			return;
		}
	}

	time(&end_time);
	cout << "Elapsed time in solving the linear system: " << difftime(end_time, begin_time) << " seconds\n";
	cout << "Total elapsed time in the algorithm: " << difftime(end_time, first_time) << " seconds\n";
	cout << "Yeahh, the sparse linear system was solved.\n";

	writeInterpolatedMesh(resultMeshName);

}

void Transformation::writeInterpolatedMesh(string resultMeshName) {

	int nVerts = sourceMesh->getNumOfVerts(), nTris = sourceMesh->getNumOfTris();

	ofstream offFile(resultMeshName);
	if (offFile.is_open()) {

		// Initialize wrtiting, write "OFF", number of vertices and triangles to .off file
		offFile << "OFF" << endl;
		offFile << nVerts << " " << nTris << " " << "0" << endl;

		// Write the vertices into .off file
		for (int i = 0; i < nVerts; i++)
			offFile << X(i, 0) << " " << X(i, 1) << " " << X(i, 2) << " " << endl;	// x, y, z coordinates of the i^{th} vertex

																					// Write the triangles into .off file
		for (int i = 0; i < nTris; i++) {
			int* corners = sourceMesh->getTriangle(i)->corners;
			offFile << "3" << " " << corners[0] << " " << corners[1] << " " << corners[2] << endl;
		}

		offFile.close();

	}
}

void Transformation::fillVectorKs(int triangle_id, const Matrix3d & P_inverse) {

	int n = sourceMesh->getNumOfVerts();
	int* triangleCorners = sourceMesh->getTriangle(triangle_id)->corners;

	K.setZero();

	for (int j = 0; j < 3; j++)
		for (int i = 0; i < 3; i++) {
			double value = P_inverse(i, j);
			K.coeffRef(triangleCorners[i], j) = value;
			K.coeffRef(n + triangle_id, j) -= value;
		}

}

void Transformation::fillVectorC(int triangle_id, const Matrix3d & P_inverse) {

	Vector3d sourcePeak;
	for (int i = 0; i < 3; i++)
		sourcePeak(i) = sourceMesh->getTetrahedron(triangle_id)->peak->coords[i];

	Vector3d R = P_inverse * sourcePeak;

	int n = sourceMesh->getNumOfVerts();
	int* triangleCorners = sourceMesh->getTriangle(triangle_id)->corners;

	C.setZero();
	C.coeffRef(n + triangle_id) = 1;

	for (int i = 0; i < 3; i++) {
		double value = R(i);
		C.coeffRef(triangleCorners[i]) = -value;
		C.coeffRef(n + triangle_id) += value;
	}

}

/*
void Transformation::fillVectorKs(SparseMatrix<double> K, int triangle_id, const Matrix3d & P_inverse) {

int n = sourceMesh->getNumOfVerts();
int* triangleCorners = sourceMesh->getTriangle(triangle_id)->corners;

K.setZero();

for (int j = 0; j < 3; j++) {
double peakCoeff = 0;
for (int i = 0; i < 3; i++) {
double value = P_inverse(i, j);
K.row(triangleCorners[i]).col(j) = value;
peakCoeff -= value;
}
K.row(n + triangle_id).col(j) = peakCoeff;
}

}

void Transformation::fillVectorC(SparseVector<double> C, int triangle_id, const Matrix3d & P_inverse) {

Vector3d sourcePeak;
for (int i = 0; i < 3; i++)
sourcePeak(i) = sourceMesh->getTetrahedron(triangle_id)->peak->coords[i];

Vector3d R = P_inverse * sourcePeak;

int n = sourceMesh->getNumOfVerts();
int* triangleCorners = sourceMesh->getTriangle(triangle_id)->corners;

C.setZero();

double peakCoeff = 1;
for (int i = 0; i < 3; i++) {
double value = R(i);
C.row(triangleCorners[i]) = -value;
peakCoeff += value;
}
C.row(n + triangle_id) = peakCoeff;

}
*/
/* ****************************************************** */
/* MY OWN QUATERNION RELATED IMPLEMENTATIONS (not in use) */
/* ****************************************************** */

/**
* @brief Manages the interpolation of rotation matrixs by using quaternion interpolation.
* @param startPoint; 3x1 array representing the vector part of the source quaternion
* @param endPoint; 3x1 array representing the vector part of the target quaternion
* @param lambda; some moment in [0, 1]
* @return the corresponding 3x3 rotation matrix at time @lambda through the interpolation
*/
Matrix3d & Transformation::interpolateRotation(double startPoint[3], double endPoint[3], float lambda) {

	Quaterniond qStart(0, startPoint[0], startPoint[1], startPoint[2]);
	Quaterniond qEnd(0, endPoint[0], endPoint[1], endPoint[2]);

	qStart.normalize();
	qEnd.normalize();

	Quaterniond qRotation = slerp(qStart, qEnd, lambda);
	return convertQuaternionToMatrix(qRotation);

}

/**
* @brief Converts a rotation matrix to a quaternion.
* @param R; 3x3 rotation matrix
* @result the corresponding float quaternion
*/
Quaterniond & Transformation::convertMatrixToQuaternion(const Matrix3d & R) {

	// q = [w,v] where v = (x, y, z)^{t}
	double w, x, y, z;
	double w_square = (1 + R(0, 0) + R(1, 1) + R(2, 2)) / 4;

	if (w_square > DBL_EPSILON) {
		w = sqrtf(w_square);
		double w_times_4 = 4 * w;
		x = (R(2, 1) - R(1, 2)) / w_times_4;
		y = (R(0, 2) - R(2, 0)) / w_times_4;
		z = (R(1, 0) - R(0, 1)) / w_times_4;
	}
	else {
		w = 0;
		double x_square = (R(1, 1) + R(2, 2)) / -2;
		if (x_square > DBL_EPSILON) {
			x = sqrtf(x_square);
			double x_times_2 = 2 * x;
			y = R(1, 0) / x_times_2;
			z = R(2, 0) / x_times_2;
		}
		else {
			x = 0;
			double y_square = (1 - R(2, 2)) / 2;
			if (y_square > DBL_EPSILON) {
				y = sqrtf(y_square);
				z = R(2, 1) / (2 * y);
			}
			else {
				y = 0;
				z = 1;
			}
		}
	}

	Quaterniond q(w, x, y, z);
	q.normalize();
	// cout << "my conversion:" << endl << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << endl;
	return q;
}

/**
* @brief Converts a quaternion to a rotation matrix.
* @param q; float quaternion
* @result the corresponding 3x3 rotation matrix
*/
Matrix3d & Transformation::convertQuaternionToMatrix(const Quaterniond & q) {

	Matrix3d M;
	double w = q.w();
	double x = q.vec().x();
	double y = q.vec().y();
	double z = q.vec().z();

	M(0, 0) = 1 - 2 * pow(y, 2) - 2 * pow(z, 2);
	M(1, 0) = 2 * x*y + 2 * w*z;
	M(2, 0) = 2 * x*z - 2 * w*y;
	M(0, 1) = 2 * x*y - 2 * w*z;
	M(1, 1) = 1 - 2 * pow(x, 2) - 2 * pow(z, 2);
	M(2, 1) = 2 * y*z + 2 * w*x;
	M(0, 2) = 2 * x*z + 2 * w*y;
	M(1, 2) = 2 * y*z - 2 * w*x;
	M(2, 2) = 1 - 2 * pow(x, 2) - 2 * pow(y, 2);

	// cout << "my conversion: " << endl << M << endl;
	return M;
}

/**
* @brief Applies Spherical Linear Intepolation (SLERP) from a quaternion to the other and computes the quaternion at some time t.
* @param q1; source quaternion
* @param q2; target quaternion
* @param lambda; some moment in [0, 1]
* @return the corresponding quaternion at time @lambda through the interpolation from @q1 to @q2
*/
Quaterniond & Transformation::slerp(const Quaterniond & q1, const Quaterniond & q2, float lambda)
{
	/* q1 and q2 is actually the start and end points of the object to be rotated.
	Also they should be normalized. */

	double cosTheta = q1.dot(q2);
	double theta = acos(cosTheta);
	double sinTheta = sin(theta);

	double coeff1 = sin((1 - lambda)*theta) / sinTheta;
	double coeff2 = sin(lambda*theta) / sinTheta;

	Quaterniond q;
	q.w() = coeff1 * q1.w() + coeff2 * q2.w();
	q.x() = coeff1 * q1.x() + coeff2 * q2.x();
	q.y() = coeff1 * q1.y() + coeff2 * q2.y();
	q.z() = coeff1 * q1.z() + coeff2 * q2.z();

	q.normalize();
	// cout << "my slerp:" << endl << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << endl;
	return q;
}

// Do experiment of Sparse Matrix Linear System Solver
void doExperimentOfLinearSystemSolver() {

	SparseMatrix<double> tempA;	tempA.resize(3, 3);	tempA.reserve(3); tempA.setZero();
	MatrixXd tempB;	tempB.resize(3, 2);	tempB.setZero();
	MatrixXd tempX;	tempX.resize(3, 2);

	tempA.coeffRef(0, 0) = 3;	tempA.coeffRef(1, 1) = 2;	tempA.coeffRef(2, 2) = 7;
	tempB(0, 0) = 39;	tempB(1, 0) = 12;	tempB(2, 0) = 91;
	tempB(0, 1) = -3;	tempB(1, 1) = 32;	tempB(2, 1) = -21;

	cout << "Sparse Matrix A :" << endl;
	cout << tempA << endl << endl;
	cout << "Result Matrix B :" << endl;
	cout << tempB << endl << endl;

	SimplicialLDLT < SparseMatrix<double> > solver(tempA);
	for (int i = 0; i < 2; i++) {
		tempX.col(i) = solver.solve(tempB.col(i));
		if (solver.info() != Success) {
			cout << "Could not solve the sparse linear system.\n";
			return;
		}
	}

	cout << "Variable Matrix X:" << endl;
	cout << tempX << endl;
}
