#include "pch.h"
#include "Metric.h"

template <class Type>
Type findMax(const vector<Type> elementList) {

	Type maximum = elementList[0];
	for (unsigned int i = 1; i<elementList.size(); i++)
		maximum = (elementList[i]>maximum) ? elementList[i] : maximum;
	return maximum;

}

template <class Type>
Type findMin(const vector<Type> elementList) {

	Type minimum = elementList[0];
	for (unsigned int i = 1; i<elementList.size(); i++)
		minimum = (elementList[i]<minimum) ? elementList[i] : minimum;
	return minimum;

}

double measureTriangleSkewness(Triangle* triangle) {

	double maxAngle = findMax(triangle->angleList);
	double minAngle = findMin(triangle->angleList);
	double optimalAngle = _PI / 3;

	double ratio1 = (maxAngle - optimalAngle) / (_PI - optimalAngle);
	double ratio2 = (optimalAngle - minAngle) / (optimalAngle);
	double maxRatio = (ratio1 > ratio2) ? ratio1 : ratio2;

	triangle->setAngSkewness(maxRatio);
	return maxRatio;

}

void measureBySkewness(Mesh* mesh) {

	mesh->computeTrisAngles();

	// calculate the average of skewness values
	double avgSkewness = 0;
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		avgSkewness += measureTriangleSkewness(mesh->getTriangle(i));
	avgSkewness = avgSkewness / mesh->getNumOfTris();
	//cout << avgSkewness << " ";

	// calculate the standard deviation of skewness values
	double stdDeviation = 0;
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		stdDeviation += pow(mesh->getTriangle(i)->getAngSkewness() - avgSkewness, 2);
	stdDeviation = sqrt(stdDeviation / mesh->getNumOfTris());
	cout << stdDeviation << " ";

}

double measureTriangleSquish(Triangle* triangle, Mesh* mesh) {

	double maxSquish = 0;
	float centerToCenterVect[3];

	// find the neighbor triangles by 
	for (int i = 0; i < triangle->triList.size(); i++) {
		Triangle* adjTriangle = mesh->getTriangle(triangle->triList[i]);
		//compute the vector from the current triangle's center to the adjacent triangle's center
		for (int j = 0; j < 3; j++)
			centerToCenterVect[j] = adjTriangle->center[j] - triangle->center[j];
		// compute the squish value
		double length = sqrt(pow(centerToCenterVect[0], 2) + pow(centerToCenterVect[1], 2) + pow(centerToCenterVect[2], 2));
		for (int j = 0; j < 3; j++)
			centerToCenterVect[j] = centerToCenterVect[j] / length;
		double squish = 0;
		for (int j = 0; j < 3; j++)
			squish += adjTriangle->areaVect[j] * centerToCenterVect[j];
		//squish /= sqrt(pow(centerToCenterVect[0], 2) + pow(centerToCenterVect[1], 2) + pow(centerToCenterVect[2], 2));
		if (squish < 0)	squish *= -1;
		squish = 1 - squish;
		if (squish > maxSquish)
			maxSquish = squish;
	}

	triangle->setSquish(maxSquish);
	return maxSquish;

}

void measureBySquish(Mesh* mesh) {

	// calculate the average of squish values
	double avgSquish = 0;
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		avgSquish += measureTriangleSquish(mesh->getTriangle(i), mesh);
	avgSquish = avgSquish / mesh->getNumOfTris();
	//cout << avgSquish << " ";

	// calculate the standard deviation of squish values
	double stdDeviation = 0;
	for (int i = 0; i < mesh->getNumOfTris(); i++)
		stdDeviation += pow(mesh->getTriangle(i)->getSquish() - avgSquish, 2);
	stdDeviation = sqrt(stdDeviation / mesh->getNumOfTris());
	cout << stdDeviation << " ";
}

double computeFrobeniusNorm(MatrixXd & M) {
	
	double norm = 0;
	for (int i = 0; i < M.rows(); i++)
		for (int j = 0; j < M.cols(); j++)
			norm += pow(M(i, j), 2);
	return sqrt(norm);

}

void measureByTransformation(Mesh* sourceMesh, Mesh* targetMesh, Mesh* inbetweenMesh) {

	sourceMesh->tetrahedralizeSurface();
	targetMesh->tetrahedralizeSurface();
	inbetweenMesh->tetrahedralizeSurface();

	Transformation* transformation = new Transformation(sourceMesh, targetMesh);
	
	double frobNorm, avgFrobNorm = 0;
	for (int i = 0; i < inbetweenMesh->getNumOfTris(); i++) {
	
		MatrixXd K, L, M;
		transformation->transformTetToTet(sourceMesh->getTetrahedron(i), inbetweenMesh->getTetrahedron(i), K);
		transformation->transformTetToTet(inbetweenMesh->getTetrahedron(i), targetMesh->getTetrahedron(i), L);
		transformation->transformTetToTet(sourceMesh->getTetrahedron(i), targetMesh->getTetrahedron(i), M); 
		MatrixXd Diff = K*L - M;
		frobNorm = computeFrobeniusNorm(Diff);
		inbetweenMesh->getTetrahedron(i)->setTransformationDiff(frobNorm);
		avgFrobNorm += frobNorm;
		
	}

	avgFrobNorm = avgFrobNorm / sourceMesh->getNumOfTris();
	//cout << avgFrobNorm << " ";

	double stdDeviation = 0;
	for (int i = 0; i < inbetweenMesh->getNumOfTris(); i++)
		stdDeviation += pow(inbetweenMesh->getTetrahedron(i)->getTransformationDiff() - avgFrobNorm, 2);
	stdDeviation = sqrt(stdDeviation / inbetweenMesh->getNumOfTris());
	cout << stdDeviation << " ";

}

void compareQuality() {

	// TODO: Delete the below data defining part and make generic

	/* ************ SCAPE DB ************ */

	// Method linear interpolation (lerp)
	string method_lerp_folder = "method_linear_interpolation/";
	string method_lerp_shapes[7] = { "mesh016.off", "lerp0_167.off", "lerp0_33.off", "lerp0_5.off",
	"lerp0_67.off", "lerp0_83.off", "mesh007.off" };

	// Method parts
	string method_parts_folder = "method_parts/";
	string method_parts_shapes[7] = { "mesh016.off", "parts0_167.off", "parts0_33.off", "parts0_5.off",
	"parts0_67.off", "parts0_83.off", "mesh007.off" };

	// Method data driven (not partitioned)
	string method_not_partitioned_folder = "method_data_driven/";
	string method_not_partitioned_shapes[7] = { "mesh016.off", "not_partitioned0_167.off", "not_partitioned0_33.off", "not_partitioned0_5.off",
	"not_partitioned0_67.off", "not_partitioned0_83.off", "mesh007.off" };

	// Method Alexa
	string method_alexa_folder = "method_alexa/";
	string method_alexa_shapes[7] = { "mesh016.off", "mesh016to007at0.167.off", "mesh016to007at0.333.off", "mesh016to007at0.5.off",
	"mesh016to007at0.667.off", "mesh016to007at0.833.off", "mesh007.off" };

	/* ************ HAND DB ************ */

	// Method Alexa
	string method_alexa_folder_hand = "DB_Hand/method_alexa/";
	string method_alexa_shapes_hand[7] = { "21.off", "21to17at0.167.off", "21to17at0.333.off", "21to17at0.5.off", "21to17at0.667.off",
	"21to17at0.833.off", "17.off" };

	string camel_pose = "camel-poses/camel-03.obj";

	/* **********************************  */

	Mesh* meshSequence_lerp[7], * meshSequence_parts[7], * meshSequence_not_partitioned[7], * meshSequence_alexa[7];

	// Load the meshes
	for (int i = 0; i < 7; i++) {
		meshSequence_lerp[i] = new Mesh();
		meshSequence_lerp[i]->loadOff((method_lerp_folder + method_lerp_shapes[i]).c_str());
		meshSequence_parts[i] = new Mesh();
		meshSequence_parts[i]->loadOff((method_parts_folder + method_parts_shapes[i]).c_str());
		meshSequence_not_partitioned[i] = new Mesh();
		meshSequence_not_partitioned[i]->loadOff((method_not_partitioned_folder + method_not_partitioned_shapes[i]).c_str());
		meshSequence_alexa[i] = new Mesh();
		meshSequence_alexa[i]->loadOff((method_alexa_folder + method_alexa_shapes[i]).c_str());
	}
/*
	// By Skewness
	cout << "/******************************************* Skewness *******************************************\\" << endl;
	cout << "       " << "Lerp" << "   " << "Parts" << "   " << "Not_partitioned" << "   " << "Alexa" << endl;
	for (int i = 0; i < 7; i++) {
		cout << "Mesh-" << i << " ";
		measureBySkewness(meshSequence_lerp[i]);
		measureBySkewness(meshSequence_parts[i]);
		measureBySkewness(meshSequence_not_partitioned[i]);
		measureBySkewness(meshSequence_alexa[i]);
		cout << endl;
	}

	// By Squish
	cout << "/******************************************* Squish *******************************************\\" << endl;
	cout << "       " << "Lerp" << "   " << "Parts" << "   " << "Not_partitioned" << "   " << "Alexa" << endl;
	for (int i = 0; i < 7; i++) {
		cout << "Mesh-" << i << " ";
		measureBySquish(meshSequence_lerp[i]);
		measureBySquish(meshSequence_parts[i]);
		measureBySquish(meshSequence_not_partitioned[i]);
		measureBySquish(meshSequence_alexa[i]);
		cout << endl;
	}

	// By Volume
	cout << "/******************************************* Volume *******************************************\\" << endl;
	cout << "       " << "Lerp" << "   " << "Parts" << "   " << "Not_partitioned" << "   " << "Alexa" << endl;
	for (int i = 0; i < 7; i++) {
		cout << "Mesh-" << i << " ";
		meshSequence_lerp[i]->computeVolume();
		cout << meshSequence_lerp[i]->getVolume() << " ";
		meshSequence_parts[i]->computeVolume();
		cout << meshSequence_parts[i]->getVolume() << " ";
		meshSequence_not_partitioned[i]->computeVolume();
		cout << meshSequence_not_partitioned[i]->getVolume() << " ";
		meshSequence_alexa[i]->computeVolume();
		cout << meshSequence_alexa[i]->getVolume() << " ";
		cout << endl;
	}
*/
	// By Transformation
	cout << "/******************************************* Transformation *******************************************\\" << endl;
	cout << "       " << "Lerp" << "   " << "Parts" << "   " << "Not_partitioned" << "   " << "Alexa" << endl;
	for (int i = 4; i < 7; i++) {
		cout << "Mesh-" << i << " ";
		measureByTransformation(meshSequence_lerp[0], meshSequence_lerp[6], meshSequence_lerp[i]);
		measureByTransformation(meshSequence_parts[0], meshSequence_parts[6], meshSequence_parts[i]);
		measureByTransformation(meshSequence_not_partitioned[0], meshSequence_not_partitioned[6], meshSequence_not_partitioned[i]);
		measureByTransformation(meshSequence_alexa[0], meshSequence_alexa[6], meshSequence_alexa[i]);
		cout << endl;
	}

}
