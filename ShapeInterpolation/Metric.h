#pragma once

#include "Mesh.h"
#include "Transformation.h"
#define _PI 3.14159265358979323846

template <class Type>
Type findMax(vector<Type> elementList);

template <class Type>
Type findMin(vector<Type> elementList);

double measureTriangleSkewness(Triangle* triangle);

void measureBySkewness(Mesh* mesh);

double measureTriangleSquish(Triangle* triangle, Mesh* mesh);

void measureBySquish(Mesh* mesh);

void measureByTransformation(Mesh* sourceMesh, Mesh* targetMesh, Mesh* inbetweenMesh);

void compareQuality();

