#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoIndexedLineSet.h>

#include "Mesh.h"


class Painter
{
public:
	void getShapeSep(Mesh* mesh, SoSeparator* res);
	void drawTriangulation(Mesh* mesh, SoSeparator* res);
	void drawSingleTriangle(Mesh* mesh, Triangle* triangle, SoSeparator* res);
};

#pragma once
