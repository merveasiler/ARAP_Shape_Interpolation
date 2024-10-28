#include "pch.h"
#include "Painter.h"

void Painter::getShapeSep(Mesh* mesh, SoSeparator* res)
{
	// Paint all vertices with the same color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0.5, 0.5, 0.5);
	res->addChild(mat);

	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++)
		coords->point.set1Value(c, mesh->getVertex(c)->coords[0],
			mesh->getVertex(c)->coords[1],
			mesh->getVertex(c)->coords[2]);

	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->getNumOfTris(); c++)
	{
		faceSet->coordIndex.set1Value(c * 4, mesh->getTriangle(c)->corners[0]);
		faceSet->coordIndex.set1Value(c * 4 + 1, mesh->getTriangle(c)->corners[1]);
		faceSet->coordIndex.set1Value(c * 4 + 2, mesh->getTriangle(c)->corners[2]);
		faceSet->coordIndex.set1Value(c * 4 + 3, -1);
	}

	res->addChild(coords);
	res->addChild(faceSet);

}

void Painter::drawTriangulation(Mesh* mesh, SoSeparator* res) {

	SoSeparator* thickEdgeSep = new SoSeparator;
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	thickEdgeSep->addChild(sty);

	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	for (int se = 0; se < mesh->getNumOfEdges(); se++)
	{
		SbVec3f end1 = mesh->getVertex(mesh->getEdge(se)->endVerts[0])->coords;
		SbVec3f end2 = mesh->getVertex(mesh->getEdge(se)->endVerts[1])->coords;
		co->point.set1Value(2 * se, end1);
		co->point.set1Value(2 * se + 1, end2);
	}

	for (int ci = 0; ci < mesh->getNumOfEdges(); ci++)
	{
		ils->coordIndex.set1Value(3 * ci, 2 * ci);
		ils->coordIndex.set1Value(3 * ci + 1, 2 * ci + 1);
		ils->coordIndex.set1Value(3 * ci + 2, -1);
	}

	thickEdgeSep->addChild(co);
	thickEdgeSep->addChild(ils);
	res->addChild(thickEdgeSep);

}

void Painter::drawSingleTriangle(Mesh* mesh, Triangle* triangle, SoSeparator* res) {

	SoSeparator* thickEdgeSep = new SoSeparator;
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);	// blue
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	thickEdgeSep->addChild(sty);

	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	for (int se = 0; se < 3; se++) {
		Edge* edge = mesh->getEdge(triangle->edgeList[se]);
		SbVec3f end1 = mesh->getVertex(edge->endVerts[0])->coords;
		SbVec3f end2 = mesh->getVertex(edge->endVerts[1])->coords;
		co->point.set1Value(2 * se, end1);
		co->point.set1Value(2 * se + 1, end2);

		ils->coordIndex.set1Value(3 * se, 2 * se);
		ils->coordIndex.set1Value(3 * se + 1, 2 * se + 1);
		ils->coordIndex.set1Value(3 * se + 2, -1);
	}

	thickEdgeSep->addChild(co);
	thickEdgeSep->addChild(ils);
	res->addChild(thickEdgeSep);

}