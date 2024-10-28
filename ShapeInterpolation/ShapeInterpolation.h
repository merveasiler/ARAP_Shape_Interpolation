#pragma once

#include "Scene.h"
#include "Mesh.h"
#include "Metric.h"
#include "Painter.h"

#include <Inventor/nodes/SoEventCallBack.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/events/SoMouseButtonEvent.h>

vector< Mesh* > meshSequence;

vector<string> split(string inputStr, string delimiter) {

	vector<string> pieces;
	string data = inputStr;

	while (data != "") {
		int ind = data.find_first_of(delimiter);
		pieces.push_back(data.substr(0, ind));
		data = data.substr(ind + 1);
	}

	return pieces;

}

void drawMeshToScene(string meshName) {

	Mesh* mesh = new Mesh;
	mesh->loadOff(meshName.c_str());

	Scene* scene = new Scene();
	Painter* painter = new Painter();
	SoSeparator* res = new SoSeparator();
	painter->getShapeSep(mesh, res);
	scene->makeScene(res);

	delete scene;
	delete painter;
	delete mesh;

}

void interpolate(string sourceMeshName, string targetMeshName, float lambda, string resultMeshName) {

	Mesh* sourceMesh = new Mesh();
	sourceMesh->loadOff(sourceMeshName.c_str());

	Mesh* targetMesh = new Mesh();
	targetMesh->loadOff(targetMeshName.c_str());

	Transformation* transformation = new Transformation(sourceMesh, targetMesh);
	transformation->interpolateARAPByAlexa(lambda, resultMeshName);

	delete sourceMesh;
	delete targetMesh;
	delete transformation;

	drawMeshToScene(resultMeshName);

}

/*
void setupEvents() {

	SoSeparator* selectionRoot;
	SoEventCallback* eventCB = new SoEventCallback;
	eventCB->addEventCallback(SoMouseButtonEvent::getClassTypeId(), identifyEventCB, selectionRoot);
	selectionRoot->addChild(eventCB);

}

void identifyEventCB(void* userData, SoEventCallback* eventCB)
{
	SoSelection *selection = (SoSelection *)userData;
	const SoEvent *event = eventCB->getEvent();
	SoBut

	if (SO_MOUSE_PRESS_EVENT(event, LEFT)) {

	}

	// Check for the Up and Down arrow keys being pressed.
	if (SO_KEY_PRESS_EVENT(event, UP_ARROW)) {
		myScaleSelection(selection, 1.1);
		eventCB->setHandled();
	}
	else if (SO_KEY_PRESS_EVENT(event, DOWN_ARROW)) {
		myScaleSelection(selection, 1.0 / 1.1);
		eventCB->setHandled();
	}
}
*/
