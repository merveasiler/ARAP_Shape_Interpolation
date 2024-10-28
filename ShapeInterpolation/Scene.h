#pragma once

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/SoPickedPoint.h>

SoPath* pickFilterCB(void *, const SoPickedPoint* pick)
{
	// See which child of selection got picked
	SoPath* p = pick->getPath();
	int i;
	for (i = 0; i < p->getLength() - 1; i++) {
		SoNode *n = p->getNode(i);
		if (n->isOfType(SoSelection::getClassTypeId()))
			break;
	}

	// Copy 2 nodes from the path:
	// selection and the picked child
	return p->copy(i, 2);
}

class Scene {

	HWND window;
	SoWinExaminerViewer * viewer;
	SoSeparator * root;
	SoSelection * selection;

public:

	Scene() {

		window = SoWin::init("Shape Interpolation");
		viewer = new SoWinExaminerViewer(window);
		root = new SoSeparator;
		root->ref();

	}

	void attachToRoot(SoSeparator* res) {
		root->addChild(res);
	}

	void configureSelection(SoSeparator* res) {

		selection = new SoSelection;
		selection->addChild(res);
		selection->setPickFilterCallback(pickFilterCB);

	}

	void configureViewer() {

		viewer->setBackgroundColor(SbColor(1, 1, 1));
		viewer->setSize(SbVec2s(640, 480));
		viewer->setSceneGraph(root);
		//viewer->setSceneGraph(selection);
	}

	void play() {

		viewer->show();
		SoWin::show(window);
		SoWin::mainLoop();

	}

	~Scene() {
		root->unref();
		//selection->removeAllChildren();
		delete viewer;
	}

	void makeScene(SoSeparator* res) {
		
		attachToRoot(res);
		//configureSelection(res);
		configureViewer();
		play();

	}	

};

#pragma once
