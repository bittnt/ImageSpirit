#include "StdAfx.h"
#include "ActionData.h"
#include "PaintWidget.h"

vecS ActionData::_objects;
vecS ActionData::_objectsNoTrain;
vecS ActionData::_materials;
vecS ActionData::_colors;
vecS ActionData::_positions;
vecS ActionData::_actions;
CmColorIdx ActionData::_colorIdx;

CStr ActionData::MOVE_TYPE[] = {"down", "up", "left", "right"};
CStr ActionData::DEFORM_TYPE[] = {"lower", "taller", "smaller", "larger"};


void ObjDes::add(CStr &word)
{
	int idx = ActionData::findFromList(word, ActionData::_objects);
	if (idx != -1){
		obj = idx;
		return;
	}
	idx = ActionData::findFromList(word, ActionData::_materials);
	if (idx != -1){
		mat = idx;
		return;
	}
	idx = ActionData::findFromList(word, ActionData::_positions);
	if (idx != -1){
		pos = idx;
		return;
	}
	idx = ActionData::findFromList(word, ActionData::_colors);
	if (idx != -1)
		clr = idx;
}

void ObjDes::print() const
{
	if (clr != -1)
		printf("%s ", _S(ActionData::_colors[clr]));
	if (mat != -1)
		printf("%s ", _S(ActionData::_materials[mat]));
	if (obj != -1)
		printf("%s ", _S(ActionData::_objects[obj]));
	if (pos != -1)
		printf("in %s ", _S(ActionData::_positions[pos]));
}

bool ActionData::loadGlobalSettings(CStr &ymlFileName)
{
	FileStorage fs(ymlFileName, FileStorage::READ);
	if (!fs.isOpened())
		return false;
	fs["Objects"] >> _objects;
	fs["ObjectsNoTrain"] >> _objectsNoTrain;
	fs["Materials"] >> _materials;
	fs["Colors"] >> _colors;
	fs["Positions"] >> _positions;
	fs["Actions"] >> _actions;
	_colorIdx = CmColorIdx(_objects.size() + _objectsNoTrain.size());
	printf("Setting summary: #bj = %d, #ObjNoTrain = %d, #Materials = %d, #Colors = %d, #Pos = %d, Actions = %d\n",
		_objects.size(), _objectsNoTrain.size(), _materials.size(), _colors.size(), _positions.size(), _actions.size());

	return true;
}

QString ActionData::showScript(CStr &script)
{
	vector<string> words;
	splitStr(script, " ", words); //CmLog::split(script, " ", words);
	string resStr = "<p style=\"font-weight:600; margin-right:5px;\">";
	const char* color = NULL;
	for (size_t i = 0; i < words.size(); i++){
		if (findFromList(words[i], _actions) != -1)
			color = "ff0000";
		else if (findFromList(words[i], _objects) != -1 || findFromList(words[i], _objectsNoTrain) != -1)
			color = "0000ff";
		else if (findFromList(words[i], _materials) != -1 || findFromList(words[i], _colors) != -1 || findFromList(words[i], _positions) != -1)
			color = "00acef";
		else
			color = "000000";
		resStr += format(" <span style=\"color:#%s\">%s</span>", color, _S(words[i]));
	}
	resStr += ". </p>";
	return QString(_S(resStr));
}

void ActionData::getObjDescription(const vecS &sentence, int startIdx, int endIdx, ObjDes &objDes)
{
	for (int i = startIdx; i < endIdx; i++)
		objDes.add(sentence[i]);
}


int ActionData::findFirstKeyword(const vecS &sentence, CStr keywords[], int numKeyword, int &keyIdx){
	int num = (int)sentence.size();
	int res = num;
	for (int i = 2; i < num; i++){
		for (int j = 0; j < numKeyword; j++){
			if (sentence[i] == keywords[j]){
				res = i;
				keyIdx = j;
				break;
			}
		}
		if (res != num)
			break;
	}
	return res;
}


ActionType ActionData::parseCmd(CStr &script, ObjDes &mainObj, ObjDes &obj2)
{
	// Find action type
	printf("`%s':\n", _S(script));
	ActionType at;
	vecS words;
	int endIdx;{
		splitStr(script, " ", words); 
		endIdx = words.size();
		if (words[0] == "Activate")
			at.t = ActionType::ACTIVATE;
		else if (words[0] == "Make"){
			at.t = ActionType::DEFORM;
			endIdx = findFirstKeyword(words, DEFORM_TYPE, sizeof(DEFORM_TYPE)/sizeof(string), at.p);
		}
		else if (words[0] == "Change"){
			int mIdx = findFromList(words[words.size() - 1], _materials);
			if (mIdx != -1)	{ 
				at.t = ActionType::CHANGE_MAT;
				at.p = mIdx;
			}
			else{
				int cIdx = findFromList(words[words.size() - 1], _colors);
				if (cIdx != -1)	{
					at.t = ActionType::CHANGE_CLR;
					at.p = cIdx;
				}
			}
			if (at.t == -1)
				return at;
			mIdx = findFromList("from", words);
			if (mIdx != -1)
				mainObj.add(words[mIdx + 1]);
		}
		else if (words[0] == "Move"){
			at.t = ActionType::MOVE;
			endIdx = findFirstKeyword(words, MOVE_TYPE, sizeof(MOVE_TYPE)/sizeof(string), at.p);
		}
		else if (words[0] == "Repeat"){
			at.t = ActionType::REPEAT;
			findFirstKeyword(words, MOVE_TYPE, sizeof(MOVE_TYPE)/sizeof(string), at.p);
		}
		else {
			printf("Unrecognized action\n");
			return at;
		}
	}

	// Parse the object descriptions
	{
		if (endIdx == words.size())	{
			CStr stopWords[] = {"from", "to", "and"};
			endIdx = findFirstKeyword(words, stopWords, sizeof(stopWords)/sizeof(string), dummyI);
		}
		getObjDescription(words, 2, endIdx, mainObj);
		int alongPos = findFromList("along", words);
		if (alongPos != -1)
			getObjDescription(words, alongPos + 2, words.size(), obj2);
	}
	return at;
}

void ActionData::printCmd(const ActionType &at, const ObjDes &obj, const ObjDes &obj2)
{
	CStr actionWords[] = {"Activate", "Make", "Change", "Change", "Move", "Repeat"};
	printf("%s the ", _S(actionWords[at.t]));
	obj.print();
	switch (at.t){
	case ActionType::DEFORM: printf("%s ", _S(DEFORM_TYPE[at.p])); break;
	case ActionType::CHANGE_CLR: printf("to %s ", _S(_colors[at.p])); break;
	case ActionType::CHANGE_MAT: printf("to %s ", _S(_materials[at.p])); break;
	case ActionType::MOVE: printf("%s ", _S(MOVE_TYPE[at.p])); break;
	case ActionType::REPEAT: printf("and move %s ", _S(MOVE_TYPE[at.p])); break;
	}
	if (!obj2.isEmpty()){
		printf("along the ");
		obj2.print();
	}
	printf("\n\n");
}

string ActionData::getObjName(int objIdx)
{
	if (objIdx < _objects.size())
		return _objects[objIdx];
	
	objIdx -= _objects.size();
	if (objIdx < _objectsNoTrain.size())
		return _objectsNoTrain[objIdx];
	return string();
}


int Action::LoadActions(PaintWidget &pw, vector<Action> &acts){
	for (size_t i = 0; i < acts.size(); i++)
		pw.deleteTexture(acts[i]._objTexture);
	acts.clear();
	acts.reserve(20);

	string pathNE = pw._dirData;
	FileNode n = pw._fs["Actions"];
	FileNodeIterator it = n.begin(), end = n.end();
	for (int idx = 0; it != end; ++it, idx++){
		acts.push_back(Action());
		Action &act = acts[idx];
		{// Action type
			string type;
			(*it)["type"] >> type;
			if (type == "Deform") 
				act._actType = DEFORM;
			else if (type == "Show")
				act._actType = DISPLAY;
			else {
				act._actType = NONE;
			}
		}
		{ // General information
			(*it)["name"] >> act._fName;
			act._fName = pathNE + act._fName;
			(*it)["maxT"] >> act._maxT;
			//(*it)["waitTime"] >> act._sleepTime;
			(*it)["timerInterval"] >> act._timerInterval;
			(*it)["showBefore"] >> act._showBeforeAction;	
			(*it)["skipWhenLater"] >> act._skipWhenLater;
			act._t = 0;
			(*it)["script"] >> act._scriptStr;
			act._script = ActionData::showScript(act._scriptStr);
		}
		{ // Action coefficients
			string coeff;
			(*it)["coeff"] >> coeff;
			pw._fs[coeff] >> act._coeff;	
			CV_Assert(act._coeff.type() == CV_64FC4);
			act._accuShifts.resize(act._coeff.cols);
			memset(&act._accuShifts[0], 0, sizeof(Point2d) * act._coeff.cols);
		}
		{ // Load textures
			Mat texImg;
			if (CmFile::GetExtention(act._fName) == ".png")
				texImg = imread(act._fName, CV_LOAD_IMAGE_UNCHANGED);
			else{
				CV_Assert(act._actType = DISPLAY);
				act._vc.open(act._fName);
				CV_Assert_(act._vc.isOpened(), ("Video not opened: %s", _S(act._fName)));
				Mat videoImg;

				string fName = CmFile::GetPathNE(act._fName) + "/001.jpg";
				if (CmFile::FileExist(fName)){
					videoImg = imread(fName);
					videoImg = videoImg(Rect(380, 0, 640, 376));
				}
				else
					act._vc >> videoImg;

				CV_Assert(videoImg.data != NULL && videoImg.type() == CV_8UC3);
				act._alphaMat = Mat::zeros(videoImg.size(), CV_8UC1);
				act._alphaMat(Rect(3, 3, videoImg.cols - 6, videoImg.rows - 6)).setTo(255);
				blur(act._alphaMat, act._alphaMat, Size(7, 7));
				texImg = CmCv::Merge(videoImg, act._alphaMat);
			}
			CV_Assert_(texImg.type() == CV_8UC4, ("Failed to load %s\n", _S(act._fName)));
			act._w = texImg.cols;
			act._h = texImg.rows;
			QImage qt_img(texImg.data, texImg.cols, texImg.rows, QImage::Format_ARGB32); // = ( QImage ( cv_img->dataIm, cv_img->width, cv_img->height, cv_img->QImage::Format_RGB888 ) ).rgbSwapped();
			act._objTexture = pw.bindTexture(QPixmap::fromImage(qt_img), GL_TEXTURE_2D);//QPixmap pm(_QS(act._fName));		
		}

		{// Load Polygons
			string polyName;
			(*it)["polyName"] >> polyName;
			if (polyName.size() == 0){
				act._polygons.resize(1);
				vecI &ploy = act._polygons[0];
				ploy.resize(act._coeff.cols);
				for (size_t p = 0; p < ploy.size(); p++)
					ploy[p] = p;
			}
			else{
				FileNode nPoly = pw._fs[polyName];
				FileNodeIterator itPoly = nPoly.begin(), itPolyEnd = nPoly.end();
				for (; itPoly != itPolyEnd; itPoly++)	{
					vecI crntPoly;
					(*itPoly)["poly"] >> crntPoly;
					act._polygons.push_back(crntPoly);
				}
			}

		}
	}
	return acts.size();
}

void Action::Draw()
{
	const Point2d corners[4] = {Point2d(0, 1), Point2d(0, 0), Point2d(1, 0), Point2d(1, 1)};
	if (_actType == DISPLAY)
		CV_Assert(_polygons.size() == 1 && _polygons[0].size() == 4);

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, _objTexture);
	for (size_t p = 0; p < _polygons.size(); p++) {
		vecI &ploy = _polygons[p];
		glBegin(GL_POLYGON);
		for (size_t i = 0; i < ploy.size(); i++){
			Vec4d &cof = _coeff.at<Vec4d>(0, ploy[i]);
			if (_actType == DISPLAY)
				glTexCoord2d(corners[i].x, corners[i].y);
			else
				glTexCoord2d(cof[0]/_w, 1-cof[1]/_h);
			//glVertex2d(cof[0] + _t * cof[2], cof[1] + _t * cof[3]); 
			Point2d shiftVal = _accuShifts[ploy[i]];
			glVertex2d(cof[0] + shiftVal.x, cof[1] + shiftVal.y); 
		}
		glEnd();
		glClear(GL_DEPTH_BUFFER_BIT);
	}

	/*
	glDisable(GL_TEXTURE_2D);
	for (size_t i = 0; i < pnts.size(); i++){
		Vec3f clr = CmShow::gColors[i];
		clr /= 255;
		glColor3f(clr[0], clr[1], clr[2]) ;
		glBegin(GL_LINES) ; 
		glVertex2d(pnts[i].x, pnts[i].y);
		int idx = (i + 1) % pnts.size();
		glVertex2d(pnts[idx].x, pnts[idx].y);
		glEnd() ;
	}
	glEnable(GL_TEXTURE_2D);
	glClear(GL_DEPTH_BUFFER_BIT); //*/
}

void Action::Next(PaintWidget &pw)
{
	if (_actType == NONE)
		return;
	_t++;

	if (_t > _maxT){
		_t= _maxT;
		if (_maxT > 0)
			pw.toolFinished();
		return;
	}

	for (int i = 0; i < _coeff.cols; i++){
		Vec4d val = _coeff.at<Vec4d>(0, i);
		_accuShifts[i] += Point2d(val[2], val[3]);
	}

	if (_actType == DISPLAY){
		Mat videoImg, texImg;
		string fName = CmFile::GetPathNE(_fName) + format("/%03d.jpg", _t);
		if (CmFile::FileExist(fName)){
			videoImg = imread(fName);
			videoImg = videoImg(Rect(380, 0, 640, 376));
		}
		else
			_vc >> videoImg;

		if (videoImg.data != NULL){
			texImg = CmCv::Merge(videoImg, _alphaMat);
			QImage qt_img(texImg.data, texImg.cols, texImg.rows, QImage::Format_ARGB32); // = ( QImage ( cv_img->dataIm, cv_img->width, cv_img->height, cv_img->QImage::Format_RGB888 ) ).rgbSwapped();
			pw.deleteTexture(_objTexture);
			_objTexture = pw.bindTexture(QPixmap::fromImage(qt_img), GL_TEXTURE_2D);//QPixmap pm(_QS(act._fName));
		}
	}
}

void ActionData::testScriptParsing()
{
	CV_Assert(loadGlobalSettings("Setting.yml"));

	const char* scripts[] = {"Activate the black car in bottom-right", 
		"Make the picture in top-middle larger",
		"Make the cabinet lower",
		"Make the wood cabinet in bottom-middle lower",
		"Make the picture in top-left larger",
		"Change the blinds in center-right from wood to cotton",
		"Change the black chair in bottom-right from leather to wood",
		"Move the picture in top-right right",
		"Repeat the chair in bottom-left and move right along the white floor",
	};
	const int NUM = sizeof(scripts)/sizeof(char*);
	for (int i = 0; i < NUM; i++){
		ObjDes objDes, objDes2;
		ActionType at = parseCmd(scripts[i], objDes, objDes2);
		printCmd(at, objDes, objDes2);
	}

}