#pragma once

struct ActionType{
	int p; // Parameter of this action. E.g. representing 'down' or 'up' for move action
	int t; // Type of the action, should be one of the constant listed bellow 
	enum {NONE = -1, ACTIVATE, DEFORM, CHANGE_MAT, CHANGE_CLR, MOVE, REPEAT};

	ActionType(){p = t = NONE;}
};

class ObjDes{
public:
	int obj, clr, pos, mat; // Each represent the object, color, position, and material
	ObjDes() {obj = clr = pos = mat = -1;};

	void add(CStr &word);
	void print() const ;
	bool isEmpty() const {return obj == -1 && clr == -1 && pos == -1 && mat == -1;}
};


struct ActionData{
	friend class PaintWidget;
	friend class ObjDes;
	
public: // Global Data terms and load functions
	static vecS _objects, _objectsNoTrain; // Object names that are trained or not trained
	static vecS _materials; // Material attributes names
	static vecS _colors; // Color attributes names
	static vecS _positions; // Position attributes names
	static vecS _actions; // Name of different actions
	static CmColorIdx _colorIdx;

	static CStr MOVE_TYPE[], DEFORM_TYPE[];

	static bool loadGlobalSettings(CStr &ymlFileName);

	static QString showScript(CStr &script);

	static ActionType parseCmd(CStr &script, ObjDes &mainObj, ObjDes &obj2);

	static void printCmd(const ActionType &at, const ObjDes &obj, const ObjDes &obj2);

	static string getObjName(int objIdx);

	static void testScriptParsing();

	// Return -1 if not in the list
	static inline int findFromList(CStr &word, vecS &strList) {int idx = find(strList.begin(), strList.end(), word) - strList.begin(); return idx < strList.size() ? idx : -1;};
	//static inline int findFromList(CStr &word, CStr strList[], int num) {for (int i = 0; i < num; i++) if (word == strList[i]) return i; return -1;}

private:
	static int findFirstKeyword(const vecS &sentence, CStr keywords[], int numKeyword, int &keyIdx);
	
	static void getObjDescription(const vecS &sentence, int startIdx, int endIdx, ObjDes &objDes);
};


class Action{
	friend class PaintWidget;
public:

	void Draw();

	void Next(PaintWidget &pw);

	const QString& Script() {return _script;}
	const char* script() {return _S(_scriptStr); }
	bool IsAnimate() { return _maxT > 0;}

	static int LoadActions(PaintWidget &pw, vector<Action> &acts);

	BOOL _showBeforeAction, _skipWhenLater;


private:
	GLuint _objTexture; // GL texture
	Mat _coeff; // Coefficients for the action
	vector<Point2d> _accuShifts; // Accumulated shifts


	string _fName; // File name of the texture
	double _w, _h; // Size of the texture
	QString _script; // Action script
	string _scriptStr; // 
	enum ActionType{NONE, DEFORM, DISPLAY} _actType; 
	int _t, _maxT, _timerInterval;

	// For displaying video
	VideoCapture _vc; // video capture
	Mat _alphaMat; // Alpha mat for displaying video

	// For complicated polygons
	vector<vecI> _polygons;
};

