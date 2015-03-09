#ifndef PAINTWIDGET_H
#define PAINTWIDGET_H

#include <QObject>
#include "ActionData.h"
#include "ImageSpirit.h"

class ImageSpirit;

class PaintWidget : public QGLWidget
{
	Q_OBJECT
public:
	PaintWidget(ImageSpirit *parent);
	~PaintWidget();
	friend class Action;
	friend class PaintWidget;

private: // GUI related parameters and files to be loaded 
	ImageSpirit &_uiMain;
	Ui::ImageSpiritClass &_ui;
	QLabel &_statusMsg, &_statusPnt;
	ObjColorList &_objColorList;
	CmMatShowGL _bgShow3u;
	Mat _staticShow3u;
	string _dir, _nameNE, _dirData;
	FileStorage _fs;
	int _crntViewAtri, _viewType, _crntTool;
	enum ViewType{VIEW_SRC, VIEW_OBJ, VIEW_SOFT_TGT, VIEW_ATTRI, VIEW_RES, VIEW_MARK};
	enum ToolType{TOOL_LINE, TOOL_BRUSH, TOOL_FILL};
	Point _downPnt, _prePnt;

private: // Files to be loaded
	vector<Mat> _attrMasks; // current atribute exist (!= 0) or not (== 0) 
	Mat _objLabel1u, _img3u, _bg3u, _depth1u, _img3f;
	Mat _inpLabel1u; // object label for background in-painting.
	int _w, _h, _bSize;

private: // Variables related with manipulation actions
	QTimer _timerAnimate; 
	vector<Action> _acts;
	int _crntAct;

	ActionType _actType; // Current action type
	ObjDes _objDesM, _objDesA; // The main and annex object descriptor for current act
	Mat _tgtAlpha1u, _annexMask; // Alpha mat of the main target and the mask of the annex object

private:
	bool loadData(CStr &fPath);
	void loadAttriTexture();
	void loadIniObjLabel(CStr &pathNE);

	enum{SHOW_NONE = 0, SHOW_BRUSH, SHOW_LINE};
	void viewStatic(Point pnt = Point(), int t = SHOW_NONE|SHOW_BRUSH);
	void paint2D();
	void paint3D();

	void setToolType(ToolType t = TOOL_LINE);
	void setViewType(ViewType t = VIEW_SRC);

	void inpaintingBg(bool checkChanged = false);
	bool checkMouseButton(Qt::MouseButtons &buttons);

	void analysisTargets(); // Get target object alpha

private: // GUI related functions
	void leaveEvent(QEvent *evt);
	void mouseMoveEvent(QMouseEvent *evt);
	void mousePressEvent(QMouseEvent *evt);
	void mouseReleaseEvent(QMouseEvent *evt);
	void initializeGL();
	void paintGL();

	void resizeGL(int w, int h);
	void connectSignals();
	bool popWarning(const QString &wType, const QString &str);

	void setupViewport(int width, int height);

private slots:
	void fileOpen();
	void fileSave();

	void toolEdit();
	void toolFinished();
	void toolInpaint();
	void toolExportOldData();
	
	void toolLine();
	void toolBrush();
	void toolFill();

	void viewImage();
	void viewSoftTarget();
	void viewObjects();
	void viewResults();
	void viewAttributes(int idx);
	void viewMark();

	void animateObjs();
};

#endif // PAINTWIDGET_H
