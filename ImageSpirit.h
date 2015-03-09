#ifndef IMAGESPIRIT_H
#define IMAGESPIRIT_H

#include <QtGui/QMainWindow>
#include "ui_ImageSpirit.h"
#include "ActionData.h"


// Corresponding of different objects and colors
class ObjColorList{
public:
	ObjColorList(QWidget *parent);
	QListWidget *_list; // Object list before and after interaction
	
	//void SetObjColorList(CStr &objName, vecI &objIdx);
	void setObjColorList(CMat &objLabel1u);

	void Initial(CStr &settingFileName);

	void show(bool isShow = true);

	static inline int Color2Int(Vec3b &color){return color[2] * 1000000 + color[1] * 1000 + color[0];}

	vecI _objIdx; //_objIdxB, 

	vector<QColor> _objColors;
	vecS _objNames;
	vecI _objColorIdx;
};

class ImageSpirit : public QMainWindow
{
	Q_OBJECT
	friend class PaintWidget;

public:
	ImageSpirit(QWidget *parent = 0, Qt::WFlags flags = 0);
	~ImageSpirit();

private:
	Ui::ImageSpiritClass ui;
	PaintWidget *paintWidget;
	QLabel _statusPnt, _statusMsg;
	ObjColorList _objColorList;
};


#endif // IMAGESPIRIT_H
