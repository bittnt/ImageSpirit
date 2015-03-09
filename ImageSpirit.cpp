#include "stdafx.h"
#include "ImageSpirit.h"
#include "PaintWidget.h"

ImageSpirit::ImageSpirit(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
	, _objColorList(parent)
{
	ui.setupUi(this);
	ui.mainToolBar->addWidget(ui.comboBoxAttributes);
	_statusPnt.setMinimumWidth(120);
	statusBar()->addWidget(&_statusPnt);
	statusBar()->addWidget(&_statusMsg);

	// Main layout
	QHBoxLayout *mainLayout = new QHBoxLayout();// QHBoxLayout();
	QScrollArea *paintScroll = new QScrollArea;
	paintScroll->setMinimumHeight(400);
	paintWidget = new PaintWidget(this);
	paintScroll->setWidget(paintWidget);
	mainLayout->addWidget(paintScroll);
	ui.centralWidget->setLayout(mainLayout);
	mainLayout->addWidget(_objColorList._list);
}

ImageSpirit::~ImageSpirit()
{
}


ObjColorList::ObjColorList(QWidget *parent)
{
	_list =	new QListWidget(parent);
	_list->setObjectName(QString::fromUtf8("Object color list widget"));
	_list->setMaximumWidth(100);
	_list->setSortingEnabled(false);
	_list->hide();
}


//void ObjColorList::SetObjColorList(CStr &objName, vecI &objIdx)
//{
//	Mat objImg = imread(objName);
//	if (objImg.data == NULL){
//		printf("Can't load file: %s\n", _S(objName));
//		return;
//	}
//	map<int, int> codeCount;
//	for (int r = 0; r < objImg.rows; r++){
//		Vec3b *val = objImg.ptr<Vec3b>(r);
//		for (int c = 0; c < objImg.cols; c++){
//			int code = Color2Int(val[c]); 
//			codeCount[code] ++;
//		}
//	}
//	objIdx.clear();
//	for (size_t i = 0; i < _objColorIdx.size(); i++)
//		if (codeCount[_objColorIdx[i]] > 100)	
//			objIdx.push_back(i);
//}

void ObjColorList::setObjColorList(CMat &objLabel1u)
{
	map<int, int> codeCount;
	for (int r = 0; r < objLabel1u.rows; r++){
		const byte *val = objLabel1u.ptr<byte>(r);
		for (int c = 0; c < objLabel1u.cols; c++)
			codeCount[val[c]] ++;
	}
	_objIdx.clear();
	for (auto it = codeCount.begin(); it != codeCount.end(); it++)
		if (it->second > 100)
			_objIdx.push_back(it->first);
}

void ObjColorList::Initial(CStr &settingFileName)
{
	int objNum = (int)ActionData::_objects.size(), objNumNT = (int)ActionData::_objectsNoTrain.size();
	_objColors.resize(objNum + objNumNT);
	_objNames.resize(objNum + objNumNT);
	_objColorIdx.resize(objNum + objNumNT);
	for (int i = 0; i < objNum + objNumNT; i++){
		Vec3b color = ActionData::_colorIdx[i];
		_objColors[i] = QColor(color[2], color[1], color[0]);
		_objNames[i] = i < objNum ? ActionData::_objects[i] : ActionData::_objectsNoTrain[i - objNum];
		_objColorIdx[i] = Color2Int(color);
	}
}

void ObjColorList::show(bool isShow)
{
	const vecI &nodes = _objIdx; 
	if (isShow){
		_list->clear();
		for (size_t i = 0; i < nodes.size(); i++){
			int idx = nodes[i];
			QColor clr = _objColors[idx];
			int bkcolr = clr.red() + clr.green() + clr.blue() < 128*3 ? 255 : 0;
			QBrush brush(QColor(bkcolr, bkcolr, bkcolr));
			brush.setStyle(Qt::SolidPattern);
			QListWidgetItem *qlistwidgetitem = new QListWidgetItem(_list);
			qlistwidgetitem->setForeground(brush);
			brush.setColor(clr);
			qlistwidgetitem->setBackground(brush);
			_list->item(i)->setText(QString(_S(_objNames[idx])));
		}
		_list->show();
	}
	else
		_list->hide();
}

