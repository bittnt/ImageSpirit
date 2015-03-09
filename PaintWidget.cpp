#include "StdAfx.h"
#include "PaintWidget.h"
#include "ImageSpirit.h"


PaintWidget::PaintWidget(ImageSpirit *parent)
	: QGLWidget(parent)
	, _crntAct(-1)
	, _uiMain(*parent)
	, _ui(parent->ui)
	, _statusMsg(parent->_statusMsg)
	, _statusPnt(parent->_statusPnt)
	, _objColorList(parent->_objColorList)
	, _crntViewAtri(0)
	, _bSize(3) // brush size
{
	setMinimumSize(620, 460);
	setBackgroundRole(QPalette::Base);
	setAttribute(Qt::WA_PaintOutsidePaintEvent);
	setContextMenuPolicy(Qt::ActionsContextMenu);
	setMouseTracking(true);
	setFocusPolicy(Qt::StrongFocus);
	connectSignals();
	setToolType(TOOL_LINE);
	setViewType(VIEW_SRC);

	if (!CmFile::FileExist("Setting.yml"))
		popWarning("Can't load file", "Setting.yml");
	else{
		ActionData::loadGlobalSettings("Setting.yml");
		_objColorList.Initial("Setting.yml");
	}
}

PaintWidget::~PaintWidget()
{
}

void PaintWidget::connectSignals()
{
	connect(_ui.actionOpen, SIGNAL(triggered()), this, SLOT(fileOpen()));
	connect(_ui.actionSave, SIGNAL(triggered()), this, SLOT(fileSave()));
	
	connect(_ui.actionImage, SIGNAL(triggered()), this, SLOT(viewImage()));
	connect(_ui.actionSoftTarget, SIGNAL(triggered()), this, SLOT(viewSoftTarget()));
	connect(_ui.actionObjects, SIGNAL(triggered()), this, SLOT(viewObjects()));
	connect(_ui.actionResults, SIGNAL(triggered()), this, SLOT(viewResults()));
	connect(_ui.actionUserMark, SIGNAL(triggered()), this, SLOT(viewMark()));
	connect(_ui.comboBoxAttributes, SIGNAL(currentIndexChanged(int)), this, SLOT(viewAttributes(int)));
	
	connect(_ui.actionEdit, SIGNAL(triggered()), this, SLOT(toolEdit()));
	connect(_ui.actionFinished, SIGNAL(triggered()), this, SLOT(toolFinished()));
	connect(_ui.actionInPainting, SIGNAL(triggered()), this, SLOT(toolInpaint()));
	connect(_ui.actionExportOldData, SIGNAL(triggered()), this, SLOT(toolExportOldData()));

	connect(_ui.actionLine, SIGNAL(triggered()), this, SLOT(toolLine()));
	connect(_ui.actionBrush, SIGNAL(triggered()), this, SLOT(toolBrush()));
	connect(_ui.actionFill, SIGNAL(triggered()), this, SLOT(toolFill()));

	connect(&_timerAnimate, SIGNAL(timeout()), this, SLOT(animateObjs()));
}

void PaintWidget::leaveEvent(QEvent *evt)
{
	_uiMain._statusPnt.setText(tr(""));
}

bool PaintWidget::checkMouseButton(Qt::MouseButtons &buttons)
{
	if (buttons & Qt::LeftButton)
		setToolType(TOOL_LINE);
	else if (buttons & Qt::MidButton)
		setToolType(TOOL_BRUSH);
	else if (buttons & Qt::RightButton)
		setToolType(TOOL_FILL);
	else 
		return false;
	return true;
}

void PaintWidget::analysisTargets()
{
	compare(_objLabel1u, _objDesM.obj, _tgtAlpha1u, CMP_EQ);
	imshow("Target object region", _tgtAlpha1u);

	if (_objDesA.obj != -1){
		compare(_objLabel1u, _objDesA.obj, _annexMask, CMP_EQ);
		imshow("Annex object mask", _annexMask);
		waitKey(1);		
	}
}

void PaintWidget::mouseMoveEvent(QMouseEvent *evt)
{
	if (_objLabel1u.data == NULL)
		return;

	// Display information
	Point newPnt(evt->x(), evt->y());
	string objName = ActionData::getObjName(_objLabel1u.at<byte>(newPnt));
	string msg = cv::format("(%03d, %03d): %s ", newPnt.x, newPnt.y, _S(objName));
	_uiMain._statusPnt.setText(_QS(msg)); 

	// Draw user mark for object class or in-painting class illustration
	if (_prePnt == newPnt || (_viewType != VIEW_MARK && _viewType != VIEW_OBJ))
		return;

	// Setting for the mouse operations, return if no button pressed
	if (!checkMouseButton(evt->buttons()))
		return;

	if (_crntTool == TOOL_LINE || _crntTool == TOOL_FILL)
		viewStatic(newPnt, SHOW_LINE);
	else {
		Mat objLabel1u = _viewType == VIEW_OBJ ? _objLabel1u : _inpLabel1u;
		line(objLabel1u, _prePnt, newPnt, _objLabel1u.at<byte>(_downPnt), _bSize, _bSize == 1 ? 8 : -1);
		viewStatic(newPnt, SHOW_BRUSH);
	}

	_prePnt = newPnt;
}

void PaintWidget::mousePressEvent(QMouseEvent *evt)
{
	_downPnt = Point(evt->x(), evt->y());
	_prePnt = _downPnt;
	checkMouseButton(evt->buttons());
}

void PaintWidget::mouseReleaseEvent(QMouseEvent *evt)
{
	if (_img3u.data == NULL || (_viewType != VIEW_OBJ && _viewType != VIEW_MARK))
		return;
	Mat objLabel1u = _viewType == VIEW_OBJ ? _objLabel1u : _inpLabel1u;
	if (objLabel1u.data == NULL)
		return;
	Point newPnt(evt->x(), evt->y());	
	checkMouseButton(Qt::MouseButtons(evt->button()));

	int markLabel = objLabel1u.at<byte>(_downPnt);
	if (_crntTool == TOOL_FILL){
		Mat fillMask = Mat::zeros(_h + 2, _w + 2, CV_8U);
		floodFill(objLabel1u, fillMask, newPnt, 255, NULL, 0, 0, FLOODFILL_MASK_ONLY);
		objLabel1u.setTo(markLabel, fillMask(Rect(1, 1, _w, _h)));
	}
	else if (_crntTool == TOOL_LINE)
		line(objLabel1u, _downPnt, newPnt, markLabel, _bSize, -1);
	viewStatic();
}

void PaintWidget::initializeGL()
{
	glClearColor(0, 0, 0, 1);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	setAutoFillBackground(false);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}


void PaintWidget::paintGL()
{
	if (_img3u.data == NULL)
		return;

	_objColorList.show(_viewType == VIEW_OBJ);
	
	// Drawing 2D scene
	glMatrixMode (GL_PROJECTION); 
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, _w, _h, 0, -1, 1);
	//_bgShow.setDispaly();
	paint2D();
	glPopMatrix();

	/* Drawing 3D scene
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_TEXTURE_2D);
	glClear(GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	paint3D();
	glPopMatrix();
	renderText(100, 100, tr("OK"));//*/	
}

void PaintWidget::paint2D()
{
	// Clear buffers and show background
	bool showBefore;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // 
	glDisable(GL_DEPTH_TEST);
	_bgShow3u.display();

	if (_viewType == VIEW_RES){
		for (int i = 0; i <= _crntAct; i++)
			if (i == _crntAct || _acts[i]._skipWhenLater == 0)
				_acts[i].Draw();
		for (int i = _crntAct + 1; i < _acts.size(); i++)
			if (_acts[i]._showBeforeAction)
				_acts[i].Draw();
	}
}

void PaintWidget::viewStatic(Point pnt, int t)
{
	if (_img3u.data == NULL || _viewType == VIEW_SRC || _viewType == VIEW_RES)
		return ;

	double a = 0.7, aI = 0.3; // Ratio of mark and image illustration
	Mat objIdx1u = _viewType == VIEW_ATTRI ? _attrMasks[_crntViewAtri] : (_viewType == VIEW_MARK ? _inpLabel1u : _objLabel1u);
	Mat uMarkH1u = _viewType == VIEW_ATTRI ? _attrMasks[_crntViewAtri] : _tgtAlpha1u;
	CV_Assert(objIdx1u.size == _img3u.size && objIdx1u.type() == CV_8UC1);
	ActionData::_colorIdx.idx2Color(objIdx1u, _staticShow3u);
	bool nvObj = _viewType != VIEW_OBJ;

#pragma omp parallel for
	for (int r = 0; r < _h; r++){
		Vec3b *show3u = _staticShow3u.ptr<Vec3b>(r);
		const Vec3b *img3u = _img3u.ptr<Vec3b>(r);
		const byte* hole1u = uMarkH1u.ptr<byte>(r);
		for (int c = 0; c < _w; c++)
			show3u[c] = aI * img3u[c] + a * (hole1u[c] && nvObj ? img3u[c] : show3u[c]);
	}

	Scalar colr = Scalar(ActionData::_colorIdx[objIdx1u.at<byte>(_downPnt)]);// _uMark3u.at<Vec3b>(_downPnt);
	if ((t & SHOW_BRUSH) && pnt.x != 0 && pnt.y != 0)
		circle(_staticShow3u, pnt, _bSize, colr);
	if (t & SHOW_LINE)
		line(_staticShow3u, _downPnt, pnt, colr, _bSize, _bSize == 1 ? 8 : -1);
	
	_bgShow3u.setImage(_staticShow3u);
	update();
}

void PaintWidget::setToolType(ToolType t)
{
	_ui.actionLine->setChecked(t == TOOL_LINE);
	_ui.actionBrush->setChecked(t == TOOL_BRUSH);
	_ui.actionFill->setChecked(t == TOOL_FILL);
	_crntTool = t;
}

void PaintWidget::setViewType(ViewType t)
{
	_ui.actionImage->setChecked(t == VIEW_SRC);
	_ui.actionObjects->setChecked(t == VIEW_OBJ);
	_ui.actionSoftTarget->setChecked(t == VIEW_SOFT_TGT);
	_ui.actionResults->setChecked(t == VIEW_RES);
	_ui.actionUserMark->setChecked(t == VIEW_MARK);

	if (_viewType == VIEW_MARK && t != VIEW_MARK)
		inpaintingBg(true);

	_viewType = t;
}

void PaintWidget::inpaintingBg(bool checkChanged)
{
	if (checkChanged){
		Mat dif1u;
		compare(_inpLabel1u, _objLabel1u, dif1u, CMP_NE);
		if (sum(dif1u).val[0] < 255*10)
			return; // Don't redo in-painting if user mark is little
	}
	printf("Background in-painting to be implemented\n");
}

void PaintWidget::paint3D()
{
	glTranslatef(0, 0, 1);
	glRotatef(0.5, 1, 1, 1);
	glColor3f(1, 0.6, 0);
	glutSolidTeapot(1); 
}

void PaintWidget::resizeGL(int w, int h)
{
	glViewport (0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90, (float)_w/_h, 0.01, 100.0);
	gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

void PaintWidget::fileOpen()
{
	string path = CmFile::BrowseFile("Images(*.jpg)\0*.jpg\0");
	if (!loadData(path))
		return;
	
	setMinimumSize(_w, _h);
	setMaximumSize(_w, _h);
	setViewType(VIEW_SRC);
	update();
}

bool PaintWidget::popWarning(const QString &wType, const QString &str)
{
	QMessageBox::warning(this, wType, str, QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Discard);
	return false;
}


bool PaintWidget::loadData(CStr &fPath)
{
	_timerAnimate.stop();
	_ui.actionFinished->setEnabled(false);

	// Load source and background image
	_nameNE = CmFile::GetNameNE(fPath);
	_dir = CmFile::GetFolder(fPath);
	_dirData = _dir + "Data/";
	CStr pathNE = _dirData + _nameNE;
	_img3u = imread(_dir + _nameNE + ".jpg");
	if (_img3u.data == NULL)
		return popWarning(tr("Fail to load file"), _QS(_dir + _nameNE + ".jpg"));
	_img3u.convertTo(_img3f, CV_32FC3, 1.0/255);
	_bgShow3u.setImage(_img3u);
	_w = _img3u.cols, _h = _img3u.rows;
	if (_w % 4 != 0)
		popWarning("waring", "Displaying images with width not multiple of 4 might have some problem.");
	_bg3u = imread(pathNE + "_B.jpg");
	if (_bg3u.data == NULL)
		_bg3u = _img3u;
	loadIniObjLabel(pathNE);
	_objColorList.setObjColorList(_objLabel1u);
	_tgtAlpha1u = Mat::zeros(_img3u.size(), CV_8U);
	_depth1u = imread(_dir + _nameNE + ".png", CV_LOAD_IMAGE_GRAYSCALE);
	loadAttriTexture();

	// Load textures in yml file
	glEnable(GL_TEXTURE_2D);
	_fs.open(pathNE + ".yml", FileStorage::READ);
	Action::LoadActions(*this, _acts);
	_crntAct = -1;
	glDisable(GL_TEXTURE_2D);

	_statusMsg.setText(tr("File loaded."));
	_uiMain.setWindowTitle(_QS(_nameNE) + ".jpg - Image Spirit");
	return true;
}

void PaintWidget::loadAttriTexture()
{
	_attrMasks.resize(1);
	_attrMasks[0] = Mat::ones(_h, _w, CV_8UC1);
	_ui.comboBoxAttributes->clear();
	_ui.comboBoxAttributes->addItem("");
	int atriSz = (int)ActionData::_materials.size();
	for (int i = 1; i <= atriSz; i++) {
		string atrName = _dirData + _nameNE + cv::format("_%d.png", i-1);
		if (!CmFile::FileExist(atrName))
			continue;
		_ui.comboBoxAttributes->addItem(QString::fromStdString(ActionData::_materials[i]));
		_attrMasks.push_back(imread(atrName, CV_LOAD_IMAGE_GRAYSCALE));
	}
}


void PaintWidget::loadIniObjLabel(CStr &pathNE)
{
	Mat objLabel3u = imread(pathNE + "_OA.png");
	if (objLabel3u.data == NULL){
		popWarning(tr("Can't load automatic object class map"), _QS(pathNE + "_OA.png"));
		return;
	}
	ActionData::_colorIdx.color2Idx(objLabel3u, _objLabel1u);
}

void PaintWidget::fileSave()
{
	_statusMsg.setText(tr("save"));
	paintGL();
	glFlush();
	QImage dispImg = grabFrameBuffer(true);
	dispImg.save(_QS(_dir + _nameNE + "_Res.png"), "png");
}

void PaintWidget::viewImage()
{
	setViewType(VIEW_SRC);
	_bgShow3u.setImage(_img3u);
	update();
}

void PaintWidget::viewObjects()
{
	setViewType(VIEW_OBJ);
	viewStatic();
}

void PaintWidget::viewSoftTarget()
{
	setViewType(VIEW_SOFT_TGT);
	printf("To be implemented\n");
}

void PaintWidget::viewResults()
{
	setViewType(VIEW_RES);
	_bgShow3u.setImage(_bg3u);
	update();
}

void PaintWidget::viewMark()
{
	setViewType(VIEW_MARK);
	_objLabel1u.copyTo(_inpLabel1u);
	viewStatic();
}

void PaintWidget::viewAttributes(int idx)
{
	if (idx >= 0)
		setViewType(VIEW_ATTRI);
	else 
		idx = 0;
	CV_Assert(idx < (int)(_attrMasks.size()));
	_crntViewAtri = idx;
	viewStatic();
}

void PaintWidget::toolEdit()
{
	if (_acts.size() == 0)
		return;
	_crntAct++;
	if (_crntAct < _acts.size()){
		_actType = ActionData::parseCmd(_acts[_crntAct]._scriptStr, _objDesM, _objDesA); //ActionData::printCmd(_actType, _objDesM, _objDesA);
		analysisTargets();

		_statusMsg.setText(_acts[_crntAct].Script());
		Sleep(800 + rand()%400);
		if (_acts[_crntAct]._actType == Action::NONE)
			return;

		if (_acts[_crntAct].IsAnimate()){
			_timerAnimate.start(_acts[_crntAct]._timerInterval);
			_ui.actionFinished->setEnabled(true);
		}
		else
			animateObjs();
	}
	else 
		_crntAct = -1;
}

void PaintWidget::toolFinished()
{
	if (_crntAct >= 0 && _crntAct < _acts.size())
		printf("Finished action `%s' at time %d\n", _acts[_crntAct].script(), _acts[_crntAct]._t);
	_timerAnimate.stop();
	//_statusMsg.setText(tr("Finished. Waiting for next speech command.\t"));
	_ui.actionFinished->setEnabled(false);
}

void PaintWidget::toolInpaint()
{
	string nameNE = _dirData + _nameNE;
	Mat img = imread(nameNE + "_B.png");
	if (img.data != NULL)	{
		Sleep(3000);
		flip(img, img, 0);
		img.copyTo(_staticShow3u);
		_viewType = VIEW_MARK;
		updateGL();
	}
}

void PaintWidget::toolExportOldData()
{
	string dir = _dir + "Data/";
	string outDir = dir + "New/";
	printf("Exporting old data in %s\n", _S(dir));
	vecS names;
	int imgNum = CmFile::GetNamesNE(dir + "*_O?.png", names);

	vecS &objNames = _objColorList._objNames;
	vecI &objColorIdx = _uiMain._objColorList._objColorIdx;
	vector<Vec3b> newColrs(objNames.size()); 

	CV_Assert(objNames.size() == ActionData::_objects.size());

	CmColorIdx clrIdx;
	map<int, int> oldNewIdx;
	for (int i = 0; i < objNames.size(); i++){
		newColrs[i] = clrIdx[i]; 
		CV_Assert(objNames[i] == ActionData::_objects[i]);
		oldNewIdx[objColorIdx[i]] = i;
	}

	for (int i = 0; i < imgNum; i++){
		Mat img = imread(dir + names[i] + ".png");
		for (int r = 0; r < img.rows; r++){
			Vec3b* dataP = img.ptr<Vec3b>(r);
			for (int c = 0; c < img.cols; c++){
				int dataIdx = dataP[c][2] * 1000000 + dataP[c][1] * 1000 + dataP[c][0];
				dataP[c] = newColrs[oldNewIdx[dataIdx]];
			}
		}
		imwrite(outDir + names[i] + ".png", img);
	}
	printf("finished\n");
}

void PaintWidget::toolLine()
{
	setToolType(TOOL_LINE);
}

void PaintWidget::toolBrush()
{
	setToolType(TOOL_BRUSH);

}

void PaintWidget::toolFill()
{
	setToolType(TOOL_FILL);
}

void PaintWidget::animateObjs()
{
	if (_crntAct >= 0 && _crntAct < _acts.size())
		_acts[_crntAct].Next(*this);
	//for (size_t i = 0; i < _acts.size(); i++)
	//	if (_acts[i]._actType == Action::DISPLAY && i != _crntAct)
	//		_acts[i].Next(*this);

	viewResults();
}



//void PaintWidget::paintEvent(QPaintEvent* event)
//{
//
//	if (_imgSrc.isNull())
//		return;
//
//	_uiMain._objColorList.Show(_viewType == VIEW_OBJ);
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//	glEnable(GL_DEPTH_TEST);
//
//	// Drawing 2D scene
//	glMatrixMode (GL_PROJECTION); 
//	glPushMatrix();
//	glLoadIdentity();
//	glOrtho(0, _w, _h, 0, -1, 1);
//	Paint2D();
//	glPopMatrix();
//
//	/* Drawing 3D scene
//	glMatrixMode(GL_MODELVIEW);
//	glPushMatrix();
//	Paint3D();
//	glPopMatrix();//*/
//
//	// Drawing text descriptions
//	QPainter painter(this);
//	//painter.setRenderHint(QPainter::Antialiasing);
//	//PaintText(painter);
//	painter.end(); 
//}

//void PaintWidget::drawInstructions(QPainter *painter, const QString &str, QPoint pos, int winLen)
//{
//	QFontMetrics metrics = QFontMetrics(font());
//	int border = metrics.leading();
//	QRect rect = metrics.boundingRect(pos.x(), pos.y(), winLen - 2*border, 50, Qt::AlignCenter | Qt::TextWordWrap, str);
//	painter->setRenderHint(QPainter::TextAntialiasing);
//	painter->fillRect(QRect(pos.x(), pos.y(), winLen, rect.height() + 2*border), QColor(0, 0, 0, 127));
//	painter->setPen(Qt::white);
//	painter->drawText(pos.x() + (winLen - rect.width())/2, pos.y() + border, rect.width(), rect.height(), Qt::AlignCenter | Qt::TextWordWrap, str);
//}

//void PaintWidget::PaintText(QPainter &painter)
//{
//	QString text = tr("Click and drag \n with the left");
//	drawInstructions(&painter, text, QPoint(10, 100), 200);
//}
