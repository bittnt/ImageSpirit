#include "stdafx.h"
#include "ImageSpirit.h"


#include "JointCRF.h"

void demoAtriCRF();

int main(int argc, char *argv[])
{
	return JointCRF::demo(argc, argv);

	glutInit(&argc, argv);
	QApplication a(argc, argv);
	ImageSpirit w;
	w.setWindowTitle(QString("Image Spirit"));
	w.setGeometry(8, 30, 784, 562);
	w.show();
	return a.exec();
}
