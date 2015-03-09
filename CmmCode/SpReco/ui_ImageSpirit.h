/********************************************************************************
** Form generated from reading UI file 'ImageSpirit.ui'
**
** Created: Tue Apr 9 15:21:43 2013
**      by: Qt User Interface Compiler version 4.8.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGESPIRIT_H
#define UI_IMAGESPIRIT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImageSpiritClass
{
public:
    QAction *actionOpen;
    QAction *actionSave;
    QAction *actionImage;
    QAction *actionSoftTarget;
    QAction *actionObjects;
    QAction *actionResults;
    QAction *actionEdit;
    QAction *actionFinished;
    QAction *actionUserMark;
    QAction *actionInPainting;
    QAction *actionExportOldData;
    QAction *actionLine;
    QAction *actionBrush;
    QAction *actionFill;
    QWidget *centralWidget;
    QComboBox *comboBoxAttributes;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuView;
    QMenu *menuTools;
    QMenu *menuTools_2;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ImageSpiritClass)
    {
        if (ImageSpiritClass->objectName().isEmpty())
            ImageSpiritClass->setObjectName(QString::fromUtf8("ImageSpiritClass"));
        ImageSpiritClass->resize(600, 400);
        actionOpen = new QAction(ImageSpiritClass);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/ImageSpirit/Resources/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpen->setIcon(icon);
        actionSave = new QAction(ImageSpiritClass);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/ImageSpirit/Resources/save.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSave->setIcon(icon1);
        actionImage = new QAction(ImageSpiritClass);
        actionImage->setObjectName(QString::fromUtf8("actionImage"));
        actionImage->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Source.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionImage->setIcon(icon2);
        actionSoftTarget = new QAction(ImageSpiritClass);
        actionSoftTarget->setObjectName(QString::fromUtf8("actionSoftTarget"));
        actionSoftTarget->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Layout.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSoftTarget->setIcon(icon3);
        actionObjects = new QAction(ImageSpiritClass);
        actionObjects->setObjectName(QString::fromUtf8("actionObjects"));
        actionObjects->setCheckable(true);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Objects.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionObjects->setIcon(icon4);
        actionResults = new QAction(ImageSpiritClass);
        actionResults->setObjectName(QString::fromUtf8("actionResults"));
        actionResults->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Attri.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionResults->setIcon(icon5);
        actionEdit = new QAction(ImageSpiritClass);
        actionEdit->setObjectName(QString::fromUtf8("actionEdit"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Edit.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionEdit->setIcon(icon6);
        actionFinished = new QAction(ImageSpiritClass);
        actionFinished->setObjectName(QString::fromUtf8("actionFinished"));
        actionFinished->setEnabled(false);
        QIcon icon7;
        icon7.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Stop.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFinished->setIcon(icon7);
        actionUserMark = new QAction(ImageSpiritClass);
        actionUserMark->setObjectName(QString::fromUtf8("actionUserMark"));
        actionUserMark->setCheckable(true);
        QIcon icon8;
        icon8.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Mark.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionUserMark->setIcon(icon8);
        actionInPainting = new QAction(ImageSpiritClass);
        actionInPainting->setObjectName(QString::fromUtf8("actionInPainting"));
        QIcon icon9;
        icon9.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Analysis.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionInPainting->setIcon(icon9);
        actionExportOldData = new QAction(ImageSpiritClass);
        actionExportOldData->setObjectName(QString::fromUtf8("actionExportOldData"));
        actionLine = new QAction(ImageSpiritClass);
        actionLine->setObjectName(QString::fromUtf8("actionLine"));
        actionLine->setCheckable(true);
        QIcon icon10;
        icon10.addFile(QString::fromUtf8(":/ImageSpirit/Resources/Line.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLine->setIcon(icon10);
        actionBrush = new QAction(ImageSpiritClass);
        actionBrush->setObjectName(QString::fromUtf8("actionBrush"));
        actionBrush->setCheckable(true);
        QIcon icon11;
        icon11.addFile(QString::fromUtf8(":/ImageSpirit/Resources/brush.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionBrush->setIcon(icon11);
        actionFill = new QAction(ImageSpiritClass);
        actionFill->setObjectName(QString::fromUtf8("actionFill"));
        actionFill->setCheckable(true);
        QIcon icon12;
        icon12.addFile(QString::fromUtf8(":/ImageSpirit/Resources/fill.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFill->setIcon(icon12);
        centralWidget = new QWidget(ImageSpiritClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        comboBoxAttributes = new QComboBox(centralWidget);
        comboBoxAttributes->setObjectName(QString::fromUtf8("comboBoxAttributes"));
        comboBoxAttributes->setGeometry(QRect(10, 10, 121, 22));
        comboBoxAttributes->setMinimumSize(QSize(100, 0));
        ImageSpiritClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ImageSpiritClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 600, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuView = new QMenu(menuBar);
        menuView->setObjectName(QString::fromUtf8("menuView"));
        menuTools = new QMenu(menuBar);
        menuTools->setObjectName(QString::fromUtf8("menuTools"));
        menuTools_2 = new QMenu(menuBar);
        menuTools_2->setObjectName(QString::fromUtf8("menuTools_2"));
        ImageSpiritClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ImageSpiritClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        ImageSpiritClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ImageSpiritClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        ImageSpiritClass->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuTools->menuAction());
        menuBar->addAction(menuTools_2->menuAction());
        menuBar->addAction(menuView->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionSave);
        menuView->addAction(actionImage);
        menuView->addAction(actionObjects);
        menuView->addAction(actionSoftTarget);
        menuView->addAction(actionResults);
        menuView->addAction(actionUserMark);
        menuTools->addAction(actionEdit);
        menuTools->addAction(actionFinished);
        menuTools->addAction(actionInPainting);
        menuTools->addAction(actionExportOldData);
        menuTools_2->addAction(actionLine);
        menuTools_2->addAction(actionBrush);
        menuTools_2->addAction(actionFill);
        mainToolBar->addAction(actionOpen);
        mainToolBar->addAction(actionSave);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionEdit);
        mainToolBar->addAction(actionFinished);
        mainToolBar->addAction(actionInPainting);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionLine);
        mainToolBar->addAction(actionBrush);
        mainToolBar->addAction(actionFill);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionImage);
        mainToolBar->addAction(actionObjects);
        mainToolBar->addAction(actionSoftTarget);
        mainToolBar->addAction(actionResults);
        mainToolBar->addAction(actionUserMark);

        retranslateUi(ImageSpiritClass);

        QMetaObject::connectSlotsByName(ImageSpiritClass);
    } // setupUi

    void retranslateUi(QMainWindow *ImageSpiritClass)
    {
        ImageSpiritClass->setWindowTitle(QApplication::translate("ImageSpiritClass", "ImageSpirit", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("ImageSpiritClass", "Open", 0, QApplication::UnicodeUTF8));
        actionOpen->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("ImageSpiritClass", "Save", 0, QApplication::UnicodeUTF8));
        actionSave->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionImage->setText(QApplication::translate("ImageSpiritClass", "Image", 0, QApplication::UnicodeUTF8));
        actionImage->setShortcut(QApplication::translate("ImageSpiritClass", "Alt+I", 0, QApplication::UnicodeUTF8));
        actionSoftTarget->setText(QApplication::translate("ImageSpiritClass", "Soft target", 0, QApplication::UnicodeUTF8));
        actionSoftTarget->setShortcut(QApplication::translate("ImageSpiritClass", "Alt+S", 0, QApplication::UnicodeUTF8));
        actionObjects->setText(QApplication::translate("ImageSpiritClass", "Objects", 0, QApplication::UnicodeUTF8));
        actionObjects->setShortcut(QApplication::translate("ImageSpiritClass", "Alt+O", 0, QApplication::UnicodeUTF8));
        actionResults->setText(QApplication::translate("ImageSpiritClass", "Results", 0, QApplication::UnicodeUTF8));
        actionResults->setShortcut(QApplication::translate("ImageSpiritClass", "Alt+A", 0, QApplication::UnicodeUTF8));
        actionEdit->setText(QApplication::translate("ImageSpiritClass", "Edit", 0, QApplication::UnicodeUTF8));
        actionEdit->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+E", 0, QApplication::UnicodeUTF8));
        actionFinished->setText(QApplication::translate("ImageSpiritClass", "Finished", 0, QApplication::UnicodeUTF8));
        actionFinished->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+F", 0, QApplication::UnicodeUTF8));
        actionUserMark->setText(QApplication::translate("ImageSpiritClass", "User Mark", 0, QApplication::UnicodeUTF8));
        actionUserMark->setShortcut(QApplication::translate("ImageSpiritClass", "Alt+M", 0, QApplication::UnicodeUTF8));
        actionInPainting->setText(QApplication::translate("ImageSpiritClass", "In-Painting", 0, QApplication::UnicodeUTF8));
        actionInPainting->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+I", 0, QApplication::UnicodeUTF8));
        actionExportOldData->setText(QApplication::translate("ImageSpiritClass", "ExportOldData", 0, QApplication::UnicodeUTF8));
        actionLine->setText(QApplication::translate("ImageSpiritClass", "Line", 0, QApplication::UnicodeUTF8));
        actionLine->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+L", 0, QApplication::UnicodeUTF8));
        actionBrush->setText(QApplication::translate("ImageSpiritClass", "Brush", 0, QApplication::UnicodeUTF8));
        actionBrush->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+B", 0, QApplication::UnicodeUTF8));
        actionFill->setText(QApplication::translate("ImageSpiritClass", "Fill", 0, QApplication::UnicodeUTF8));
        actionFill->setShortcut(QApplication::translate("ImageSpiritClass", "Ctrl+F", 0, QApplication::UnicodeUTF8));
        comboBoxAttributes->clear();
        comboBoxAttributes->insertItems(0, QStringList()
         << QString()
        );
        menuFile->setTitle(QApplication::translate("ImageSpiritClass", "File", 0, QApplication::UnicodeUTF8));
        menuView->setTitle(QApplication::translate("ImageSpiritClass", "View", 0, QApplication::UnicodeUTF8));
        menuTools->setTitle(QApplication::translate("ImageSpiritClass", "Edit", 0, QApplication::UnicodeUTF8));
        menuTools_2->setTitle(QApplication::translate("ImageSpiritClass", "Tools", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ImageSpiritClass: public Ui_ImageSpiritClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGESPIRIT_H
