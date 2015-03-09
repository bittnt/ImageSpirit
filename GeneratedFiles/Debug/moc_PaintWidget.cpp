/****************************************************************************
** Meta object code from reading C++ file 'PaintWidget.h'
**
** Created: Wed Apr 10 10:54:00 2013
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "StdAfx.h"
#include "../../PaintWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PaintWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PaintWidget[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      24,   12,   12,   12, 0x08,
      35,   12,   12,   12, 0x08,
      46,   12,   12,   12, 0x08,
      61,   12,   12,   12, 0x08,
      75,   12,   12,   12, 0x08,
      95,   12,   12,   12, 0x08,
     106,   12,   12,   12, 0x08,
     118,   12,   12,   12, 0x08,
     129,   12,   12,   12, 0x08,
     141,   12,   12,   12, 0x08,
     158,   12,   12,   12, 0x08,
     172,   12,   12,   12, 0x08,
     190,  186,   12,   12, 0x08,
     210,   12,   12,   12, 0x08,
     221,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_PaintWidget[] = {
    "PaintWidget\0\0fileOpen()\0fileSave()\0"
    "toolEdit()\0toolFinished()\0toolInpaint()\0"
    "toolExportOldData()\0toolLine()\0"
    "toolBrush()\0toolFill()\0viewImage()\0"
    "viewSoftTarget()\0viewObjects()\0"
    "viewResults()\0idx\0viewAttributes(int)\0"
    "viewMark()\0animateObjs()\0"
};

void PaintWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        PaintWidget *_t = static_cast<PaintWidget *>(_o);
        switch (_id) {
        case 0: _t->fileOpen(); break;
        case 1: _t->fileSave(); break;
        case 2: _t->toolEdit(); break;
        case 3: _t->toolFinished(); break;
        case 4: _t->toolInpaint(); break;
        case 5: _t->toolExportOldData(); break;
        case 6: _t->toolLine(); break;
        case 7: _t->toolBrush(); break;
        case 8: _t->toolFill(); break;
        case 9: _t->viewImage(); break;
        case 10: _t->viewSoftTarget(); break;
        case 11: _t->viewObjects(); break;
        case 12: _t->viewResults(); break;
        case 13: _t->viewAttributes((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: _t->viewMark(); break;
        case 15: _t->animateObjs(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData PaintWidget::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject PaintWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_PaintWidget,
      qt_meta_data_PaintWidget, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PaintWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PaintWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PaintWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PaintWidget))
        return static_cast<void*>(const_cast< PaintWidget*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int PaintWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
