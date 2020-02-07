#-------------------------------------------------
#
# Project created by QtCreator 2020-01-26T13:20:04
#
#-------------------------------------------------

QT       -= gui
QT       += concurrent

TARGET = runIons
TEMPLATE = lib

DEFINES += RUNIONS_LIBRARY

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    pa.cpp \
    simionfield.cpp \
    timefield.cpp \
    phasestateintegrator.cpp \
    trackion.cpp \
    pythonfields.cpp \
    pyrunions.cpp \
    pyexport.cpp

HEADERS += \
    pa.h \
    task.h \
    util.h \
    maybe.h \
    field.h \
    simionfield.h \
    timefield.h \
    phasestateintegrator.h \
    trackion.h \
    pythonfields.h \
    pyrunions.h

INCLUDEPATH += \
    /usr/include/python3.4m

LIBS += -lboost_python-py34

unix {
    target.path = /usr/lib
    INSTALLS += target
}

QMAKE_CXXFLAGS += -std=c++11
