#-------------------------------------------------
#
# Project created by QtCreator 2010-08-23T12:40:57
#
#-------------------------------------------------

TEMPLATE = app


QT       += core
QT       -= gui

TARGET = Brown_retracker  # the name of the executable that will be produced
CONFIG += console


unix:!macx{
    INCLUDEPATH += \
                   "$$PWD" \


    LIBS += \
            -L"$$PWD/lib/levmar" -lm -llevmar \
            #-L"$$PWD/lib/lapack" -llapack \
            #-L"$$PWD/lib/blas" -lblas \
            -L"$$PWD/lib/openblas" -lblas -llapack \

}


macx{
    CONFIG += x86_64 Cocoa

    INCLUDEPATH += /usr/local/include
    LIBS += -framework Accelerate -L/usr/local/lib
}

HEADERS = \
    retrack_func/retracker_functions.h \
    misc/misc.h \
    Brown/Brown.h

SOURCES = \
    examples/main.cpp \
    Brown/Brown.cpp \
    retrack_func/retracker_functions.cpp \
    misc/misc.cpp \

    
    

