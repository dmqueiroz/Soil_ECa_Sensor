#-------------------------------------------------
#
# Project created by QtCreator 2016-10-28T14:25:45
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = My_Sensor_205
TEMPLATE = app


SOURCES += main.cpp\
           mysensor.cpp \          
           BlackCore.cpp \          
           BlackPWM.cpp \
           gps_receiver.cpp \
    Kriging_Lib.cpp

HEADERS  += mysensor.h \
            BlackLib.h \
            BlackCore.h \
            BlackDef.h \
            BlackErr.h \            
            BlackPWM.h \
            gps_receiver.h \
    Kriging_Lib.h

FORMS    += mysensor.ui

QMAKE_CXXFLAGS += -std=c++0x
