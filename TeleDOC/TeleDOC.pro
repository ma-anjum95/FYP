#-------------------------------------------------
#
# Project created by QtCreator 2017-04-16T22:12:45
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = TeleDOC
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp \
    PPG_C++/point.cpp \
    PPG_C++/ppg_analysis.cpp \
    PPG_C++/ppg_lines.cpp \
    PPG_C++/kiss_fft.c \
    ppgworker.cpp \
    MAX117/max30102.cpp

HEADERS  += \
    mainwindow.h \
    qcustomplot.h \
    PPG_C++/ppg_analysis.h \
    PPG_C++/ppg_lines.h \
    PPG_C++/point.h \
    PPG_C++/kiss_fft.h \
    PPG_C++/_kiss_fft_guts.h \
    ppgworker.h \
    MAX117/max30102.h

FORMS    += mainwindow.ui

RESOURCES += \
    images.qrc

LIBS += -lbcm2835 -lwiringPi

QT += testlib

QMAKE_CXXFLAGS += -std=gnu++11
