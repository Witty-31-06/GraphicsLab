QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
DEFINES += QT_NO_KEYWORDS
SOURCES += \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    mainwindow.h

FORMS += \
    mainwindow.ui

# Intel MKL and Intel TBB paths (adjust paths according to your installation)
MKLROOT = /opt/intel/oneapi/mkl/latest
TBBROOT = /usr/local/
# Include paths for Intel MKL and TBB

INCLUDEPATH += /usr/include/eigen3
# Link paths for Intel MKL and TBB
LIBS += -L$$MKLROOT/lib/intel64 \
        -L$$TBBROOT/lib/

# Link against Intel MKL and TBB libraries
LIBS += -lmkl_rt -ltbb

# Additional linker flags (important for MKL)
LIBS += -lpthread -lm -ldl

# Enable OpenMP for parallelism
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

# Default rules for deployment
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
