APP_NAME      := calib

SOURCE_FILES  := main.cc  FitModule.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst -lSpectrum
APP_CXX_FLAGS += -O0 -std=c++11  -g -w

.PHONY:  default
default: build install

include $(HADDIR)/hades.app.mk

