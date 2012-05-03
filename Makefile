## For MacOS, uncomment OSTYPE, otherwise a *nix environment is assumed
##    This is in place to deal with the fact that Macs don't support
##    clock_gettime(CLOCK_REALTIME, &timeout);
## OSTYPE=macos

## To enable multicore:
## 1.  Explicity remove the next comment to enable the Butenhof workqueue
## 2.  Download http://homepage.mac.com/dbutenhof/Threads/code/workq.c and
##              http://homepage.mac.com/dbutenhof/Threads/code/workq.h 
##     and put them in the ./src directory
## 3.  apply the patch by running "make patch"
##
## NOTE:  there may be legal implications to using the Butenhof code
##        without permission
##
## BUTENHOF=-DENABLE_MCWQ\ -DMAX_THREADS=6

ifeq ($(OSTYPE),macos)
  LIB=-lgsl -lgslcblas -lpthread
else
  LIB=-lgsl -lgslcblas -lpthread -lrt
endif

## NOTE:  set this to the directory that contains Rinternals.h
INC=-I$(shell R RHOME)/include

## Compiler & flags
CC=g++
CPPFLAGS=-g $(INC) -O2 -fPIC

## Binary file names
RAPINAME=FastMap

files := src/dataobject.cc \
src/fastmap.cc \
src/generaltree.cc \
src/graphalg.cc \
src/graph.cc \
src/gslmat.cc \
src/gslvec.cc \
src/phylo_rbtalg.cc \
src/phylo_util.cc \
src/rootedbtree.cc \
src/R_interface.cc \
src/SingleCoreWorkQueue.cc \
src/stat_util.cc \
src/test.cc \
src/unrootedbtree.cc \
src/WorkQueue.cc

ifdef BUTENHOF
files+=src/MultiCoreWorkQueue.cc src/workq.c
endif

all: clean
	MAKEFLAGS='CXXFLAGS=-Wall\ -g\ $(BUTENHOF)' R CMD SHLIB $(files) -o src/$(RAPINAME).so $(LIB)

patch:
	patch src/workq.c < src/patch.workq.c
	patch src/workq.h < src/patch.workq.h

unpatch:
	patch -R src/workq.c < src/patch.workq.c
	patch -R src/workq.h < src/patch.workq.h

clean:
	rm -rf $(BINNAME) src/*.o src/*.d src/*\.d\.* src/*.so src/*~ .RData

release_clean: clean
	rm -rf src/workq.c src/workq.h
