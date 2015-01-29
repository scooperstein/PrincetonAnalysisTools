

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
HOMEDIR = ${PWD}

all: AnalysisDict.cxx AnalysisDict.so

AnalysisDict.cxx: HelperClasses/SampleContainer.h HelperClasses/InfoStructs.h HelperClasses/BDTInfo.h AnalysisManager.h plugins/VHbbAnalysis.h plugins/VHbbTrigger.h LinkDef.h
	rootcint -f $@ -c $^

AnalysisDict.so: AnalysisDict.cxx HelperClasses/SampleContainer.cc HelperClasses/BDTInfo.h AnalysisManager.cc plugins/VHbbAnalysis.cc plugins/VHbbTrigger.cc
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -lTMVA -I$(ROOTSYS)/include -I${HOMEDIR} $^

clean: 
	rm AnalysisDict.cxx
	rm AnalysisDict.h
	rm AnalysisDict.so

