

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
HOMEDIR = ${PWD}

all: AnalysisDict.cxx AnalysisDict.so

AnalysisDict.cxx: HelperClasses/SampleContainer.h HelperClasses/InfoStructs.h HelperClasses/BDTInfo.h AnalysisManager.h plugins/VHbbAnalysis.h plugins/VHbbTrigger.h plugins/WorkspaceAnalysis.h LinkDef.h
	rootcint -f $@ -c $^

AnalysisDict.so: AnalysisDict.cxx HelperClasses/SampleContainer.cc HelperClasses/BDTInfo.cc AnalysisManager.cc plugins/VHbbAnalysis.cc plugins/VHbbTrigger.cc plugins/WorkspaceAnalysis.cc
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -lTMVA -lRooFitCore -I${HOMEDIR} $^

clean: 
	rm AnalysisDict.cxx
	rm AnalysisDict.h
	rm AnalysisDict.so


