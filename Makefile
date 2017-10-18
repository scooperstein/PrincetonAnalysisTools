

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
HOMEDIR = ${PWD}

all: AnalysisDict.cxx AnalysisDict.so

AnalysisDict.cxx: HelperClasses/SampleContainer.h HelperClasses/InfoStructs.h HelperClasses/BDTInfo.h HelperClasses/BDTVariable.h HelperClasses/SystematicContainer.h HelperClasses/SFContainer.h HelperClasses/RoccoR.h AnalysisManager.h plugins/VHbbAnalysis.h plugins/VHbbTrigger.h plugins/EWKAnalysis.h plugins/EWKZAnalysis.h LinkDef.h
	rootcint -f $@ -c $^

AnalysisDict.so: AnalysisDict.cxx HelperClasses/SampleContainer.cc HelperClasses/BDTInfo.cc HelperClasses/BDTVariable.cc HelperClasses/SystematicContainer.cc HelperClasses/SFContainer.cc HelperClasses/RoccoR.cc AnalysisManager.cc plugins/VHbbAnalysis.cc plugins/VHbbTrigger.cc plugins/EWKAnalysis.cc plugins/EWKZAnalysis.cc HelperClasses/EquationSolver.h
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -lTMVA -I${HOMEDIR} $^

clean: 
	rm AnalysisDict.cxx
	rm AnalysisDict.h
	rm AnalysisDict.so


