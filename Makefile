

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
HOMEDIR = ${PWD}

all: AnalysisDict.cxx AnalysisDict.so

AnalysisDict.cxx: HelperClasses/SampleContainer.h HelperClasses/InfoStructs.h HelperClasses/BDTInfo.h HelperClasses/BDTVariable.h HelperClasses/SystematicContainer.h HelperClasses/SFContainer.h AnalysisManager.h plugins/VHbbAnalysis.h plugins/VHbbTrigger.h HelperClasses/BTagCalibrationStandalone.h LinkDef.h
	rootcint -f $@ -c $^


AnalysisDict.so: AnalysisDict.cxx HelperClasses/SampleContainer.cc HelperClasses/BDTInfo.cc HelperClasses/BDTVariable.cc HelperClasses/SystematicContainer.cc HelperClasses/SFContainer.cc AnalysisManager.cc plugins/VHbbAnalysis.cc plugins/VHbbTrigger.cc HelperClasses/EquationSolver.h HelperClasses/BTagCalibrationStandalone.cpp
	g++ -shared -fPIC -Wall -Werror -o $@ ${ROOTFLAGS} ${ROOTLIBS} -lTMVA -I${HOMEDIR} $^

clean: 
	rm AnalysisDict.cxx
	rm AnalysisDict.h
	rm AnalysisDict.so


