

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
HOMEDIR = /uscms_data/d3/sbc01/HbbAnalysis13TeV/PrincetonAnalysisTools

all: HelperClasses/SampleContainer_cc.so HelperClasses/InfoStructs_h.so HelperClasses/BDTInfo_h.so AnalysisManager_cc.so plugins/VHbbAnalysis_cc.so

HelperClasses/SampleContainer_cc.cxx: HelperClasses/SampleContainer.h HelperClasses/SampleContainer_LinkDef.h
	rootcint -f $@ -c $^

HelperClasses/SampleContainer_cc.so: HelperClasses/SampleContainer_cc.cxx HelperClasses/SampleContainer.cc
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm HelperClasses/SampleContainer_cc.cxx

HelperClasses/InfoStructs_h.cxx: HelperClasses/InfoStructs.h HelperClasses/InfoStructs_LinkDef.h
	rootcint -f $@ -c  -p $^

HelperClasses/InfoStructs_h.so: HelperClasses/InfoStructs_h.cxx HelperClasses/InfoStructs.h
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm HelperClasses/InfoStructs_h.cxx

HelperClasses/BDTInfo_h.cxx: HelperClasses/BDTInfo.h HelperClasses/BDTInfo_LinkDef.h
	rootcint -f $@ -c  -p $^

HelperClasses/BDTInfo_h.so: HelperClasses/BDTInfo_h.cxx HelperClasses/BDTInfo.h
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm HelperClasses/BDTInfo_h.cxx

AnalysisManager_cc.cxx: AnalysisManager.h AnalysisManager_LinkDef.h
	rootcint -f $@ -c  -p $^

AnalysisManager_cc.so: AnalysisManager_cc.cxx AnalysisManager.cc
	g++ -shared -fPIC -o $@ $(ROOTFLAGS) ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm AnalysisManager_cc.cxx

#plugins/VHbbAnalysis_h.so: plugins/VHbbAnalysis.h
#	g++ -shared -fPIC -o $@ plugins/VHbbAnalysis.h ${ROOTFLAGS} ${ROOTLIBS}

plugins/VHbbAnalysis_cc.cxx: plugins/VHbbAnalysis.h plugins/VHbbAnalysis_LinkDef.h
	rootcint -f $@ -c  -p $^

plugins/VHbbAnalysis_cc.so: plugins/VHbbAnalysis_cc.cxx plugins/VHbbAnalysis.cc
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm plugins/VHbbAnalysis_cc.cxx

#plugins/VHbbTrigger_h.so: plugins/VHbbTrigger.h
#	g++ -shared -fPIC -o $@ plugins/VHbbTrigger.h ${ROOTFLAGS} ${ROOTLIBS}

plugins/VHbbTrigger_cc.cxx: plugins/VHbbTrigger.h plugs/VHbbTrigger_LinkDef.h
	rootcint -f $@ -c  -p $^

plugins/VHbbTrigger_cc.so: plugins/VHbbTrigger_cc.cxx plugins/VHbbTrigger.cc
	g++ -shared -fPIC -o $@ ${ROOTFLAGS} ${ROOTLIBS} -I$(ROOTSYS)/include -I${HOMEDIR} $^
	rm plugins/VHbbTrigger_cc.cxx

clean:
	rm HelperClasses/SampleContainer_cc.so
	rm HelperClasses/SampleContainer_cc.cxx
	rm HelperClasses/InfoStructs_h.so
	rm HelperClasses/InfoStructs_h.cxx
	rm HelperClasses/BDTInfo_h.so
	rm HelperClasses/BDTInfo_h.cxx
	rm AnalysisManager_cc.so
	rm AnalysisManager_cc.cxx
	rm plugins/VHbbAnalysis_cc.so
	rm plugins/VHbbAnalysis_cc.cxx

