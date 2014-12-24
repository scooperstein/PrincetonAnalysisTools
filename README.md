PrincetonAnalysisTools
======================

source env.sh


#FIXME - introduce a Makefile to build these libraries
For now, make .so libraries by hand in ROOT

.L HelperClasses/SampleContainer.cc++

.L HelperClasses/InfoStructs.h++

.L HelperClasses/BDTInfo.h++


.L AnalysisManager.cc++

.L plugins/VHbbAnalysis.cc++

.L plugins/VHbbAnalysis.h++


cd VHbbAnalysis

python RunAnalysis.py
