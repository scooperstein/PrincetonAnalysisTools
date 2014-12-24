PrincetonAnalysisTools
======================

source env.sh


#FIXME - introduce a Makefile to build these libraries
For now, make .so libraries by hand in ROOT

.L HelperClasses/SampleContainers.cc++

.L HelperClasses/BaseAnalysis.h++

.L HelperClasses/BDTInfo.h++


.L AnalysisManager.cc++

.L plugins/VHbbAnalysis.cc++


cd VHbbAnalysis

python RunAnalysis.py
