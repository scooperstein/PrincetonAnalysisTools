#!/bin/bash
##
## Script to run Analysis Manager jobs on Condor from the LPC
## Author: Stephane Cooperstein
##
## Argument 1: Analysis config file
## Argument 2: Sample name to run on
##

# Set up environment
echo "setting up the environment"
#pushd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_2_0_pre6/src
pushd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_4_7/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
pushd /uscms_data/d3/sbc01/HbbAnalysis13TeV/PrincetonAnalysisTools/
source env.sh
popd
popd 
#source /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.22-cms/bin/thisroot.sh
#source /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.18-cms/bin/thisroot.sh
echo "successfully set up the enviroment"

echo "moving text files to their respective directories"
mkdir -p cfg
mv samples.txt earlybranches.txt existingbranches.txt newbranches.txt bdtsettings.txt reg1_settings.txt reg2_settings.txt settings.txt cfg
mkdir -p aux
mv new-weights-23Jan.xml TMVA_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4_BDT.weights.xml aux


echo "running RunSample.py"
python RunSample.py $1 $2 $3 $4
ls
echo "all done!"
