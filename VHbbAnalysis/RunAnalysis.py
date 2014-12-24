import ROOT
import sys
import ReadInput

# do stuff :)
# FIXME - link all libraries with proper Makefile
ROOT.gSystem.Load("SampleContainer_cc.so")
ROOT.gSystem.Load("BDTInfo_h.so")
ROOT.gSystem.Load("InfoStructs_h.so")
ROOT.gSystem.Load("VHbbAnalysis_h.so")
ROOT.gSystem.Load("AnalysisManager_cc.so")


# the file passed to the constructor is just the template for all the other files. For now all added files will need to match this structure.
vhbba = ROOT.VHbbAnalysis()
vhbba.Initialize("WH_125_lumiWeighted.root")
vhbba.debug=2


# reads samples, existing branches and new branches
ReadInput.ReadTextFile("vhbb_config.txt", "cfg", vhbba)

if(vhbba.debug>100):
    vhbba.PrintBranches()

# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
vhbba.Loop()
