import ROOT
import sys
import ReadInput

if len(sys.argv) < 2:
    print "Please give one argument:  the cfg file"
    sys.exit(0)

# do stuff :)
# FIXME - link all libraries with proper Makefile
ROOT.gSystem.Load("SampleContainer_cc.so")
ROOT.gSystem.Load("BDTInfo_h.so")
ROOT.gSystem.Load("InfoStructs_h.so")
ROOT.gSystem.Load("VHbbAnalysis_h.so")
ROOT.gSystem.Load("AnalysisManager_cc.so")


# the file passed to the constructor is just the template for all the other files. For now all added files will need to match this structure.
vhbba = ROOT.VHbbAnalysis()
vhbba.Initialize("WH_HToBB_WToLNu_M-125_13TeV_V4_tree.root")
vhbba.debug=2


# reads samples, existing branches and new branches
ReadInput.ReadTextFile(sys.argv[1], "cfg", vhbba)

if(vhbba.debug>100):
    vhbba.PrintBranches()

# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
vhbba.Loop()
