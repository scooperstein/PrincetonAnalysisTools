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

# reads samples, existing branches and new branches
am=ReadInput.ReadTextFile(sys.argv[1], "cfg")
am.debug=2

if(am.debug>100):
    am.PrintBranches()

# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
am.Loop()

