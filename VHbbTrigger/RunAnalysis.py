import ROOT
import sys
import ReadInput

if len(sys.argv) < 2:
    print "Please give one argument:  the cfg file"
    sys.exit(0)

ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
am=ReadInput.ReadTextFile(sys.argv[1], "cfg")
am.debug=2

if(am.debug>100):
    am.PrintBranches()

# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
am.Loop()

