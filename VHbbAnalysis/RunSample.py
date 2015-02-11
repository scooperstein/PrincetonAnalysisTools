import ROOT
import sys
import ReadInput

if len(sys.argv) < 3:
    print "Please give two arguments:  the cfg file and the sample name"
    sys.exit(0)

# do stuff :)
ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
am=ReadInput.ReadTextFile(sys.argv[1], "cfg")
am.debug=2

if(am.debug>100):
    am.PrintBranches()

# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
#am.Loop()
am.LoopSample(sys.argv[2])
