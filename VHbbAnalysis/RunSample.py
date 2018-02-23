import ROOT
import sys
import ReadInput

if (len(sys.argv) != 3 and len(sys.argv) != 5):
    print "Please give two arguments:  the cfg file and the sample name"
    print "Or give four arguments: the cfg file, the sample name, a comma-separated list of input files, and the name of the output root file"
    sys.exit(0)

# do stuff :)
ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
samplesToRun = []
samplesToRun.append(sys.argv[2])
filesToRun = []
if len(sys.argv) == 5:
    for item in sys.argv[3].split(','):
        filesToRun.append(item)

am=ReadInput.ReadTextFile(sys.argv[1], "cfg", samplesToRun, filesToRun)
am.debug=2

print "Read in the input files, now let's run it!"
if(am.debug>100):
    am.PrintBranches()

print "Done printing branches, now to loop"
# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
#aim.Loop()

if (len(sys.argv) == 3):
    am.Loop(sys.argv[2])
else:
    am.Loop(sys.argv[2], sys.argv[3], sys.argv[4])
