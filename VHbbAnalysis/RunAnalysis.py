import ROOT
import sys

sys.path.insert(0, "python")
import ReadInput

# do stuff :)
ROOT.gSystem.Load("../HelperClasses/SampleContainer_cc.so")
ROOT.gSystem.Load("../AnalysisManager_cc.so")

# the file passed to the constructor is just the template for all the other files. For now all added files will need to match this structure.
am = ROOT.AnalysisManager("WH_125_lumiWeighted.root")
am.debug=2

# read branches already in the tree
branchdic=ReadInput.ReadTextFile("existingbranches.txt","branchlist")

for branch in branchdic:
    am.SetupBranch(branch,branchdic[branch][0], branchdic[branch][1])

sampledic=ReadInput.ReadTextFile("samples.txt","samplefile")

# ConfigureOutputTree() needs to be called between adding the current branches and adding the new branches
am.ConfigureOutputTree()
# read new branches that we want to add to the tree
newbranchdic=ReadInput.ReadTextFile("newbranches.txt","branchlist")
#
for branch in newbranchdic:
    print branch
    am.SetupNewBranch(branch,newbranchdic[branch][0], newbranchdic[branch][1])

for sample in sampledic:
    #print sample, sampledic[sample]
    samplecon = ROOT.SampleContainer()
    samplecon.sampleNum = sampledic[sample][0]
    samplecon.xsec      = sampledic[sample][1]
    samplecon.kfactor   = sampledic[sample][2]
    samplecon.scale     = sampledic[sample][3]
    for filename in sampledic[sample][4]:
        samplecon.AddFile(filename)
    am.AddSample(samplecon)


# setup analysis - in the future...
#am.analysis=ROOT.WHbbAnalysis()
#am.analysis=ROOT.Analysis("anlaysisname")
if(am.debug>100):
    am.PrintBranches()

# loop over all the samples. We'll probably want to add a method to loop over only a specific subset of samples
am.Loop()
