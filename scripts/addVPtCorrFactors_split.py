import sys
import ROOT
import numpy

def getCorrFactorTT(V_pt):
    return (1.064 - 0.000380*V_pt)
def getCorrFactorTTUp(V_pt):
    return (1.064 - 0.000469*V_pt)
def getCorrFactorTTDown(V_pt):
    return (1.064 - 0.000291*V_pt)
def getCorrFactorWLF(V_pt):
    return (1.097 - 0.000575*V_pt)
def getCorrFactorWLFUp(V_pt):
    return (1.097 - 0.000621*V_pt)
def getCorrFactorWLFDown(V_pt):
    return (1.097 - 0.000529*V_pt)
def getCorrFactorWHF(V_pt):
    return (1.259 - 0.00167*V_pt)
def getCorrFactorWHFUp(V_pt):
    return (1.259 - 0.00180*V_pt)
def getCorrFactorWHFDown(V_pt):
    return (1.259 - 0.00154*V_pt)

ifile = ROOT.TFile(sys.argv[1])
ofile = ROOT.TFile(sys.argv[2],"RECREATE")

tree = ifile.Get("tree")
nentries = tree.GetEntries()
#nentries = 1000

otree = tree.CloneTree(0)
VPtCorrFactor = numpy.zeros(1,dtype=float)
VPtCorrFactorUp = numpy.zeros(1,dtype=float)
VPtCorrFactorDown = numpy.zeros(1,dtype=float)
otree.Branch("VPtCorrFactorSplit3",VPtCorrFactor,"VPtCorrFactorSplit/D")
otree.Branch("VPtCorrFactorSplit3Up",VPtCorrFactorUp,"VPtCorrFactorSplit3Up/D")
otree.Branch("VPtCorrFactorSplit3Down",VPtCorrFactorDown,"VPtCorrFactorSplit3Down/D")

print "Processing %i entries" % nentries
for i in range(nentries):
    if (i%10000 == 0): print "Processing entry: %i" % i
    tree.GetEntry(i)
    si = tree.sampleIndex
    V_pt = tree.V_pt
    VPtCorrFactor[0] = 1.0
    VPtCorrFactorUp[0] = 1.0
    VPtCorrFactorDown[0] = 1.0
    if (si == 120 or si == 50 or si == 51 or si == 52):
        VPtCorrFactor[0] = getCorrFactorTT(V_pt)
        VPtCorrFactorUp[0] = getCorrFactorTTUp(V_pt)
        VPtCorrFactorDown[0] = getCorrFactorTTDown(V_pt)
    elif (((si>=4100 and  si<=4902) or (si>=2200 and si<=2202) or (si>=48100 and si<=49102)) and si%100 == 0):
        VPtCorrFactor[0] = getCorrFactorWLF(V_pt)
        VPtCorrFactorUp[0] = getCorrFactorWLFUp(V_pt)
        VPtCorrFactorDown[0] = getCorrFactorWLFDown(V_pt)
    elif (((si>=4100 and  si<=4902) or (si>=2200 and si<=2202) or (si>=48100 and si<=49102)) and (si%100 == 1 or si%100 == 2)):
        VPtCorrFactor[0] = getCorrFactorWHF(V_pt)
        VPtCorrFactorUp[0] = getCorrFactorWHFUp(V_pt)
        VPtCorrFactorDown[0] = getCorrFactorWHFDown(V_pt)
    elif (si == 16 or si == 17 or si == 20 or si == 21):
        VPtCorrFactor[0] = getCorrFactorWHF(V_pt)
        VPtCorrFactorUp[0] = getCorrFactorWHFUp(V_pt)
        VPtCorrFactorDown[0] = getCorrFactorWHFDown(V_pt)
    #print si,V_pt,VPtCorrFactor[0]
    #VPtCorrFactor[0] = getCorrFactor(V_pt)
    #VPtCorrFactorUp[0] = getCorrFactorUp(V_pt)
    #VPtCorrFactorDown[0] = getCorrFactorDown(V_pt)
    #VPtCorrFactor[0] = getCorrFactor(V_pt,1.0)
    otree.Fill()
    
ofile.cd()
otree.Write()

for key in ifile.GetListOfKeys():
     if key.GetName() == 'tree':
         continue
     obj = key.ReadObj()
     obj.Write(key.GetName())
