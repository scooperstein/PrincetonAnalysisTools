import ROOT
import sys
from math import sqrt,exp
import numpy

def getWeight(pt1,pt2):
    sf1 = exp(0.0615 - 0.0005*pt1)
    sf2 = exp(0.0615 - 0.0005*pt2)
    return sqrt(sf1*sf2)
    

ifile = ROOT.TFile(sys.argv[1])
ofile = ROOT.TFile(sys.argv[2],"RECREATE")

tree = ifile.Get("tree")
nentries = tree.GetEntries()
#nentries = 100000

otree = tree.CloneTree(0)
topPtWeight = numpy.zeros(1,dtype=float)
otree.Branch("topPtWeight",topPtWeight,"topPtWeight/D")

print "total entries: ",nentries
for i in range(nentries):
    if (i%10000 == 0): print "entry: %i" %i
    tree.GetEntry(i)
    topPtWeight[0] = getWeight(tree.GenTop_pt[0],tree.GenTop_pt[1])
    otree.Fill()

ofile.cd()
otree.Write()

for key in ifile.GetListOfKeys():
     if key.GetName() == 'tree':
         continue
     obj = key.ReadObj()
     obj.Write(key.GetName())
