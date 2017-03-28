from ROOT import *
import sys
import numpy
from math import sqrt,pow
from array import array

def deltaR(phi1,phi2,eta1,eta2):
    dphi = abs(phi1-phi2)
    deta = abs(eta1-eta2)
    dr = sqrt(pow(dphi,2) + pow(deta,2))
    return dr

ifile = TFile(sys.argv[1],"r")
tree = ifile.Get("tree")

ofile = TFile(sys.argv[2],"RECREATE")
otree = tree.CloneTree(0)

for key in ifile.GetListOfKeys():
     if key.GetName() == 'tree':
         continue
     obj = key.ReadObj()
     obj.Write(key.GetName())

maxnum = 10
n = array( 'i', [ 0 ] )
maxBBInd = numpy.zeros(1, dtype=int)
#nFatjetAK08ungroomed = 0
#otree.SetBranchAddress("nFatjetAK08ungroomed",nFatjetAK08ungroomed)
FatjetAK08ungroomed_genHDR = numpy.zeros(10, dtype=float)
#FatjetAK08ungroomed_genHDR = array( 'f', maxnum*[ 0. ] )

otree.Branch( 'mynum', n, 'mynum/I' )
otree.Branch("maxBBInd",maxBBInd,"maxBBInd/I")
otree.Branch("FatjetAK08ungroomed_genHDR",FatjetAK08ungroomed_genHDR,"FatjetAK08ungroomed_genHDR[10]/D")
#otree.Branch("FatjetAK08ungroomed_genHDR",FatjetAK08ungroomed_genHDR,"FatjetAK08ungroomed_genHDR[mynum]/F")

nentries =  tree.GetEntries()
#nentries = 100
print "total entries: %i" % nentries

for i in range(nentries):
    if (i%10000==0):
        print "processing entry: %i" % i
    tree.GetEntry(i)
    maxBBInd[0] = -1
    maxBBVal = -10.
    for j in range(tree.nFatjetAK08ungroomed):
        val = tree.FatjetAK08ungroomed_bbtag[j]
        if (val > maxBBVal):
            maxBBVal = val
            maxBBInd[0] = j
        FatjetAK08ungroomed_genHDR[j] = deltaR(tree.FatjetAK08ungroomed_phi[j],tree.GenBJJ_phi,tree.FatjetAK08ungroomed_eta[j],tree.GenBJJ_eta)
        #print tree.FatjetAK08ungroomed_phi[j],tree.GenBJJ_phi,tree.FatjetAK08ungroomed_eta[j],tree.GenBJJ_eta,FatjetAK08ungroomed_genHDR[j]
    #print FatjetAK08ungroomed_genHDR
    ofile.cd()
    otree.Fill()
ofile.cd()
otree.Scan("FatjetAK08ungroomed_genHDR")
otree.Write()
