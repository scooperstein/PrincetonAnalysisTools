import ROOT
import sys
import numpy

ifile = ROOT.TFile(sys.argv[1])
tree = ifile.Get("tree")

TT = 1.01
Wj0b = 1.15
Wj1b = 2.08
Wj2b = 1.05

nentries = tree.GetEntries()
ofile = ROOT.TFile(sys.argv[2],"RECREATE")
otree = tree.CloneTree(0)

CS_SF_new = numpy.zeros(1,float)
otree.Branch("CS_SF_new",CS_SF_new,"CS_SF_new/D")

#nentries = 10
print "nentries = %i" % nentries
for i in range(nentries):
    if (i % 10000 == 0): print "processing entry: %i" % i
    tree.GetEntry(i)
    si = tree.sampleIndex
    if (si==2201 or si==4401 or si==4501 or si==4601 or si==4701 or si==4801 or si==4901):
        CS_SF_new[0] = Wj1b
    elif (si==2202 or si==4402 or si==4502 or      si==4602 or si==4702 or si==4802 or si==4902):
        CS_SF_new[0] = Wj2b
    elif (si==2200 or si==4400 or si==4500 or si==4600 or si==4700 or si==4800 or si==4900):
        CS_SF_new[0] = Wj0b
    elif (si==120 or si==50 or si==51 or si==52):
        CS_SF_new[0] = TT
    else:
        CS_SF_new[0] = 1.0
    ofile.cd()
    otree.Fill()

otree.Write()

     
