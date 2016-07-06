## Add weight which accounts for data/MC discrepancies in TT Powheg nAddJet distribution
## 
## Author: Stephane Cooperstein
##

import ROOT
import sys
import numpy

ifile = ROOT.TFile(sys.argv[1])
ofile = ROOT.TFile(sys.argv[2],"RECREATE")

tree = ifile.Get("tree")
otree = tree.CloneTree(0)

nentries = tree.GetEntries()
print "Total Entries: %i" % nentries

weight_map = {}
weight_map[0] = 1.053647
weight_map[1] = 1.042559
weight_map[2] = 1.032228
weight_map[3] = 1.048734
weight_map[4] = 0.933228
weight_map[5] = 0.844447
weight_map[6] = 0.992374
weight_map[7] = 0.857796
weight_map[8] = 0.796930
weight_map[9] = 0.937715

nAddJet_reweightUp = numpy.zeros(1, dtype=float)
nAddJet_reweightDown = numpy.zeros(1, dtype=float)
otree.Branch("nAddJet_reweightUp",nAddJet_reweightUp,"nAddJet_reweightUp/D")
otree.Branch("nAddJet_reweightDown",nAddJet_reweightDown,"nAddJet_reweightDown/D")

for i in range(nentries):
    if (i % 10000 == 0): 
        print "Processing entry: %i" % i
    tree.GetEntry(i)
    if (tree.sampleIndex != 120 or tree.nAddJet_f >=10):
        nAddJet_reweightUp[0] = 1.0
        nAddJet_reweightDown[0] = 1.0
    else:
        # sample is ttbar powheg
        nAddJet_reweightUp[0] = max(weight_map[int(tree.nAddJet_f)], 2.0 - weight_map[int(tree.nAddJet_f)])
        nAddJet_reweightDown[0] = min(weight_map[int(tree.nAddJet_f)], 2.0 - weight_map[int(tree.nAddJet_f)])
        #print weight_map[int(tree.nAddJet_f)]
    otree.Fill()

ofile.cd()
otree.Write()
