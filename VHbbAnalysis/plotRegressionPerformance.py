from ROOT import *
import sys

ifile = TFile(sys.argv[1])
tree = ifile.Get("tree")

#nentries = tree.GetEntries()
nentries = 10000
print "Processing %i events" % nentries

hist = TH1F("hist","hist",100,-1,1)
hist_reg = TH1F("hist_reg","hist_reg",100,-1,1)

ofile = TFile("hists.root","RECREATE")

for ientry in range(nentries):
    if (ientry % 10000 == 0): print "Processeing entry %i" % ientry
    tree.GetEntry(ientry)
    genjets = []
    for i in range(tree.nGenJet):
        genjet = TLorentzVector()
        genjet.SetPtEtaPhiM(tree.GenJet_wNuPt[i],tree.GenJet_wNuEta[i],tree.GenJet_wNuPhi[i],tree.GenJet_wNuM[i])
        genjets.append(genjet)
    #for i in range(tree.nJet):
    jet1 = TLorentzVector()
    if tree.nJet < 2: continue
    ind1 = 0
    #ind1 = tree.hJetInd1
    jet1.SetPtEtaPhiM(tree.Jet_pt[ind1],tree.Jet_eta[ind1],tree.Jet_phi[ind1],tree.Jet_mass[ind1])
    jet2 = TLorentzVector()
    ind2 = 0
    #ind2 = tree.hJetInd2
    jet2.SetPtEtaPhiM(tree.Jet_pt[ind2],tree.Jet_eta[ind2],tree.Jet_phi[ind2],tree.Jet_mass[ind2])
    
    minDR1 = 999
    genmatched1 = TLorentzVector()
    minDR2 = 999
    genmatched2 = TLorentzVector()
 
    for genjet in genjets:
        dR = jet1.DeltaR(genjet)
        if (dR < minDR1 and dR <= 0.5):
            genmatched1 = genjet
            minDR1 = dR
    if (minDR1 < 999):
        genjets.remove(genmatched1)
    for genjet in genjets:
        dR = jet2.DeltaR(genjet)
        if (dR < minDR2 and dR <= 0.5):
            genmatched2 = genjet
            minDR2 = dR

    #print jet1.Pt(),genmatched1.Pt(),( ( jet1.Pt() - genmatched1.Pt() )  / genmatched1.Pt() )
    if (genmatched1.Pt() > 0):
        hist.Fill( ( jet1.Pt() - genmatched1.Pt() )  / genmatched1.Pt() )
        hist_reg.Fill( ( tree.Jet_pt_reg[ind1] - genmatched1.Pt() )  / genmatched1.Pt() )
    if (genmatched2.Pt() > 0):
        hist.Fill( ( jet2.Pt() - genmatched2.Pt() )  / genmatched2.Pt() )
        hist_reg.Fill( ( tree.Jet_pt_reg[ind2] - genmatched2.Pt() )  / genmatched2.Pt() )
        
ofile.cd()
hist.Write()
hist_reg.Write()
ofile.Close()
