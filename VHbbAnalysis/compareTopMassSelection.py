from ROOT import *
import sys

ifilename = sys.argv[1]
ifile = TFile(ifilename,"r")
tree = ifile.Get("tree")

hist = TH1F("hist","hist",100,0.,5.)

#nentries = tree.GetEntries()
nentries = 5000
print "total entries: %i" % nentries

den = 0
num = 0
for i in range(nentries):
    if (i % 1000 == 0):
        print "processing entry: %i" %i
    tree.GetEntry(i)
    if tree.nGenLep != 1: continue
    if tree.nGenBQuarkFromTop != 2: continue
    
    den += 1
    gj1 = TLorentzVector()
    gj1.SetPtEtaPhiM(tree.GenBQuarkFromTop_pt[0],tree.GenBQuarkFromTop_eta[0],tree.GenBQuarkFromTop_phi[0],tree.GenBQuarkFromTop_mass[0])
    gj2 = TLorentzVector()
    gj2.SetPtEtaPhiM(tree.GenBQuarkFromTop_pt[1],tree.GenBQuarkFromTop_eta[1],tree.GenBQuarkFromTop_phi[1],tree.GenBQuarkFromTop_mass[1])
   
    gl = TLorentzVector()
    gl.SetPtEtaPhiM(tree.GenLep_pt[0],tree.GenLep_eta[0],tree.GenLep_phi[0],tree.GenLep_mass[0])
     
    dR1 = gl.DeltaR(gj1)
    dR2 = gl.DeltaR(gj2)

    genjet = gj1
    if (dR2 < dR1):
        genjet = gj2

    lep = TLorentzVector()
    lep.SetPtEtaPhiM(tree.selLeptons_pt[0],tree.selLeptons_eta[0],tree.selLeptons_phi[0],tree.selLeptons_mass[0])
 
    minDR = 999
    closestJet = TLorentzVector()
    for i in range(tree.nJet):
        if tree.Jet_pt[i] < 30 or tree.Jet_btagCSV[i] < 0.5: continue
        jet = TLorentzVector()
        jet.SetPtEtaPhiM(tree.Jet_pt[i],tree.Jet_eta[i],tree.Jet_phi[i],tree.Jet_mass[i])
       
        dR = jet.DeltaR(lep)
        if (dR < minDR):
            minDR = dR
            closestJet = jet

    jet1 = TLorentzVector()
    jet1.SetPtEtaPhiM(tree.Jet_pt[tree.hJetInd1],tree.Jet_eta[tree.hJetInd1],tree.Jet_phi[tree.hJetInd1],tree.Jet_mass[tree.hJetInd1])
    jet2 = TLorentzVector()
    jet2.SetPtEtaPhiM(tree.Jet_pt[tree.hJetInd2],tree.Jet_eta[tree.hJetInd2],tree.Jet_phi[tree.hJetInd2],tree.Jet_mass[tree.hJetInd2])
    
    #closestJet = jet1

    jetpt = TLorentzVector()
    jetpt.SetPtEtaPhiM(tree.Jet_pt[0],tree.Jet_eta[0],tree.Jet_phi[0],tree.Jet_mass[0])
   
    closestJet = jetpt 

    dr = closestJet.DeltaR(genjet)
    hist.Fill(dr)
    if (dr < 0.5):
        num += 1

ofile = TFile("out.root","RECREATE")
ofile.cd()
hist.Write()
print num,den
print ((1.0*num)/den)

