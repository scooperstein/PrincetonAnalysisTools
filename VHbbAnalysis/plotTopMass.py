from ROOT import *
import sys
from math import pow

ifile = TFile(sys.argv[1],"r")
tree = ifile.Get("tree")


nentries = tree.GetEntries()
print "total entries: %i" % nentries

hres = TH1F("hres","hres",50,-1,1)

for i in range(nentries):
    if (i%1000 == 0):
        print "processing entry: %i" % i
    tree.GetEntry(i)

    if tree.GenTop_decayMode[0] != 0: continue
    #if tree.nGenLep != 0: continue
    #gl = TLorentzVector()
    #gl.SetPtEtaPhiM(tree.GenLep_pt[0]
 
    topmass = tree.Top1_mass_fromLepton_regPT_w4MET
    topmass_gen = tree.GenTop_mass[0]

    mW = 80.38
    Lepton = TLorentzVector()
    Lepton.SetPtEtaPhiM(tree.selLeptons_pt[0],tree.selLeptons_eta[0],tree.selLeptons_phi[0],tree.selLeptons_mass[0])
    MET = TLorentzVector()
    MET.SetPtEtaPhiM(tree.met_pt,0.,tree.met_phi,0.)

    MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py()) 
    mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py()
    a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz())
    a2 = pow(a,2)
    b  = (pow(Lepton.Energy(),2.)*(MisET2) - pow(mu,2.))/(pow(Lepton.Energy(),2) - pow(Lepton.Pz(),2))

    mT = 173.21
    topmass_gen = mT

    #if (a2-b > 0): continue # skip real solutions
    #if (a2-b < 0): continue # skip complex solutions

    #hres.Fill((topmass-topmass_gen)/topmass_gen)
    hres.Fill((topmass-topmass_gen)/topmass_gen)



#hres = TH1F("hres","hres",50,-1,1)
canv = TCanvas("canv","canv")
#tree.Draw("(Top1_mass_fromLepton_regPT_w4MET-GenTop_mass[0])/GenTop_mass[0]>>hres","GenTop_decayMode[0]==0")
hres.SetTitle("Reco. Top Mass Resolution")
hres.GetXaxis().SetTitle("(mTop(reco)-mTop(gen))/mTop(gen)")
hres.Draw()
#raw_input()
canv.SaveAs("hres.C")

