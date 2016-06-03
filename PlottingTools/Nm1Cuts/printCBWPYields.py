## Author: Stephane Cooperstein
##
## Print weighted signal and background yields for a range of cut-based working points

import ROOT
import sys

ifile = ROOT.TFile(sys.argv[1])
tree = ifile.Get("tree")

wps = ["Analysis3Em3", "Analysis3Em3p5","Analysis4Em3","Analysis5Em3","Analysis6Em3","Analysis8Em3","Analysis1Em2","Analysis1p5Em2","Analysis2Em2","Analysis3Em2","Analysis4Em2","Analysis5Em2"]

#weight_string = "weight*weight_PU*bTagWeight"
weight_string = "weight*weight_PU*bTagWeight*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)*CS_SF"

for wp in wps:
    print wp
    hsig = ROOT.TH1F("hsig","hsig",60,90,150)
    hbkg = ROOT.TH1F("hbkg","hbkg",60,90,150)
    tree.Draw("H_mass>>hsig","(Pass_%s==1&&sampleIndex==-12501)*(%s)" % (wp,weight_string) )
    tree.Draw("H_mass>>hbkg","(Pass_%s==1&&sampleIndex>0&&sampleIndex!=12&&sampleIndex!=120)*(%s)" % (wp,weight_string) )
    print "Sig Yield = %f" % hsig.Integral()
    print "Bkg Yield = %f" % hbkg.Integral()

