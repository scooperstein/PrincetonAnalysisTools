from ROOT import *
import sys
import numpy

## Usage: python addSOverBWeight.py [inputntuple] [outputntuple]
## Adds a branch, "sb_weight", with the per-event S/S+B weight for
## the events corresponding bin in the SR BDT Score distribution
## To be applied to all MC and data.


#ROOT.gSystem.SetBatch(True)

ifile = TFile(sys.argv[1])
tree = ifile.Get("tree")
ofile = TFile(sys.argv[2],"RECREATE")
#ofile = TFile("ofile.root","RECREATE")

# input for calculating S / S+B weight
tree_weight = TChain("tree")
tree_weight.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_mc.root")
tree_weight.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_signal.root")
tree_weight.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_ttpowheg.root")

## copy all the histograms and stuff do..


sel = "selLeptons_relIso_0<0.06&&H_mass>0&&H_mass<255&&cutFlow>=10&&V_pt>100&&((Vtype==2&&isWmunu)||(Vtype==3&&isWenu))" #selection on top of input tree that defines signal region
var = "CMS_vhbb_BDT_Wln_13TeV" # variable to use to categorize by S/B
#var = "BDT_V24_Jan31_noMbb_400_5_VVSig"
#varMap = {}
#varMap[var] = tree.CMS_vhbb_BDT_Wln_13TeV
xmin = -1.0
xmax = 1.0 

sig_sel = "(sampleIndex==-12500||sampleIndex==-12501)||(sampleIndex>=3500&&sampleIndex<=3702&&Sum$(abs(GenWZQuark_pdgId)==5)>=2)" # additional requirement to define signal
bkg_sel = "" # additional requirement to define background

#nentries = 1000
nentries = tree.GetEntries()

# WHbb has variable binning due to re-binner, need to manually input the binning to use for the BDT shape
# Important point is that the BDT histograms hSig and hBkg have the same BDT binning as the nominal analysis.

import numpy
binedges = [-1.,    -0.886, -0.746, -0.606, -0.466, -0.326, -0.186, -0.046,  0.094,  0.234,  0.374,  0.514,  0.654,  0.794,  0.934,  1.   ] 
binedge_array = numpy.zeros(len(binedges),dtype=float)
for i in range(len(binedges)):
    binedge_array[i] = binedges[i]

#hSig = TH1F("hSig","hSig",300,xmin,xmax)
#hBkg = TH1F("hBkg","hBkg",300,xmin,xmax)

hSig = TH1F("hSig","hSig",len(binedges)-1,binedge_array)
hBkg = TH1F("hBkg","hBkg",len(binedges)-1,binedge_array)

tree_weight.Draw("%s>>hSig" %var,"((%s)&&(%s))*weight" % (sel,sig_sel) )
tree_weight.Draw("%s>>hBkg" %var,"((%s)&&!(%s))*weight" % (sel,sig_sel) )

## NB: You may find it easier to save in a separate output file the hSig and hBkg histograms
## corresponding to the BDT hist in the SR for signal and summed background so that you can just 
## that in here instead of all the code above. i.e. something like
## ifile = TFile("hists.root")
## hSig = ifile.Get("hSig")
## hBkg = ifile.Get("hBkg")


otree = tree.CloneTree(0)
sb_weight = numpy.zeros(1,dtype=float)
otree.Branch("sb_weight3",sb_weight,"sb_weight3/D")

print "total entries: %i " % nentries
for ientry in range(nentries):
    tree.GetEntry(ientry)
    #val = varMap[var]
    val = tree.CMS_vhbb_BDT_Wln_13TeV # remember to change bdtname here as well
    #val = tree.BDT_V24_Jan31_noMbb_400_5_VVSig
    val_bin = hSig.FindBin(val)
    s = hSig.GetBinContent(val_bin)
    b = hBkg.GetBinContent(val_bin)
    sb_weight[0] = s / (s+b)
    #print val,sb_weight[0]
    otree.Fill()
    if (ientry % 10000 == 0): 
        print "processing entry: %i" % ientry
        print val,sb_weight[0]

otree.Write()

# plot weights for validity checks
canv = TCanvas("canv","canv")
hTot = hSig.Clone()
hTot.Add(hBkg)
hSig.Divide(hTot) # S / S+B
hSig.Draw("hist")
canv.SaveAs("sb_weights.pdf")
canv.SetLogy(True)
canv.SaveAs("sb_weights_log.pdf")


