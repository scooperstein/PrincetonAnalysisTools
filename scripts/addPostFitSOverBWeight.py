from ROOT import *
import sys
import numpy

## Usage: python addPostFitSOverBWeight.py [inputntuple] [outputntuple] [mlfit.root]
## Adds a branch, "sb_weight", with the per-event S/S+B weight for
## the events corresponding bin in the SR BDT Score distribution
## To be applied to all MC and data.


ifile = TFile.Open(sys.argv[1])
ofile = TFile.Open(sys.argv[2],"RECREATE")
tree = ifile.Get("tree")

# input for calculating S / S+B weight

fitfile = TFile.Open(sys.argv[3])
#hist_postfit_sig_bn = fitfile.Get("shapes_fit_s/WmnHighPt/total_signal")
hist_postfit_sig_bn = fitfile.Get("shapes_prefit/WmnHighPt/total_signal")
hist_postfit_bkg_bn = fitfile.Get("shapes_fit_s/WmnHighPt/total_background")

fitfile_vv = TFile.Open(sys.argv[4])
hist_postfit_sig_vv_bn = fitfile_vv.Get("shapes_fit_s/WmnHighPt/total_signal")
hist_postfit_bkg_vv_bn = fitfile_vv.Get("shapes_fit_s/WmnHighPt/total_background")

xmin = 0
xmax = 15 

#nentries = 1000
nentries = tree.GetEntries()

# WHbb has variable binning, need to manually input the binning to use for the BDT shape
# Important point is that the BDT histograms hist_postfit_sig and hist_postfit_bkg have the same BDT binning as the nominal analysis.

import numpy
binedges = [0.3,0.4,0.5,0.6,0.65,0.70,0.75,0.80,0.85,0.90,1.0]
binedges_vv = [0.2,0.3,0.4,0.5,0.6,0.65,0.70,0.75,0.80,0.85,0.90,1.0]
#binedges = [-1.,    -0.886, -0.746, -0.606, -0.466, -0.326, -0.186, -0.046,  0.094,  0.234,  0.374,  0.514,  0.654,  0.794,  0.934,  1.   ] 
binedge_array = numpy.zeros(len(binedges),dtype=float)
binedge_array_vv = numpy.zeros(len(binedges_vv),dtype=float)
for i in range(len(binedges)):
    binedge_array[i] = binedges[i]
for i in range(len(binedges_vv)):
    binedge_array_vv[i] = binedges_vv[i]

#hist_postfit_sig = TH1F("hist_postfit_sig","hist_postfit_sig",300,xmin,xmax)
#hist_postfit_bkg = TH1F("hist_postfit_bkg","hist_postfit_bkg",300,xmin,xmax)

hist_postfit_sig = TH1F("hist_postfit_sig","hist_postfit_sig",len(binedges)-1,binedge_array)
hist_postfit_bkg = TH1F("hist_postfit_bkg","hist_postfit_bkg",len(binedges)-1,binedge_array)
hist_postfit_sig_vv = TH1F("hist_postfit_sig_vv","hist_postfit_sig_vv",len(binedges_vv)-1,binedge_array_vv)
hist_postfit_bkg_vv = TH1F("hist_postfit_bkg_vv","hist_postfit_bkg_vv",len(binedges_vv)-1,binedge_array_vv)

# combine gives post-fit histograms with bin # as x-axis, need to correct first
for i in range(1,hist_postfit_sig.GetNbinsX()+1):
    hist_postfit_sig.SetBinContent(i,hist_postfit_sig_bn.GetBinContent(i))
    hist_postfit_bkg.SetBinContent(i,hist_postfit_bkg_bn.GetBinContent(i))
for i in range(1,hist_postfit_sig_vv.GetNbinsX()+1):
    hist_postfit_sig_vv.SetBinContent(i,hist_postfit_sig_vv_bn.GetBinContent(i))
    hist_postfit_bkg_vv.SetBinContent(i,hist_postfit_bkg_vv_bn.GetBinContent(i))


## NB: You may find it easier to save in a separate output file the hist_postfit_sig and hist_postfit_bkg histograms
## corresponding to the BDT hist in the SR for signal and summed background so that you can just 
## that in here instead of all the code above. i.e. something like
## ifile = TFile("hists.root")
## hist_postfit_sig = ifile.Get("hist_postfit_sig")
## hist_postfit_bkg = ifile.Get("hist_postfit_bkg")

ofile.cd()
otree = tree.CloneTree(0)
sb_weight = numpy.zeros(1,dtype=float)
sb_weight_vv = numpy.zeros(1,dtype=float)
sb_weight_comb = numpy.zeros(1,dtype=float)
otree.Branch("sb_weight",sb_weight,"sb_weight/D")
otree.Branch("sb_weight_vv",sb_weight_vv,"sb_weight_vv/D")
otree.Branch("sb_weight_comb",sb_weight_comb,"sb_weight_comb/D")

print "total entries: %i " % nentries
for ientry in range(nentries):
    tree.GetEntry(ientry)
    val = tree.CMS_vhbb_BDT_Wln_13TeV # branch for VH BDT
    val_vv = tree.CMS_vvbb_BDT_Wln_13TeV # branch for VV BDT
    val_bin = hist_postfit_sig.FindBin(val)
    val_bin_vv = hist_postfit_sig_vv.FindBin(val_vv)
    if val_bin == 0:
        val_bin = 1 # underflow
    if val_bin_vv == 0:
        val_bin_vv = 1 # underflow
    s = hist_postfit_sig.GetBinContent(val_bin)
    b = hist_postfit_bkg.GetBinContent(val_bin)
    s_vv = hist_postfit_sig_vv.GetBinContent(val_bin_vv)
    b_vv = hist_postfit_bkg_vv.GetBinContent(val_bin_vv)
    
    #if ((s+b)>0):
    #    sb_weight[0] = s / (s+b)
    #else:
    #    sb_weight[0] = 0.
 
    sb_weight[0] = 0.
    if (b > 0):
        #sb_weight[0] = (s + s_vv) / (b)
        #sb_weight[0] = (s + s_vv) / (b + b_vv)
        #if val > val_vv:
        #    sb_weight[0] = s / b
        #else:
        #    sb_weight[0] = s_vv / b_vv
        sb_weight_comb[0] = max( (s/b), (s_vv/b_vv) )
        sb_weight_vv[0] =  s_vv / b_vv
        sb_weight[0] =  s / (s + b)
 
    #print val,sb_weight[0]
    otree.Fill()
    if (ientry % 1000 == 0 or sb_weight[0] > 10.): 
        print "processing entry: %i" % ientry
        print s,s_vv,b,b_vv
        print val,val_vv,sb_weight[0]

ofile.cd()
otree.Write()

## plot weights for validity checks
#canv = TCanvas("canv","canv")
#hTot = hist_postfit_sig.Clone()
#hTot.Add(hist_postfit_bkg)
#hist_postfit_sig.Divide(hTot) # S / S+B
#hist_postfit_sig.Draw("hist")
#canv.SaveAs("sb_weights.pdf")
#canv.SetLogy(True)
#canv.SaveAs("sb_weights_log.pdf")


