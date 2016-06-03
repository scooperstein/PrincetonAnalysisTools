## Compare systematic shape variation for given shape systematic vs. nominal for BDT
##
## Author: Stephane Cooperstein
##

import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser("Make systematic shape comparison plot")
parser.add_argument('-i', '--inputfile', type=str, default="", help="The input root file w/ shapes")
#parser.add_argument('-s', '--systematic', type=str, default="", help="the name of the shape systematic")
#parser.add_argument('-v','--variable', type=str, default="CMS_vhbb_BDT_Wln_13TeV", help="the variable to plot")
args = parser.parse_args()

samples = ["WH","TT","Wj0b","Wj1b","Wj2b","s_Top"]
#samples = ["Wj1b"]
#systematics = ["CMS_vhbb_statWj1b_Comb_bin16_13TeV"]
systematics = ["CMS_vhbb_scale_j","CMS_vhbb_res_j","CMS_vhbb_bTagWeightHF","CMS_vhbb_bTagWeightLF","CMS_vhbb_bTagWeightLFStats1","CMS_vhbb_bTagWeightHFStats1","CMS_vhbb_bTagWeightLFStats2","CMS_vhbb_bTagWeightHFStats2","CMS_vhbb_btagWeightcErr1","CMS_vhbb_btagWeightcErr2","CMS_vhbb_TTModel_Wln_13TeV"]

ifile = ROOT.TFile(args.inputfile, "r")
canv = ROOT.TCanvas("canv","canv")
for systematic in systematics:
    for sample in samples:
        if (systematic == "CMS_vhbb_TTModel_Wln_13TeV" and sample != "TT"): continue
        hnominal = ifile.Get(sample)
        hUp = ifile.Get("%s_%sUp" % (sample,systematic))
        hDown = ifile.Get("%s_%sDown" % (sample,systematic))

        hUp.SetLineColor(ROOT.kRed)
        hDown.SetLineColor(ROOT.kBlue)
        hnominal.SetLineColor(ROOT.kBlack)
 
        hUp.SetTitle("Shape Systematic: %s, Sample: %s" % (systematic,sample))

        hUp.SetMaximum(1.2*max(hUp.GetMaximum(),hDown.GetMaximum(),hnominal.GetMaximum()))

        hUp.Draw("hist")
        hnominal.Draw("ep0 same")
        hDown.Draw("hist same")
        canv.Update()
        print sample,hnominal.Integral()
        leg = ROOT.TLegend(0.7,1.0,0.7,1.0)
        leg.AddEntry(hnominal,"nominal")
        leg.AddEntry(hUp, "up")
        leg.AddEntry(hDown, "down")
     
        leg.Draw("same")
        canv.Update()
        canv.SaveAs("%s_%s.png" % (systematic,sample))
        canv.Update()
        #canv.Close()
