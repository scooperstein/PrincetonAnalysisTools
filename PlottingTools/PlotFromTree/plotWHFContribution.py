import ROOT
import sys
from math import sqrt

ROOT.gROOT.SetBatch(True)

ifile = ROOT.TFile(sys.argv[1], "read")

tree = ifile.Get("tree")

nbins = 100
BDTmin = -1.0
BDTmax = 1.0
BDTname = "BDT_13TeV_H125Sig_0b1b2bWjetsTTbarBkg_Mjj"

binsize = (BDTmax - BDTmin)/nbins

#h1 = ROOT.TH1F("h1","h1",nbins,BDTmin, BDTmax)
graph = ROOT.TGraphErrors()

for i in range(nbins):
    h2 = ROOT.TH1F("h2","h2",nbins,BDTmin,BDTmax)
    h3 = ROOT.TH1F("h3","h3",nbins,BDTmin,BDTmax)
    cutVal = BDTmin + i*binsize
    tree.Draw("%s>>h2" % BDTname,"(%s > %f && sampleIndex>100&&(sampleIndex%%10==1 || sampleIndex%%10==2))*weight" % (BDTname,cutVal))
    tree.Draw("%s>>h3" % BDTname,"(%s > %f && sampleIndex>0)*weight" % (BDTname,cutVal))
    if (h2.Integral() != 0 and h3.Integral() != 0):
        graph.SetPoint(i, cutVal, h2.Integral()/h3.Integral() )
        graph.SetPointError(i, sqrt((1./h2.Integral()) + (1./h3.Integral())) * (h2.Integral()/h3.Integral()) )
    else:
        graph.SetPoint(i, cutVal, 0.)
        graph.SetPointError(i, 0.)

canv = ROOT.TCanvas("canv","canv")
graph.GetXaxis().SetTitle("BDT cut value")
graph.GetYaxis().SetTitle("WHF/BKG")
graph.Draw("APE*")
canv.SaveAs("WHFContribution.png")

ofile = ROOT.TFile("WHFContribution.root", "RECREATE")
ofile.cd()
graph.Write("WHFContribution")
ofile.Close()
    
