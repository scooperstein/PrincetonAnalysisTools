#   
#   This program takes two ROOT files containing
#   trees with the same name.  The variables in the 
#   trees are assumed to be named the same as well.
#   
#   plotlist contains the plotting information
#   for each desired histogram.
#   
#   The objective of the script is to validate 
#   variables in a new sample w.r.t. a previously
#   validated sample.
#
#   Written by Chris Palmer (Princeton)
#   

import ROOT
import sys

ROOT.gROOT.SetBatch(True)

filename1=sys.argv[1]
filename2=sys.argv[2]


tfile1=ROOT.TFile(filename1)
tfile2=ROOT.TFile(filename2)

treename="tree"

tree1=tfile1.Get(treename)
tree2=tfile2.Get(treename)

can=ROOT.TCanvas("can","",800,500)

can.cd()

presel = "(1==1)*sign(genWeight)"

plotlist={}

plotlist["JetPt"]              = ["Jet_pt",presel,100,0,200]
plotlist["JetEta"]             = ["Jet_eta",presel,50,-2.5,2.5]
plotlist["JetPhi"]             = ["Jet_phi",presel,50,-3.2,3.2]
plotlist["JetMass"]            = ["Jet_mass",presel,100,0,200]
plotlist["JetCSV"]             = ["Jet_btagCSV",presel,100,0,1.]
plotlist["JetPuId"]            = ["Jet_puId",presel,2,0,2]
plotlist["JetId"]              = ["Jet_id",presel,5,0,5]
plotlist["LepPt"]              = ["selLeptons_pt",presel,100,0,200]
plotlist["LepEta"]             = ["selLeptons_eta",presel,50,-2.5,2.5]
plotlist["LepPhi"]             = ["selLeptons_phi",presel,50,-3.2,3.2]
plotlist["LepMass"]            = ["selLeptons_mass",presel,100,0,200]
plotlist["LepPdgId"]           = ["abs(selLeptons_pdgId)",presel,25,0,25]
plotlist["LepId"]              = ["selLeptons_tightId",presel,3,0,3]
plotlist["MetPt"]              = ["met_pt",presel,100,0,400]
plotlist["MetPhi"]             = ["met_phi",presel,50,-3.2,3.2]
plotlist["PuppiMetPt"]         = ["metPuppi_pt",presel,100,0,400]
plotlist["HPt"]                = ["H_pt",presel,100,0,400]
plotlist["HEta"]               = ["H_eta",presel,50,-2.5,2.5]
plotlist["HPhi"]               = ["H_phi",presel,50,-3.2,3.2]
plotlist["HMass"]              = ["H_mass",presel,100,0,500]
plotlist["WPt"]                = ["V_pt",presel,100,0,400]
plotlist["NAddJet"]            = ["Sum$(Jet_pt>25&&Jet_puId==1&&abs(Jet_eta)<2.9)",presel,10,0,10]
plotlist["NAddLep"]            = ["Sum$(selLeptons_pt>15&&abs(selLeptons_eta)<2.5)",presel,10,0,10]
plotlist["HVDPhi"]             = ["HVdPhi",presel,50,-3.2,3.2]

plots=plotlist.keys()
plots.sort()

for plot in plots:
    print plot
    temp1=ROOT.TH1F("temp1",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    tree1.Draw(plotlist[plot][0]+">>temp1",plotlist[plot][1])
    temp2=ROOT.TH1F("temp2",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1])

    temp1.SetLineColor(ROOT.kRed)
    temp2.SetLineColor(ROOT.kBlue)

    temp1.Scale(1./temp1.Integral())
    temp2.Scale(1./temp2.Integral())

    maxval=max(temp1.GetMaximum(),temp2.GetMaximum())*1.2

    print "max1, max2,",temp1.GetMaximum(),temp2.GetMaximum()
    print "integral1, integral2,",temp1.Integral(),temp2.Integral()
    temp1.SetMaximum(maxval)
    temp1.Draw()
    temp2.Draw("sames")

    leg=ROOT.TLegend(0.7,0.7,1.,1.)
    leg.AddEntry(temp1,"V12")
    leg.AddEntry(temp2,"V13")
    leg.Draw()

    can.Update()
    can.SaveAs(plot+".png")
    #raw_input()
