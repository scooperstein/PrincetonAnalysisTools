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


plotlist={}

#plotlist["hltElePtsig"]       =["hltElePt",   "(sampleIndex<0)*weight",40,0,80]
#plotlist["hltElePtbkg"]       =["hltElePt",   "(sampleIndex>0)*weight",40,0,80]
#plotlist["phijjsig"]       =["phijj",   "(sampleIndex<0)*weight",30,0,3.15]
#plotlist["phijjbkg"]       =["phijj",   "(sampleIndex>0)*weight",30,0,3.15]
#plotlist["pass_eleHLTsig"]       =["pass_eleHLT",   "(sampleIndex<0)*weight",2,0,2]
#plotlist["pass_eleHLTbkg"]       =["pass_eleHLT",   "(sampleIndex>0)*weight",2,0,2]
#plotlist["pass_elePlusMetHLTsig"]=["pass_elePlusMetHLT",   "(sampleIndex<0)*weight",2,0,2]
#plotlist["pass_elePlusMetHLTbkg"]=["pass_elePlusMetHLT",   "(sampleIndex>0)*weight",2,0,2]
#plotlist["elIsoPtsig"]  =["elIsoPt",   "(sampleIndex<0)*weight",70,0,70]
#plotlist["elIsoPtbkg"]  =["elIsoPt",   "(sampleIndex>0)*weight",70,0,70]
#plotlist["elPtsig"]     =["elPt",   "(sampleIndex<0)*weight",70,0,70]
#plotlist["elPtbkg"]     =["elPt",   "(sampleIndex>0)*weight",70,0,70]
#plotlist["jet1Ptsig"]   =["jetPt1", "(sampleIndex<0)*weight",60,0,300]
#plotlist["jet1Ptbkg"]   =["jetPt1", "(sampleIndex>0)*weight",60,0,300]
#plotlist["jet2Ptsig"]   =["jetPt2", "(sampleIndex<0)*weight",50,-3,100]
#plotlist["jet2Ptbkg"]   =["jetPt2", "(sampleIndex>0)*weight",50,-3,100]
#plotlist["notPassBaselinesig"]   =["notPassBaseline", "(sampleIndex<0)*weight",3,-1,2]
#plotlist["notPassBaselinebkg"]   =["notPassBaseline", "(sampleIndex>0)*weight",3,-1,2]
#plotlist["passCurrentL1sig"]   =["passCurrentL1", "(sampleIndex<0)*weight",3,-1,2]
#plotlist["passCurrentL1bkg"]   =["passCurrentL1", "(sampleIndex>0)*weight",3,-1,2]

plotlist["jet1Pt"]             =["Jet_pt[0]","",100,0,400]
plotlist["jet2Pt"]             =["Jet_pt[1]","",100,0,400]

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
    leg.AddEntry(temp1,"Old Files")
    leg.AddEntry(temp2,"New Files")
    leg.Draw()

    can.Update()
    can.SaveAs(plot+".png")
    #raw_input()
