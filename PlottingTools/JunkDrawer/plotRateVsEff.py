#
#   This program takes a ROOT file with a
#   small tree containing efficiencies and
#   rates at various thresholds in a few
#   variables.
#   
#   The objective of the program is to 
#   plot the phase space of viable working
#   points as a function of variable
#   thresholds.
#
#   Written by Chris Palmer (Princeton)
#

import ROOT
import numpy

filename="testing_phicut_new.root"
tfile=ROOT.TFile.Open(filename)
tree=tfile.Get("newtree")

j2PtCuts=numpy.arange(10,101,10)
j2PtCuts[0]=-10

hists=[]

can=ROOT.TCanvas("can","",2800,1200)
#can.cd()
#can.Divide(1,3)
#variables=["eff","effhpt100","rate"]
#title=["Efficiency of W#rightarrowl#nu Hbb","Efficiency of W#rightarrowl#nu Hbb (P_{T,H}>100 GeV)","Rate"]
#maxs=[0.016,0.008,1200]
#can.Divide(2,2)
#variables=["eff","effhpt100","effhpt200","rate"]
#title=["Efficiency of W#rightarrowl#nu Hbb","Efficiency of W#rightarrowl#nu Hbb (P_{T,H}>100 GeV)","Efficiency of W#rightarrowl#nu Hbb P_{T,H}>200 GeV","New Rate"]
#maxs=[0.076,0.02,0.002,160000]
can.Divide(2,1)
variables=["effhpt100","rate"]
title=["Efficiency of W#rightarrowl#nu Hbb, P_{T,H}>100 GeV","New Rate, MinBias, Hz"]
maxs=[0.00414,1200]
hists={}
for j2PtCut in j2PtCuts:
    for ivar in xrange(len(variables)):
        pad=can.cd(ivar+1)
        hists["temphist"+str(variables[ivar])]=ROOT.TH2F("temphist"+str(variables[ivar]),str(title[ivar]),9,40,130,10,10,35)
        #hists["temphist"+str(variables[ivar])]=ROOT.TH2F("temphist"+str(variables[ivar]),str(title[ivar])+", Jet2 P_{T}="+str(j2PtCut)+"GeV",9,40,130,10,10,35)
        hists["temphist"+str(variables[ivar])].GetXaxis().SetTitle("Lead Jet P_{T} (GeV)")
        hists["temphist"+str(variables[ivar])].GetXaxis().SetTitleSize(0.05)
        hists["temphist"+str(variables[ivar])].GetYaxis().SetTitleSize(0.05)
        hists["temphist"+str(variables[ivar])].GetXaxis().SetTitleOffset(0.9)
        hists["temphist"+str(variables[ivar])].GetYaxis().SetTitleOffset(0.9)
        hists["temphist"+str(variables[ivar])].GetYaxis().SetTitle("EG Seed P_{T} (GeV)")
        hists["temphist"+str(variables[ivar])].GetZaxis().SetRangeUser(0.0001,maxs[ivar])
        #pad.SetLogz()
        print "filling hist"
        cutstring="(jetPt2=="+str(j2PtCut)+"&&rate<1200&&effhpt100>0.00103)*"+str(variables[ivar])
        #cutstring="(jetPt2=="+str(j2PtCut)+")*"+str(variables[ivar])
        #cutstring="(jetPt2=="+str(j2PtCut)+"&&rate<1200)*"+str(variables[ivar])
        #cutstring="(jetPt2=="+str(j2PtCut)+"&&rate<1200&&(eff>0.005))*"+str(variables[ivar])
        print cutstring
        tree.Draw("elPt:jetPt1>>temphist"+str(variables[ivar]),cutstring,"colz")
        print "filled hist"
        
    can.Update()
    #raw_input()a
    can.SaveAs("effAndRate_"+str(j2PtCut)+".png")
