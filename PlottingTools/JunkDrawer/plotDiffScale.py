#   
#   From a flat ROOT tree this program   
#   plots rate and efficiency as a function
#   of selected variables.  These rates
#   and efficiencies are plotted in a two
#   dimensional histogram as well as in one
#   dimensional histograms on the same 
#   canvas with two different scales listed.
#   
#   Written by Chris Palmer (Princeton)
#   Partially based on 
#   https://root.cern.ch/root/html/tutorials/hist/transpad.C.html
#

import ROOT
import math
import array

xlen=1400
ylen=1100


# FIXME this function crashes the first time
def transpad(myfile,name1,name2):
    c1 = ROOT.TCanvas("c1","transparent pad",200,10,xlen,ylen)
    pad1 = ROOT.TPad("pad1","",0,0,1,1)
    pad2 = ROOT.TPad("pad2","",0,0,1,1)
    pad2.SetFillStyle(4000) #will be transparent
    pad1.Draw()
    pad1.cd()

    h1 = myfile.Get(name1)
    h1.GetYaxis().SetTitle("Rate of MinBias (Hz)")
    h1.GetYaxis().SetTitleOffset(1.5)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetLineWidth(2)
    h1.SetTitle("")
    h2 = myfile.Get(name2)
    h1.Draw()
    pad1.Update() #this will force the generation of the "stats" box
    #TPaveStats *ps1 = (TPaveStats*)h1.GetListOfFunctions().FindObject("stats")
    #ps1.SetX1NDC(0.4) ps1.SetX2NDC(0.6)
    pad1.Modified()
    c1.cd()

    #compute the pad range with suitable margins
    ymin = 0
    ymax = h2.GetMaximum()*1.1
    dy = (ymax-ymin)/0.8 #10 per cent margins top and bottom
    xmin = h2.GetBinLowEdge(1)
    xmax = h2.GetBinLowEdge(h1.GetNbinsX()+1)
    dx = (xmax-xmin)/0.8 #10 per cent margins left and right
    pad2.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
    pad2.Draw()
    pad2.cd()
    h2.SetLineColor(ROOT.kRed+1)
    h2.SetLineWidth(2)
    h2.Draw("][sames")
    pad2.Update()
    #TPaveStats *ps2 = (TPaveStats*)h2.GetListOfFunctions().FindObject("stats")
    #ps2.SetX1NDC(0.65) ps2.SetX2NDC(0.85)
    #ps2.SetTextColor(ROOT.kRed)

    # draw axis on the right side of the pad
    axis = ROOT.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L")
    axis.SetLabelColor(ROOT.kRed)
    axis.SetTitle("Efficiency of WlnuHbb")
    axis.SetTitleOffset(1.4)
    axis.SetLabelColor(ROOT.kRed+1)
    axis.Draw("same")
    pad2.Update()

    c1.SaveAs(name1+"_"+name2+".png")


def RateVsEff(myfile,name1,name2):
    can=ROOT.TCanvas("can","can",xlen,ylen)
    can.cd()
    h1 = myfile.Get(name1)
    h2 = myfile.Get(name2)
    integrals1=[]
    integrals2=[]
    errors1=[]
    errors2=[]
    
    bins=int(h1.GetNbinsX())

    for bin in xrange(bins):
        integral1=h1.Integral(bin,bins)
        integral2=h2.Integral(bin,bins)

        integrals1.append(integral1)
        integrals2.append(integral2)
        
        if integral1 != 0:
            error1=math.sqrt(integral1/0.62)*0.62
        else:
            error1=0

        if integral2 != 0:
            error2=math.sqrt(integral2*200000)/200000
        else:
            error2=0

        errors1.append(error1)
        errors2.append(error2)
  
    grph=ROOT.TGraphErrors(bins,array.array('f',integrals1),array.array('f',integrals2))#,array.array('f',errors1),array.array('f',errors2))
    grph.SetTitle("")
    grph.GetXaxis().SetTitle("Rate of MinBias (Hz)")
    grph.GetYaxis().SetTitle("Efficiency of WlnuHbb")
    grph.GetYaxis().SetTitleOffset(1.5)
    grph.Draw("AL")

    can.SaveAs(name1+"effvsrate.png")

myfile=ROOT.TFile("output.root")

basenames=["dphiejet_cat0_","leadjetpt_cat0_","subjetpt_cat0_","abseisoeta_cat0_","eisoeta_cat0_","eisopt0_cat0_"]

for basename in basenames:
    name1=basename+"nugun"
    name2=basename+"WH125"

    RateVsEff(myfile,name1,name2)
    transpad(myfile,name1,name2)
