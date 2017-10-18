import os
import ROOT
import sys
from math import sqrt,pow

ROOT.gStyle.SetOptStat(0)

bins = "test"
#bins = ["Mjj","PtJJ","JJPtBal","detajj","dphijj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","leppt","lepeta","lepj1dphi","met","vpt","vmt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","dphijj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","leppt","lepeta","lepj1dphi","met","vpt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","leppt","lepeta","lepj1dphi","met","vpt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","leppt","lepeta","met","vpt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","met","vpt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","vpt","veta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","JJPtBal","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","vpt","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","vpt","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","rpt","yW","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","rpt","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1pt","jet2pt","jet1eta","jet2eta","jet1qgl","jet2qgl","lepeta","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1pt","jet1eta","jet1qgl","jet2qgl","lepeta","rpt","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1eta","jet1qgl","jet2qgl","lepeta","rpt","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1eta","jet1qgl","jet2qgl","rpt","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1qgl","jet2qgl","rpt","zW"]
#bins = ["Mjj","PtJJ","detajj","jet1qgl","jet2qgl","zW","leppt"]
#bins = ["Mjj","zW","detajj","jet1qgl","jet2qgl","PtJJ"]
pwd = os.getcwd()

#ref = 76.378243 # 7 var
#ref = 7.924923e-01 # 7 var
#ref = 75.610191 # 6 var
#ref = 0.792492 # 6 var
#ref = 80.378270 # 11 var
#ref = 0.804015 # 11 var
ref = 81.653499 # 21 var
#ref = 0.806454 # 21 var

hist = ROOT.TH1F("hist","hist",len(bins),0,len(bins))
histsig = ROOT.TH1F("histsig","histsig",len(bins),0,len(bins))

for i in range(1):
    b = bins[i]
    #textfile = fopen("%s/stdout_%s_Nm1","r")
    ifile = ROOT.TFile(sys.argv[1])
    #ifile = ROOT.TFile("%s/TMVA_13TeV_V25_EWK_June14_WJetsIncMad_21VarNm1%s.root" % (b,b))
    ##ifile = ROOT.TFile("../TMVA_13TeV_V25_EWK_June14_WJetsIncMad_21VarNm1.root")
    hroc = ifile.Get("Method_BDT/EWKWJets/MVA_EWKWJets_rejBvsS")
    roc = hroc.Integral() / hroc.GetNbinsX()

    # calc integrated S/sqrt(S+B)
    hs = ifile.Get("Method_BDT/EWKWJets/MVA_EWKWJets_S")
    hb = ifile.Get("Method_BDT/EWKWJets/MVA_EWKWJets_B")

    hs.Scale(93376.42/hs.Integral())
    hb.Scale(5776112.49/hb.Integral())

    sig = 0.
    for j in range(hs.GetNbinsX(),0,-1):
        s_bin = hs.GetBinContent(j)
        b_bin = hb.GetBinContent(j)
        #print j,sqrt(pow(s_bin,2) / (s_bin + b_bin))
        sig += pow(s_bin,2) / (s_bin + b_bin)
    sig = sqrt(sig)

 #   print maxSig,maxSigCutVal
    ifile.Close()
    #print "%s: %f" % (b,roc/ref)
    print "%s: %f : %f" % (b,roc,sig)
    hist.GetXaxis().SetBinLabel(i+1,b)
    histsig.GetXaxis().SetBinLabel(i+1,b)
    hist.Fill(i,roc)
    histsig.Fill(i,sig)

canv = ROOT.TCanvas("canv","canv")
#hist.GetYaxis().SetTitleOffset(1.25)
#hist.SetMarkerStyle(29)
#hist.GetYaxis().SetTitle("N-1 ROC Integral")
#hist.SetTitle("")
#hist.SetMaximum(1.0025*ref)
#hist.Draw("hist pl")
histsig.GetYaxis().SetTitleOffset(1.25)
histsig.SetMarkerStyle(29)
histsig.GetYaxis().SetTitle("N-1 Integrated S / sqrt(S+B)")
histsig.SetTitle("")
histsig.SetMaximum(1.02*ref)
histsig.Draw("histsig pl")
line = ROOT.TLine(0,ref,len(bins),ref)
line.SetLineColor(ROOT.kRed)
line.Draw("same");

#axis = ROOT.TGaxis(ROOT.gPad.GetUxmin(),ROOT.gPad.GetUymin(),ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymax(),0,1.1*histsig.GetMaximum(),510,"+L")
#axis.SetLineColor(ROOT.kRed)
#axis.SetTextColor(ROOT.kRed)
#axis.Draw()

#canv.SaveAs("Nm1ROC_21.pdf")

  
