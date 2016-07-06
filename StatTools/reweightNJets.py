##
## Derive estimate of data/MC scale factors for TT powheg only vs. data
## for the number of additional jets distribution.
##
## Author: Stephane Cooperstein
##

import ROOT
import sys

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

#path = sys.argv[1]
#tree = ROOT.TChain("tree")
#tree.Add("%s/*/*.root" % path)

ifile = ROOT.TFile(sys.argv[1])
tree = ifile.Get("tree")


nbins = 10
xlow = 0
xhigh = 10

sel = "Jet_btagCSV[hJetInd2]>0.46&&(H_mass<90||H_mass>150)&&controlSample==1&&Pass_nominal==1&&(isWenu||isWmunu)&&(Vtype==2||Vtype==3)"

sampleMap = {}
sampleMap["TT"] = [120]
sampleMap["s_Top"] = [16,17,20,21]
sampleMap["WH"] = [-12501]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4400,4500,4600,4700,4800,4900]
sampleMap["Wj1b"] = [2201,4401,4501,4601,4701,4801,4901]
sampleMap["Wj2b"] = [2202,4402,4502,4602,4702,4802,4902]
sampleMap["VVHF"] = [3501,3502,3601,3602,3701,3702]
sampleMap["VVLF"] = [3500,3600,3700]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
sampleMap["Zj0b"] = [2300,6100,6200,6300,6400]
sampleMap["Zj1b"] = [2301,6101,6201,6301,6401]
sampleMap["Zj2b"] = [2302,6102,6202,6302,6402]

tt = ROOT.TH1F("tt","tt",nbins,xlow,xhigh)

ifile_data = ROOT.TFile(sys.argv[2])
tree_data = ifile_data.Get("tree")
data = ROOT.TH1F("data","data",nbins,xlow,xhigh)
tree_data.Draw("nAddJet_f>>data","((%s)&&sampleIndex==0)*weight" % sel)
print "Data = ",data.Integral()
dataTot = data.Clone("dataTot")
ifile.cd()
tree.Draw("nAddJet_f>>tt","((%s)&&sampleIndex==%i)*weight" % (sel,sampleMap["TT"][0]))

nonTT = ROOT.TH1F("nonTT","nonTT",nbins,xlow,xhigh)

print "Data = ",data.Integral()
print "TT = ",tt.Integral()

for sample in sampleMap:
   if sample == "TT": continue
   cutString = "(%s)&&(" % sel
   for index in sampleMap[sample]:
       cutString += "(sampleIndex==%i)||" % index 
   cutString = cutString[0:len(cutString)-2]
   cutString += ")"
   print sample
   print cutString
   hist = ROOT.TH1F(sample,sample,nbins,xlow,xhigh)
   tree.Draw("nAddJet_f>>%s" % sample,"((%s))*weight" % cutString)
   data.Add(hist,-1)
   print hist.Integral()
   nonTT.Add(hist)

canv = ROOT.TCanvas("canv","canv")
data.SetTitle("")
data.GetXaxis().SetTitle("Number of Additional Jets")
data.Draw("ep")
tt.SetLineColor(ROOT.kBlue)
tt.Draw("hist same")
tt.Draw("ep same")

leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg.AddEntry(data, "2.32/fb 2015 Data - non-TT MC")
leg.AddEntry(tt, "TT Powheg MC")
leg.Draw("same")

canv.SaveAs("TT_vs_data_minus_nonTTMC.pdf")

nonTT.SetLineColor(ROOT.kRed)
stack = ROOT.THStack()
stack.Add(tt)
stack.Add(nonTT)

dataTot.SetTitle("")
dataTot.GetXaxis().SetTitle("Number of Additional Jets")
dataTot.Draw("ep")
stack.Draw("hist same")
stack.Draw("ep same")

leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg.AddEntry(data, "2.32/fb 2015 Data","ep")
leg.AddEntry(tt, "TT Powheg MC")
leg.AddEntry(nonTT, "All other MC")
leg.Draw("same")

canv.SaveAs("MC_vs_data_nAddJet.pdf")

data.Divide(tt)
for i in range(1,data.GetNbinsX()+1):
    print "NAddJet %i: %f" % (data.GetBinLowEdge(i), data.GetBinContent(i))

data.GetYaxis().SetTitle("Data - TTnonMC / TT")
data.Draw()
canv.SaveAs("data_minus_nonTTMC_over_TT.pdf") 
  
