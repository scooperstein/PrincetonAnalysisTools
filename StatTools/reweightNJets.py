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

#ifile = ROOT.TFile(sys.argv[1])
#tree = ifile.Get("tree")
tree = ROOT.TChain("tree")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WH125_powheg/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/DYToLL_HT100to200/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/DYToLL_HT200to400/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/DYToLL_HT400to600/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/DYToLL_HT600toInf/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/DYToLL_madgraph/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow/TT_powheg/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/TToLeptons_s/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/TToLeptons_t/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/T_tW/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/Tbar_tW/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WBJets/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT100To200/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT1200To2500/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT200To400/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT2500ToInf/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT400To600/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT600To800/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets-HT800To1200/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets_BGenFilter/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WJets_madgraph/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WW/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/WZ/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/ZZ/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT100to200/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT200to300/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT300to500/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT500to700/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT700to1000/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT1000to1500/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT1500to2000/*.root")
tree.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v3/QCD_HT2000toInf/*.root")

nbins = 10
xlow = 0
xhigh = 10

import numpy
bins = numpy.zeros(3,dtype=float)
bins[0] = -0.5
bins[1] = 1.5
bins[2] = 9.5

#sel = "Jet_btagCSV[hJetInd2]>0.46&&(H_mass<90||H_mass>150)&&controlSample==1&&Pass_nominal==1&&(isWenu||isWmunu)&&(Vtype==2||Vtype==3)"
sel = "V_pt>100&&cutFlow>=5&&H_pt>100&&Jet_btagCSV[hJetInd1]>0.935&&Jet_btagCSV[hJetInd2]>0.46&&(H_mass<90||H_mass>150)&&(isWenu||isWmunu)&&(Vtype==2||Vtype==3)"

sampleMap = {}
sampleMap["TT"] = [120]
sampleMap["s_Top"] = [16,17,20,21]
sampleMap["WH"] = [-12501]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4100,4200,4300,4400,4500,4600,4700,4800,4900]
sampleMap["Wj1b"] = [2201,4101,4201,4301,4401,4501,4601,4701,4801,4901]
sampleMap["Wj2b"] = [2202,4102,4202,4302,4402,4502,4602,4702,4802,4902]
sampleMap["VVHF"] = [3501,3502,3601,3602,3701,3702]
sampleMap["VVLF"] = [3500,3600,3700]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
sampleMap["Zj0b"] = [2300,6100,6200,6300,6400]
sampleMap["Zj1b"] = [2301,6101,6201,6301,6401]
sampleMap["Zj2b"] = [2302,6102,6202,6302,6402]

tt = ROOT.TH1F("tt","tt",nbins,xlow,xhigh)
#tt = ROOT.TH1F("tt","tt",len(bins)-1,bins)

#ifile_data = ROOT.TFile(sys.argv[2])
#tree_data = ifile_data.Get("tree")
tree_data = ROOT.TChain("tree")
tree_data.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/Run2016BToG_Ele_v2/*.root")
tree_data.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29_cutFlow_v4/Run2016BToG_Mu_v2/*.root")
data = ROOT.TH1F("data","data",nbins,xlow,xhigh)
#data = ROOT.TH1F("data","data",len(bins)-1,bins)
tree_data.Draw("nAddJet_f>>data","((%s)&&sampleIndex==0&&run<=276811)" % sel)
print "Data = ",data.Integral()
dataTot = data.Clone("dataTot")
#ifile.cd()
tree.Draw("nAddJet_f>>tt","((%s)&&sampleIndex==%i)*(weight/(1.))" % (sel,sampleMap["TT"][0]))
tree.Draw("nAddJet_f>>tt","((%s)&&sampleIndex==%i)*(0.85/0.88)*(weight/(1.))" % (sel,sampleMap["TT"][0]))
#tree.Draw("nAddJet_f>>tt","((%s)&&sampleIndex==%i)*(weight/(20.)*(1./CS_SF))" % (sel,sampleMap["TT"][0]))

nonTT = ROOT.TH1F("nonTT","nonTT",nbins,xlow,xhigh)
#nonTT = ROOT.TH1F("nonTT","nonTT",len(bins)-1,bins)

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
   #hist = ROOT.TH1F(sample,sample,len(bins)-1,bins)
   tree.Draw("nAddJet_f>>%s" % sample,"((%s))*(weight/(1.))" % cutString)
   data.Add(hist,-1)
   print hist.Integral()
   print "0-addjet bin: ",hist.GetBinContent(1)," +/- ",hist.GetBinError(1)
   nonTT.Add(hist)

print "Total non-TT MC: ",nonTT.Integral()
print "non-TT MC 0-addjet bin: ",nonTT.GetBinContent(1)," +/- ",nonTT.GetBinError(1)

data.SetMaximum(1.3*max(data.GetMaximum(),tt.GetMaximum()))
data.SetMinimum(0.0)

canv = ROOT.TCanvas("canv","canv")
data.SetTitle("")
data.GetXaxis().SetTitle("Number of Additional Jets")
data.Draw("ep")
tt.SetLineColor(ROOT.kBlue)
tt.Draw("hist same")
tt.Draw("ep same")

leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg.AddEntry(data, "12.9/fb 2016 Data - non-TT MC")
#leg.AddEntry(data, "2.32/fb 2015 Data - non-TT MC")
leg.AddEntry(tt, "TT Powheg MC")
leg.Draw("same")

canv.SaveAs("TT_vs_data_minus_nonTTMC.pdf")

nonTT.SetLineColor(ROOT.kRed)
stack = ROOT.THStack()
stack.Add(tt)
stack.Add(nonTT)

dataTot.SetMaximum(1.3*max(dataTot.GetMaximum(),(tt.GetMaximum()+nonTT.GetMaximum())))
dataTot.SetMinimum(0.0)

dataTot.SetTitle("")
dataTot.GetXaxis().SetTitle("Number of Additional Jets")
dataTot.Draw("ep")
stack.Draw("hist same")
stack.Draw("ep same")

leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg.AddEntry(data, "12.9/fb 2016 Data","ep")
leg.AddEntry(tt, "TT Powheg MC")
leg.AddEntry(nonTT, "All other MC")
leg.Draw("same")

canv.SaveAs("MC_vs_data_nAddJet.pdf")

data.Divide(tt)
for i in range(1,data.GetNbinsX()+1):
    print "NAddJet %i: %f" % (data.GetBinLowEdge(i), data.GetBinContent(i))

#data.SetMaximum(1.3*data.GetMaximum())
#data.SetMinimum(0.7*data.GetMinimum())
data.SetMaximum(1.2)
data.SetMinimum(0.5)
data.GetYaxis().SetTitle("Data - TTnonMC / TT")
data.Draw()
canv.SaveAs("data_minus_nonTTMC_over_TT.pdf") 
  
