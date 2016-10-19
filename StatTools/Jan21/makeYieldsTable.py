import ROOT
from math import sqrt
from math import pow

ifile_ttWmn = ROOT.TFile("hists_ttWmn.root")
ifile_ttWen = ROOT.TFile("hists_ttWen.root")
ifile_whfWmn = ROOT.TFile("hists_whfWmn.root")
ifile_whfWen = ROOT.TFile("hists_whfWen.root")
ifile_wlfWmn = ROOT.TFile("hists_wlfWmn.root")
ifile_wlfWen = ROOT.TFile("hists_wlfWen.root")
ifile_wen = ROOT.TFile("hists_WenHighPt.root")
ifile_wmn = ROOT.TFile("hists_WmnHighPt.root")

mcTot_0 = 0.
mcTotEnt_0 = 0
mcTot_1 = 0.
mcTotEnt_1 = 0
mcTot_2 = 0.
mcTotEnt_2 = 0
mcTot_3 = 0.
mcTotEnt_3 = 0
mcTot_4 = 0.
mcTotEnt_4 = 0
mcTot_5 = 0.
mcTotEnt_5 = 0
#for sample in ["data_obs"]:
#for sample in ["Zj0b","Zj1b","Zj2b","Wj0b","Wj1b","Wj2b","TT","s_Top","VVHF","VVLF"]:
for sample in ["Zj0b","Zj1b","Zj2b","Wj0b","Wj1b","Wj2b","TT","s_Top","VVLF","VVHF","ZH","WH"]:
    h_ttWmn = ifile_ttWmn.Get("BDT_ttWmn_%s" % sample)
    mcTot_2 += h_ttWmn.Integral()
    mcTotEnt_2 += pow(h_ttWmn.Integral() / sqrt(h_ttWmn.GetEffectiveEntries()), 2)
    h_ttWen = ifile_ttWen.Get("BDT_ttWen_%s" % sample)
    mcTot_3 += h_ttWen.Integral()
    mcTotEnt_3 += pow(h_ttWen.Integral() / sqrt(h_ttWen.GetEffectiveEntries()), 2)
    h_whfWmn = ifile_whfWmn.Get("BDT_whfWmn_%s" % sample)
    mcTot_4 += h_whfWmn.Integral()
    mcTotEnt_4 += pow(h_whfWmn.Integral() / sqrt(h_whfWmn.GetEffectiveEntries()), 2)
    h_whfWen = ifile_whfWen.Get("BDT_whfWen_%s" % sample)
    mcTot_5 += h_whfWen.Integral()
    mcTotEnt_5 += pow(h_whfWen.Integral() / sqrt(h_whfWen.GetEffectiveEntries()), 2)
    h_wlfWmn = ifile_wlfWmn.Get("BDT_wlfWmn_%s" % sample)
    mcTot_0 += h_wlfWmn.Integral()
    mcTotEnt_0 += pow(h_wlfWmn.Integral() / sqrt(h_wlfWmn.GetEffectiveEntries()), 2)
    h_wlfWen = ifile_wlfWen.Get("BDT_wlfWen_%s" % sample) 
    mcTot_1 += h_wlfWen.Integral()
    mcTotEnt_1 += pow(h_wlfWen.Integral() / sqrt(h_wlfWen.GetEffectiveEntries()), 2)
    #print "%s        & $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $ \\\\" % (sample,h_wlfWmn.Integral(),1./sqrt(h_wlfWmn.GetEffectiveEntries())*h_wlfWmn.Integral(),h_wlfWen.Integral(),1./sqrt(h_wlfWen.GetEffectiveEntries())*h_wlfWen.Integral(),h_ttWmn.Integral(),h_ttWmn.Integral()*1./sqrt(h_ttWmn.GetEffectiveEntries()),h_ttWen.Integral(),1./sqrt(h_ttWen.GetEffectiveEntries())*h_ttWen.Integral(),h_whfWmn.Integral(),1./sqrt(h_whfWmn.GetEffectiveEntries())*h_whfWmn.Integral(),h_whfWen.Integral(),1./sqrt(h_whfWen.GetEffectiveEntries())*h_whfWen.Integral())
    h_Wmn = ifile_wmn.Get("BDT_WmnHighPt_%s" % sample)
    h_Wen = ifile_wen.Get("BDT_WenHighPt_%s" % sample)
    print sample
    #print "%.2f +/= %.2f" % (h_Wmn.Integral(),h_Wmn.Integral()*1./sqrt(h_Wmn.GetEffectiveEntries()))
    print "%.2f +/= %.2f" % (h_Wen.Integral(),h_Wen.Integral()*1./sqrt(h_Wen.GetEffectiveEntries()))
    nBins = h_Wmn.GetNbinsX()
    totErr = 0.
    totSum = 0.
    for i in range(nBins-3,nBins+1):
        #print i
        totSum += h_Wmn.GetBinContent(i)
        totSum += h_Wen.GetBinContent(i)
        totErr += pow(h_Wmn.GetBinError(i),2)
    totErr = sqrt(totErr)
    #print "%.2f +/= %.2f" % (totSum, totErr)
#print "Total Background        & $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $ \\\\" % (mcTot_0,sqrt(mcTotEnt_0),mcTot_1,sqrt(mcTotEnt_1),mcTot_2,sqrt(mcTotEnt_2),mcTot_3,sqrt(mcTotEnt_3),mcTot_4,sqrt(mcTotEnt_4),mcTot_5,sqrt(mcTotEnt_5))
