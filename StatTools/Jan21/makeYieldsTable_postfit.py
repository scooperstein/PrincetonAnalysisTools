import ROOT
from math import sqrt
from math import pow

ifile = ROOT.TFile("mlfit.root")

mcTot_0 = 0.
mcTotErr_0 = 0
mcTot_1 = 0.
mcTotErr_1 = 0
mcTot_2 = 0.
mcTotErr_2 = 0
mcTot_3 = 0.
mcTotErr_3 = 0
mcTot_4 = 0.
mcTotErr_4 = 0
mcTot_5 = 0.
mcTotErr_5 = 0
#for sample in ["data_obs"]:
#for sample in ["Zj0b","Zj1b","Zj2b","Wj0b","Wj1b","Wj2b","TT","s_Top","VVHF","VVLF"]:
for sample in ["Zj0b","Zj1b","Zj2b","Wj0b","Wj1b","Wj2b","TT","s_Top","VVLF","VVHF","ZH","WH"]:
#for sample in ["Zj0b","Zj1b","Zj2b","Wj0b","Wj1b","Wj2b","TT","s_Top","VVLF","VVHF"]:
    h_ttWmn = ifile.Get("shapes_fit_s/ttWmn/%s" % sample)
    mcTot_2 += h_ttWmn.Integral()
    mcErr_2 = 0. 
    mcErr_3 = 0. 
    mcErr_4 = 0. 
    mcErr_5 = 0. 
    mcErr_0 = 0. 
    mcErr_1 = 0. 
    for i in range(1,h_ttWmn.GetNbinsX()+1):
        mcErr_2 += pow(h_ttWmn.GetBinError(i),2)
    mcTotErr_2 += mcErr_2
    mcErr_2 = sqrt(mcErr_2)
    h_ttWen = ifile.Get("shapes_fit_s/ttWen/%s" % sample)
    mcTot_3 += h_ttWen.Integral()
    for i in range(1,h_ttWen.GetNbinsX()+1):
        mcErr_3 += pow(h_ttWen.GetBinError(i),2)
    mcTotErr_3 += mcErr_3
    mcErr_3 = sqrt(mcErr_3)
    h_whfWmn = ifile.Get("shapes_fit_s/whfWmn/%s" % sample)
    mcTot_4 += h_whfWmn.Integral()
    for i in range(1,h_whfWmn.GetNbinsX()+1):
        mcErr_4 += pow(h_whfWmn.GetBinError(i),2)
    mcTotErr_4 += mcErr_4
    mcErr_4 = sqrt(mcErr_4)
    h_whfWen = ifile.Get("shapes_fit_s/whfWen/%s" % sample)
    mcTot_5 += h_whfWen.Integral()
    for i in range(1,h_whfWen.GetNbinsX()+1):
        mcErr_5 += pow(h_whfWen.GetBinError(i),2)
    mcTotErr_5 += mcErr_5
    mcErr_5 = sqrt(mcErr_5)
    h_wlfWmn = ifile.Get("shapes_fit_s/wlfWmn/%s" % sample)
    mcTot_0 += h_wlfWmn.Integral()
    for i in range(1,h_wlfWmn.GetNbinsX()+1):
        mcErr_0 += pow(h_wlfWmn.GetBinError(i),2)
    mcTotErr_0 += mcErr_0
    mcErr_0 = sqrt(mcErr_0)
    h_wlfWen = ifile.Get("shapes_fit_s/wlfWen/%s" % sample) 
    mcTot_1 += h_wlfWen.Integral()
    for i in range(1,h_wlfWen.GetNbinsX()+1):
        mcErr_1 += pow(h_wlfWen.GetBinError(i),2)
    mcTotErr_1 += mcErr_1
    mcErr_1 = sqrt(mcErr_1)
    print "%s        & $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $ \\\\" % (sample,h_wlfWmn.Integral(),mcErr_0,h_wlfWen.Integral(),mcErr_1,h_ttWmn.Integral(),mcErr_2,h_ttWen.Integral(),mcErr_3,h_whfWmn.Integral(),mcErr_4,h_whfWen.Integral(),mcErr_5)
    #print "%s        & $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $ \\\\" % (sample,h_wlfWmn.Integral(),1./sqrt(h_wlfWmn.GetEffectiveEntries())*h_wlfWmn.Integral(),h_wlfWen.Integral(),1./sqrt(h_wlfWen.GetEffectiveEntries())*h_wlfWen.Integral(),h_ttWmn.Integral(),h_ttWmn.Integral()*1./sqrt(h_ttWmn.GetEffectiveEntries()),h_ttWen.Integral(),1./sqrt(h_ttWen.GetEffectiveEntries())*h_ttWen.Integral(),h_whfWmn.Integral(),1./sqrt(h_whfWmn.GetEffectiveEntries())*h_whfWmn.Integral(),h_whfWen.Integral(),1./sqrt(h_whfWen.GetEffectiveEntries())*h_whfWen.Integral())
    #h_Wmn = ifile.Get("shapes_fit_s/WmnHighPt/%s" % sample)
    #h_Wen = ifile.Get("shapes_fit_s/WenHighPt/%s" % sample)
    #print sample
    #print "%.2f +/= %.2f" % (h_Wmn.Integral(),h_Wmn.Integral()*1./sqrt(h_Wmn.GetEffectiveEntries()))
    #print "%.2f +/= %.2f" % (h_Wen.Integral(),h_Wen.Integral()*1./sqrt(h_Wen.GetEffectiveEntries()))
    #nBins = h_Wmn.GetNbinsX()
    #totErr = 0.
    #totSum = 0.
    #for i in range(nBins-3,nBins+1):
    #    #print i
    #    totSum += h_Wmn.GetBinContent(i)
    #    totSum += h_Wen.GetBinContent(i)
    #    totErr += pow(h_Wmn.GetBinError(i),2)
    #totErr = sqrt(totErr)
    #print "%.2f +/= %.2f" % (totSum, totErr)
print "Total Background        & $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $   &  $ %.2f \pm %.2f $ \\\\" % (mcTot_0,sqrt(mcTotErr_0),mcTot_1,sqrt(mcTotErr_1),mcTot_2,sqrt(mcTotErr_2),mcTot_3,sqrt(mcTotErr_3),mcTot_4,sqrt(mcTotErr_4),mcTot_5,sqrt(mcTotErr_5))
