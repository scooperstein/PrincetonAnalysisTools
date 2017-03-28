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

#ROOT.gSystem.Load("PU_C.so")

#filename1=sys.argv[1]
#filename2=sys.argv[2]


#tfile1=ROOT.TFile(filename1)
#tfile2=ROOT.TFile(filename2)

treename="tree"

#tree1=tfile1.Get(treename)
#tree2=tfile2.Get(treename)

tree1=ROOT.TChain(treename)
tree2=ROOT.TChain(treename)

#tree1.Add(filename1)
#tree2.Add(filename2)

tree1.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_SR_March13/output_mc.root")
tree2.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_SR_March13_NLOWJets/output_wjetsnlo_4.root")
tree2.Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_SR_March13_NLOWJets/output_wjetsnlo_partonbinned_3.root")

#tree1.Add("/eos/uscms/store/user/capalmer/VHBBHeppyNtuples/V12/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V12_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150723_084710/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V13/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V13_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_084337/0000/*.root")

#tree1.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V13/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V13_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_084337/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V14Test/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_R14_TEST_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151020_145014/0000/*.root")

#tree1.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V13/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V13_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_084337/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V13/WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V13_WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_084407/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V13/WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V13_WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151002_084353/0000/*.root")

#tree1.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V14/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V14_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_220043/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V20test/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_A20_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160201_142009/0000/*.root")

#tree1.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V20/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V20_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172432/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V21/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V21_WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_150209/0000/*.root")

#tree1.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V20/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V20_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_171237/0000/*.root")
#tree2.Add("/eos/uscms/store/group/lpchbb/HeppyNtuples/V21/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V21_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_145703/0000/*.root")

#can=ROOT.TCanvas("can","",800,500)

#can.cd()

#presel = "(Vtype==3 || Vtype==2)*sign(genWeight)"
#presel = "((V_pt>100&&HCSV_pt>100)&&(Vtype==2||Vtype==3))*sign(genWeight)"
#presel = "(evt==50972414)*weight"
#presel = "Jet_btagCSV[hJetInd1]>0.935&&Pass_nominal"0presel = "(H_pt>100&&V_pt>100&&(Vtype==2||Vtype==3)&&H_mass>90&&H_mass<150&&selLeptons_relIso_0<0.06&&cutFlow>=9&&Pass_nominal&&((sampleIndex>=2200&&sampleIndex<=2202)||(sampleIndex>=4100&&sampleIndex<=4902)||(sampleIndex>=7100&&sampleIndex<=7402)||(sampleIndex>=8000&&sampleIndex<=8202&&lheV_pt<150)))*weight*bTagWeightMoriondCMVA*VPtCorrFactorSplit3*(1./CS_SF)"
#presel = "(sampleIndex%100==2&&H_pt>100&&V_pt>150&&(Vtype==2||Vtype==3)&&H_mass>90&&H_mass<150&&selLeptons_relIso_0<0.06&&cutFlow>=9&&Pass_nominal&&((sampleIndex>=2200&&sampleIndex<=2202)||(sampleIndex>=4100&&sampleIndex<=4902)||(sampleIndex>=7100&&sampleIndex<=7402)||(sampleIndex>=8000&&sampleIndex<=8202&&lheV_pt<100)))*weight*bTagWeightMoriondCMVA*VPtCorrFactorSplit3*(1./CS_SF)"
#presel = "(sampleIndex%100==2&&H_pt>100&&V_pt>100&&(Vtype==2||Vtype==3)&&H_mass>90&&H_mass<150&&selLeptons_relIso_0<0.06&&cutFlow>=9&&Pass_nominal&&((sampleIndex>=2200&&sampleIndex<=2202)||(sampleIndex>=4100&&sampleIndex<=4902)||(sampleIndex>=7100&&sampleIndex<=7402)))*weight*VPtCorrFactorSplit3*(1./CS_SF)"
presel = "(Pass_nominal&&V_pt>100&&selLeptons_relIso_0<0.06&&H_mass>90&&H_mass<150&&cutFlow>=10&&((isWmunu&&Vtype==2)||(isWenu&&Vtype==3))&&sampleIndex%100==2&&((sampleIndex>=2200&&sampleIndex<=2202)||(sampleIndex>=4100&&sampleIndex<=4902)||(sampleIndex>=7100&&sampleIndex<=7402)||(sampleIndex>=8000&&sampleIndex<=8202&&lheV_pt<100)))*weight*(1./CS_SF)*bTagWeightMoriondCMVA*VPtCorrFactorSplit3"
#presel = "(H_pt>100&&V_pt>100&&(Vtype==2||Vtype==3)&&H_mass>90&&H_mass<150&&selLeptons_relIso_0<0.06&&cutFlow>=9&&Pass_nominal)*weight"

plotlist={}

#plotlist["JetPt"]              = ["Jet_pt",presel,100,15,400]
#plotlist["Jet1Pt"]              = ["Jet_pt[hJetInd1]",presel,100,15,400]
#plotlist["Jet2Pt"]              = ["Jet_pt[hJetInd2]",presel,100,15,400]
#plotlist["Jet1Eta"]             = ["Jet_eta[hJetInd1]",presel,50,-2.5,2.5]
#plotlist["Jet2Eta"]             = ["Jet_eta[hJetInd2]",presel,50,-2.5,2.5]
#plotlist["Jet1Phi"]             = ["Jet_phi[hJetInd1]",presel,50,-3.2,3.2]
#plotlist["Jet2Phi"]             = ["Jet_phi[hJetInd2]",presel,50,-3.2,3.2]
#plotlist["Jet1Mass"]            = ["Jet_mass[hJetInd1]",presel,100,0,200]
#plotlist["Jet2Mass"]            = ["Jet_mass[hJetInd2]",presel,100,0,200]
#plotlist["JetCSV"]             = ["Jet_btagCSV",presel,100,0,1.]
#plotlist["Jet1CSV"]             = ["Jet_btagCSV[hJetInd1]",presel,20,0.0,1.]
["Jet2CSV"]             = ["Jet_btagCSV[hJetInd2]",presel,20,0.0,1.]
#plotlist["Jet1PuId"]            = ["Jet_puId[hJetInd1]",presel,10,0,10]
#plotlist["Jet2PuId"]            = ["Jet_puId[hJetInd2]",presel,10,0,10]
#plotlist["Jet1Id"]              = ["Jet_id[hJetInd1]",presel,8,0,8]
#plotlist["Jet2Id"]              = ["Jet_id[hJetInd2]",presel,8,0,8]
#plotlist["LepPt"]              = ["selLeptons_pt[0]",presel,100,22,400]
#plotlist["LepEta"]             = ["selLeptons_eta[0]",presel,50,-2.5,2.5]
#plotlist["LepPhi"]             = ["selLeptons_phi[0]",presel,50,-3.2,3.2]
#plotlist["Lep1PdgId"]           = ["abs(selLeptons_pdgId[0])",presel,25,0,25]
#plotlist["Lep1Id"]              = ["selLeptons_tightId[0]",presel,3,0,3]
["MetPt"]              = ["met_pt",presel,40,0,400]
#plotlist["MetPhi"]             = ["met_phi",presel,50,-3.2,3.2]
#plotlist["PuppiMetPt"]         = ["metPuppi_pt",presel,100,0,400]
["HPt"]                = ["H_pt",presel,40,0,400]
#plotlist["HCSVEta"]               = ["H_eta",presel,50,-2.5,2.5]
#plotlist["HCSVPhi"]               = ["HCSV_phi",presel,50,-3.2,3.2]
["HMass"]              = ["H_mass",presel,24,90,150]
["WPt"]                = ["V_pt",presel,30,100,400]
["WMt"]                = ["V_mt",presel,20,0,200]
#plotlist["nJetTot"]            = ["nJet",presel,20,0,20]
#plotlist["NJetPt25Eta2p9"]     = ["Sum$(Jet_pt>25&&Jet_puId>0&&abs(Jet_eta)<2.9)",presel,10,0,10]
["NAddJet"]     = ["nAddJet_f",presel,10,0,10]
#plotlist["NLepPt15Eta2p5"]     = ["Sum$(selLeptons_pt>15&&abs(selLeptons_eta)<2.5)",presel,10,0,10]
#plotlist["NLepPt30Eta2p5"]     = ["Sum$(selLeptons_pt>30&&abs(selLeptons_eta)<2.5)",presel,10,0,10]
#plotlist["HVDPhi"]             = ["HVdPhi",presel,50,-3.2,3.2]
["HVDPhi"]             = ["abs(HVdPhi)",presel,20,2.5,3.2]
#plotlist["MetRawPt"]           = ["met_rawpt",presel,100,0,400]
["TopMass"] = ["Top1_mass_fromLepton_regPT_w4MET",presel,20,100,500]
["lepMetDPhi"] = ["lepMetDPhi",presel,20,0,2.0]
["nSA5Jet"] = ["softActivityVH_njets5",presel,10,0,10]
plotlist["BDT"] = ["CMS_vhbb_BDT_Wln_13TeV",presel,15,-1,1]
#plotlist["HJ1_HJ2_dR"] = ["abs(HJ1_HJ2_dR)",presel,40,0,5]


plots=plotlist.keys()
plots.sort()

for plot in plots:
    print plot
    temp1=ROOT.TH1F("temp1",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    tree1.Draw(plotlist[plot][0]+">>temp1",plotlist[plot][1])
    temp2=ROOT.TH1F("temp2",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    #tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1].replace("V_pt>100","V_pt>100&&CMS_vhbb_BDT_Wln_13TeV>0.608"))
    #tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1].replace("V_pt>100","V_pt>100&&CMS_vhbb_BDT_Wln_13TeV<0.608&&CMS_vhbb_BDT_Wln_13TeV>0.376"))
    tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1])

    temp1.SetLineColor(ROOT.kRed)
    temp2.SetLineColor(ROOT.kBlue)
   
    from math import sqrt
    for i in range(1,temp1.GetNbinsX()+1):
        B = temp1.GetBinContent(i)
        dB = temp1.GetBinError(i)
        if B > 0:
            print "Bin ",i,": ",B,dB,(dB/sqrt(B))
    for i in range(1,temp2.GetNbinsX()+1):
        B = temp2.GetBinContent(i)
        dB = temp2.GetBinError(i)
        if B > 0:
            print "Bin ",i,": ",B,dB,(dB/sqrt(B))

    #if (temp1.Integral() > 0):
    #    temp1.Scale(1./temp1.Integral())
    #if (temp2.Integral() > 0):
    #    temp2.Scale(1./temp2.Integral())

    maxval=max(temp1.GetMaximum(),temp2.GetMaximum())*1.2

    canv = ROOT.TCanvas("canv","canv")
    canv.SetBottomMargin(0.3)
    xmin = float(plotlist[plot][3])
    xmax = float(plotlist[plot][4])
    frame = ROOT.TH1F("frame","",1,xmin,xmax)
    frame.SetStats(0)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleOffset(0.91);
    frame.GetYaxis().SetTitle("Events")
    #frame.GetXaxis().SetTitle("atanhBDT")
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetRangeUser(0.,maxval)
    #frame.GetYaxis().SetRangeUser(10E-03,stack.GetStack().Last().GetMaximum()*100)
    #frame.GetYaxis().SetRangeUser(10E-03,200)
    frame.Draw()


    print "max1, max2,",temp1.GetMaximum(),temp2.GetMaximum()
    print "integral1, integral2,",temp1.Integral(),temp2.Integral()
    temp1.SetMaximum(maxval)
    temp1.Draw("epsame")
    temp2.Draw("epsame")

    leg=ROOT.TLegend(0.7,0.7,1.,1.)
    leg.AddEntry(temp1,"W+jets LO Weighted")
    leg.AddEntry(temp2,"W+jets NLO")
    #leg.AddEntry(temp2,"WH(bb) Powheg BDT > 0.608")
    leg.Draw()

    pad2 = ROOT.TPad("pad2", "pad2", 0., 0., 1., 1.)
    pad2.SetTopMargin(0.73)
    #       pad2.SetBottomMargin(0.)
    #       pad2.SetRightMargin(0.)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    pad2.Draw()
    pad2.cd()
    frame2 = ROOT.TH1F("frame2","",1,xmin,xmax)
    frame2.SetMinimum(-0.5)
    frame2.SetMaximum(0.5)
    frame2.GetYaxis().SetLabelSize(0.02)
    frame2.GetXaxis().SetLabelSize(0.04)
    frame2.GetYaxis().SetTitleSize(0.04)
    #frame2.GetXaxis().SetTitle("BDT Score")
    #frame2.GetXaxis().SetTitle("Sub-leading Jet CMVA")
    frame2.SetStats(0)
    frame2.GetYaxis().SetTitle("LO/NLO - 1")
    frame2.Draw()

    data2=temp2.Clone("data2")
    mc_total=temp1.Clone("mc_total")
    data2.Add(mc_total,-1)
    data2.Divide(mc_total)

    mc_total_uncUp=temp1.Clone("mc_total")
    mc_total_uncDown=temp1.Clone("mc_total")
    for i in range(0,mc_total_uncUp.GetNbinsX()):
        e=0.
        if mc_total.GetBinContent(i+1) != 0:
                e = mc_total.GetBinError(i+1)/mc_total.GetBinContent(i+1)
        mc_total_uncUp.SetBinContent(i+1,e)
        mc_total_uncDown.SetBinContent(i+1,-e)
    mc_total_uncUp.SetLineColor(ROOT.kBlack)
    mc_total_uncUp.SetLineWidth(1)
    mc_total_uncUp.SetFillColor(ROOT.kBlack)
    mc_total_uncUp.SetFillStyle(3004)
    mc_total_uncDown.SetLineColor(ROOT.kBlack)
    mc_total_uncDown.SetLineWidth(1)
    mc_total_uncDown.SetFillColor(ROOT.kBlack)
    mc_total_uncDown.SetFillStyle(3004)

    mc_total_uncUp.Draw("HIST SAME")
    mc_total_uncDown.Draw("HIST SAME") 
    line = ROOT.TLine(xmin,0,xmax,0)
    line.SetLineStyle(3)
    line.Draw("same")

    data2.Draw("PEsame")

    chi2 = temp1.Chi2Test(temp2,"WWCHI2/NDF")
    chi2text = ROOT.TPaveText(0.25,0.8,0.35,0.85,"brNDC");
    chi2text.SetTextFont(62);
    chi2text.SetTextSize(0.04)
    chi2text.SetBorderSize(0)
    chi2text.SetLineColor(0)
    chi2text.SetLineStyle(0)
    chi2text.SetLineWidth(0)
    chi2text.SetFillColor(0)
    chi2text.SetFillStyle(0)
    chi2text.AddText("#chi^{2}_{ }#lower[0.1]{/^{}#it{dof} = %.2f}"%(chi2))
    chi2text.Draw()

    canv.Update()
    canv.SaveAs(plot+".png")
    #raw_input()
