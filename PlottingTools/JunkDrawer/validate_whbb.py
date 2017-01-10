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

filename1=sys.argv[1]
filename2=sys.argv[2]


#tfile1=ROOT.TFile(filename1)
#tfile2=ROOT.TFile(filename2)

treename="tree"

#tree1=tfile1.Get(treename)
#tree2=tfile2.Get(treename)

tree1=ROOT.TChain(treename)
tree2=ROOT.TChain(treename)

tree1.Add(filename1)
tree2.Add(filename2)

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

can=ROOT.TCanvas("can","",800,500)

can.cd()

#presel = "(Vtype==3 || Vtype==2)*sign(genWeight)"
#presel = "((V_pt>100&&HCSV_pt>100)&&(Vtype==2||Vtype==3))*sign(genWeight)*puWeight"
#presel = "Jet_btagCSV[hJetInd1]>0.935&&Pass_nominal"
presel = "(sampleIndex==-12501&&H_pt>100&&V_pt>100&&(Vtype==2||Vtype==3)&&run<=276811&&H_mass>90&&H_mass<150&&selLeptons_relIso_0<0.06&&cutFlow>=9&&Pass_nominal)*weight"

plotlist={}

plotlist["JetPt"]              = ["Jet_pt",presel,100,15,400]
plotlist["Jet1Pt"]              = ["Jet_pt[hJetInd1]",presel,100,15,400]
plotlist["Jet2Pt"]              = ["Jet_pt[hJetInd2]",presel,100,15,400]
plotlist["Jet1Eta"]             = ["Jet_eta[hJetInd1]",presel,50,-2.5,2.5]
plotlist["Jet2Eta"]             = ["Jet_eta[hJetInd2]",presel,50,-2.5,2.5]
plotlist["Jet1Phi"]             = ["Jet_phi[hJetInd1]",presel,50,-3.2,3.2]
plotlist["Jet2Phi"]             = ["Jet_phi[hJetInd2]",presel,50,-3.2,3.2]
plotlist["Jet1Mass"]            = ["Jet_mass[hJetInd1]",presel,100,0,200]
plotlist["Jet2Mass"]            = ["Jet_mass[hJetInd2]",presel,100,0,200]
plotlist["JetCSV"]             = ["Jet_btagCSV",presel,100,0,1.]
plotlist["Jet1CSV"]             = ["Jet_btagCSV[hJetInd1]",presel,20,0.935,1.]
plotlist["Jet2CSV"]             = ["Jet_btagCSV[hJetInd2]",presel,20,0.46,1.]
#plotlist["Jet1PuId"]            = ["Jet_puId[hJetInd1]",presel,2,0,2]
#plotlist["Jet2PuId"]            = ["Jet_puId[hJetInd2]",presel,2,0,2]
plotlist["Jet1Id"]              = ["Jet_id[hJetInd1]",presel,8,0,8]
plotlist["Jet2Id"]              = ["Jet_id[hJetInd2]",presel,8,0,8]
plotlist["Lep1Pt"]              = ["selLeptons_pt[0]",presel,100,22,400]
plotlist["Lep1Eta"]             = ["selLeptons_eta[0]",presel,50,-2.5,2.5]
plotlist["Lep1Phi"]             = ["selLeptons_phi[0]",presel,50,-3.2,3.2]
#plotlist["Lep1PdgId"]           = ["abs(selLeptons_pdgId[0])",presel,25,0,25]
#plotlist["Lep1Id"]              = ["selLeptons_tightId[0]",presel,3,0,3]
plotlist["MetPt"]              = ["met_pt",presel,100,0,400]
plotlist["MetPhi"]             = ["met_phi",presel,50,-3.2,3.2]
#plotlist["PuppiMetPt"]         = ["metPuppi_pt",presel,100,0,400]
plotlist["HCSVPt"]                = ["H_pt",presel,100,0,400]
#plotlist["HCSVEta"]               = ["H_eta",presel,50,-2.5,2.5]
#plotlist["HCSVPhi"]               = ["H_phi",presel,50,-3.2,3.2]
plotlist["HCSVMass"]              = ["H_mass",presel,100,0,500]
plotlist["WPt"]                = ["V_pt",presel,100,0,400]
plotlist["WMt"]                = ["V_mt",presel,100,0,400]
plotlist["nJetTot"]            = ["nJet",presel,20,0,20]
#plotlist["NJetPt25Eta2p9"]     = ["Sum$(Jet_pt>25&&Jet_puId>0&&abs(Jet_eta)<2.9)",presel,10,0,10]
plotlist["NAddJet"]     = ["nAddJet_f",presel,10,0,10]
#plotlist["NLepPt15Eta2p5"]     = ["Sum$(selLeptons_pt>15&&abs(selLeptons_eta)<2.5)",presel,10,0,10]
#plotlist["NLepPt30Eta2p5"]     = ["Sum$(selLeptons_pt>30&&abs(selLeptons_eta)<2.5)",presel,10,0,10]
#plotlist["HVDPhi"]             = ["HVdPhi",presel,50,-3.2,3.2]
plotlist["HVDPhi"]             = ["abs(HVdPhi)",presel,20,2.5,3.2]
#plotlist["MetRawPt"]           = ["met_rawpt",presel,100,0,400]
#plotlist["Jet1PuId"]           = ["Jet_puId[0]",presel,2,0,2]
#plotlist["Jet2PuId"]           = ["Jet_puId[1]",presel,2,0,2]
plotlist["TopMass"] = ["Top1_mass_fromLepton_regPT_w4MET",presel,20,100,500]
plotlist["lepMetDPhi"] = ["lepMetDPhi",presel,20,0,3.15]
plotlist["nSA5Jet"] = ["softActivityVH_njets5",presel,10,0,10]
plotlist["HJ1_HJ2_dR"] = ["abs(HJ1_HJ2_dR)",presel,40,0,5]


plots=plotlist.keys()
plots.sort()

for plot in plots:
    print plot
    temp1=ROOT.TH1F("temp1",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    tree1.Draw(plotlist[plot][0]+">>temp1",plotlist[plot][1])
    temp2=ROOT.TH1F("temp2",plot,plotlist[plot][2],plotlist[plot][3],plotlist[plot][4])
    #tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1].replace("V_pt>100","V_pt>100&&CMS_vhbb_BDT_Wln_13TeV>0.608"))
    tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1].replace("V_pt>100","V_pt>100&&CMS_vhbb_BDT_Wln_13TeV<0.608&&CMS_vhbb_BDT_Wln_13TeV>0.376"))
    #tree2.Draw(plotlist[plot][0]+">>temp2",plotlist[plot][1])

    temp1.SetLineColor(ROOT.kRed)
    temp2.SetLineColor(ROOT.kBlue)

    #temp1.Scale(1./temp1.Integral())
    #temp2.Scale(1./temp2.Integral())

    maxval=max(temp1.GetMaximum(),temp2.GetMaximum())*1.2

    print "max1, max2,",temp1.GetMaximum(),temp2.GetMaximum()
    print "integral1, integral2,",temp1.Integral(),temp2.Integral()
    temp1.SetMaximum(maxval)
    temp1.Draw()
    temp2.Draw("sames")

    leg=ROOT.TLegend(0.7,0.7,1.,1.)
    leg.AddEntry(temp1,"WH(bb) Powheg")
    leg.AddEntry(temp2,"WH(bb) Powheg 0.376 < BDT < 0.608")
    #leg.AddEntry(temp2,"WH(bb) Powheg BDT > 0.608")
    leg.Draw()

    can.Update()
    can.SaveAs(plot+".png")
    #raw_input()
