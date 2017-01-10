## Compare systematic shape variation for given shape systematic vs. nominal for BDT
##
## Author: Stephane Cooperstein
##

import argparse
import ROOT
import TdrStyles

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser("Make systematic shape comparison plot")
parser.add_argument('-i', '--inputfile', type=str, default="", help="The input root file w/ shapes")
parser.add_argument('-c', '--channel', type=str, default="", help="The input root file w/ shapes")
#parser.add_argument('-s', '--systematic', type=str, default="", help="the name of the shape systematic")
#parser.add_argument('-v','--variable', type=str, default="CMS_vhbb_BDT_Wln_13TeV", help="the variable to plot")
args = parser.parse_args()

samples = ["WH","ZH","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","Zj0b","Zj1b","Zj2b"]
#samples = ["WH","ZH","TT","Wj0b","Wj1b","Wj2b","s_Top","VVHF","VVLF","Zj0b","Zj1b","Zj2b"]
#samples = ["TT"]
#systematics = ["CMS_vhbb_statWj1b_Comb_bin16_13TeV"]
#systematics = ["CMS_vhbb_scale_j_13TeV","CMS_vhbb_res_j_13TeV"]
#systematics =  ["CMS_vhbb_puWeight"]
#systematics = ["CMS_vhbb_LHE_weights_scale_muR_Wj2b","CMS_vhbb_LHE_weights_scale_muR_Wj1b","CMS_vhbb_LHE_weights_scale_muR_Wj0b","CMS_vhbb_LHE_weights_scale_muR_TT","CMS_vhbb_LHE_weights_scale_muR_VVHF","CMS_vhbb_LHE_weights_scale_muR_VVLF","CMS_vhbb_LHE_weights_scale_muR_ZH","CMS_vhbb_LHE_weights_scale_muR_WH","CMS_vhbb_LHE_weights_scale_muR_Zj2b","CMS_vhbb_LHE_weights_scale_muR_Zj1b","CMS_vhbb_LHE_weights_scale_muR_Zj0b","CMS_vhbb_LHE_weights_scale_muF_Wj2b","CMS_vhbb_LHE_weights_scale_muF_Wj1b","CMS_vhbb_LHE_weights_scale_muF_Wj0b","CMS_vhbb_LHE_weights_scale_muF_TT","CMS_vhbb_LHE_weights_scale_muF_VVHF","CMS_vhbb_LHE_weights_scale_muF_VVLF","CMS_vhbb_LHE_weights_scale_muF_ZH","CMS_vhbb_LHE_weights_scale_muF_WH","CMS_vhbb_LHE_weights_scale_muF_Zj2b","CMS_vhbb_LHE_weights_scale_muF_Zj1b","CMS_vhbb_LHE_weights_scale_muF_Zj0b"]
#systematics = ["CMS_vhbb_LHE_weights_pdf_30"]
systematics = ["CMS_vhbb_LHE_weights_scale_muR_Wj2b"]
#systematics = ['CMS_vhbb_LHE_weights_pdf_0', 'CMS_vhbb_LHE_weights_pdf_1', 'CMS_vhbb_LHE_weights_pdf_2', 'CMS_vhbb_LHE_weights_pdf_3', 'CMS_vhbb_LHE_weights_pdf_4', 'CMS_vhbb_LHE_weights_pdf_5', 'CMS_vhbb_LHE_weights_pdf_6', 'CMS_vhbb_LHE_weights_pdf_7', 'CMS_vhbb_LHE_weights_pdf_8', 'CMS_vhbb_LHE_weights_pdf_9', 'CMS_vhbb_LHE_weights_pdf_10', 'CMS_vhbb_LHE_weights_pdf_11', 'CMS_vhbb_LHE_weights_pdf_12', 'CMS_vhbb_LHE_weights_pdf_13', 'CMS_vhbb_LHE_weights_pdf_14', 'CMS_vhbb_LHE_weights_pdf_15', 'CMS_vhbb_LHE_weights_pdf_16', 'CMS_vhbb_LHE_weights_pdf_17', 'CMS_vhbb_LHE_weights_pdf_18', 'CMS_vhbb_LHE_weights_pdf_19', 'CMS_vhbb_LHE_weights_pdf_20', 'CMS_vhbb_LHE_weights_pdf_21', 'CMS_vhbb_LHE_weights_pdf_22', 'CMS_vhbb_LHE_weights_pdf_23', 'CMS_vhbb_LHE_weights_pdf_24', 'CMS_vhbb_LHE_weights_pdf_25', 'CMS_vhbb_LHE_weights_pdf_26', 'CMS_vhbb_LHE_weights_pdf_27', 'CMS_vhbb_LHE_weights_pdf_28', 'CMS_vhbb_LHE_weights_pdf_29', 'CMS_vhbb_LHE_weights_pdf_30', 'CMS_vhbb_LHE_weights_pdf_31', 'CMS_vhbb_LHE_weights_pdf_32', 'CMS_vhbb_LHE_weights_pdf_33', 'CMS_vhbb_LHE_weights_pdf_34', 'CMS_vhbb_LHE_weights_pdf_35', 'CMS_vhbb_LHE_weights_pdf_36', 'CMS_vhbb_LHE_weights_pdf_37', 'CMS_vhbb_LHE_weights_pdf_38', 'CMS_vhbb_LHE_weights_pdf_39', 'CMS_vhbb_LHE_weights_pdf_40', 'CMS_vhbb_LHE_weights_pdf_41', 'CMS_vhbb_LHE_weights_pdf_42', 'CMS_vhbb_LHE_weights_pdf_43', 'CMS_vhbb_LHE_weights_pdf_44', 'CMS_vhbb_LHE_weights_pdf_45', 'CMS_vhbb_LHE_weights_pdf_46', 'CMS_vhbb_LHE_weights_pdf_47', 'CMS_vhbb_LHE_weights_pdf_48', 'CMS_vhbb_LHE_weights_pdf_49']
#systematics2 = ['CMS_vhbb_LHE_weights_pdf_50', 'CMS_vhbb_LHE_weights_pdf_51', 'CMS_vhbb_LHE_weights_pdf_52', 'CMS_vhbb_LHE_weights_pdf_53', 'CMS_vhbb_LHE_weights_pdf_54', 'CMS_vhbb_LHE_weights_pdf_55', 'CMS_vhbb_LHE_weights_pdf_56', 'CMS_vhbb_LHE_weights_pdf_57', 'CMS_vhbb_LHE_weights_pdf_58', 'CMS_vhbb_LHE_weights_pdf_59', 'CMS_vhbb_LHE_weights_pdf_60', 'CMS_vhbb_LHE_weights_pdf_61', 'CMS_vhbb_LHE_weights_pdf_62', 'CMS_vhbb_LHE_weights_pdf_63', 'CMS_vhbb_LHE_weights_pdf_64', 'CMS_vhbb_LHE_weights_pdf_65', 'CMS_vhbb_LHE_weights_pdf_66', 'CMS_vhbb_LHE_weights_pdf_67', 'CMS_vhbb_LHE_weights_pdf_68', 'CMS_vhbb_LHE_weights_pdf_69', 'CMS_vhbb_LHE_weights_pdf_70', 'CMS_vhbb_LHE_weights_pdf_71', 'CMS_vhbb_LHE_weights_pdf_72', 'CMS_vhbb_LHE_weights_pdf_73', 'CMS_vhbb_LHE_weights_pdf_74', 'CMS_vhbb_LHE_weights_pdf_75', 'CMS_vhbb_LHE_weights_pdf_76', 'CMS_vhbb_LHE_weights_pdf_77', 'CMS_vhbb_LHE_weights_pdf_78', 'CMS_vhbb_LHE_weights_pdf_79', 'CMS_vhbb_LHE_weights_pdf_80', 'CMS_vhbb_LHE_weights_pdf_81', 'CMS_vhbb_LHE_weights_pdf_82', 'CMS_vhbb_LHE_weights_pdf_83', 'CMS_vhbb_LHE_weights_pdf_84', 'CMS_vhbb_LHE_weights_pdf_85', 'CMS_vhbb_LHE_weights_pdf_86', 'CMS_vhbb_LHE_weights_pdf_87', 'CMS_vhbb_LHE_weights_pdf_88', 'CMS_vhbb_LHE_weights_pdf_89', 'CMS_vhbb_LHE_weights_pdf_90', 'CMS_vhbb_LHE_weights_pdf_91', 'CMS_vhbb_LHE_weights_pdf_92', 'CMS_vhbb_LHE_weights_pdf_93', 'CMS_vhbb_LHE_weights_pdf_94', 'CMS_vhbb_LHE_weights_pdf_95', 'CMS_vhbb_LHE_weights_pdf_96', 'CMS_vhbb_LHE_weights_pdf_97', 'CMS_vhbb_LHE_weights_pdf_98', 'CMS_vhbb_LHE_weights_pdf_99', 'CMS_vhbb_LHE_weights_pdf_100']
#systematics.extend(systematics2)
#systematics = ["CMS_vhbb_LHE_weights_pdf_TT","CMS_vhbb_LHE_weights_pdf_Wj0b","CMS_vhbb_LHE_weights_pdf_Wj1b","CMS_vhbb_LHE_weights_pdf_Wj2b","CMS_vhbb_LHE_weights_pdf_WH","CMS_vhbb_LHE_weights_pdf_ZH","CMS_vhbb_LHE_weights_pdf_VVHF","CMS_vhbb_LHE_weights_pdf_VVLF","CMS_vhbb_LHE_weights_pdf_Zj0b","CMS_vhbb_LHE_weights_pdf_Zj1b","CMS_vhbb_LHE_weights_pdf_Zj2b"]
#systematics = ["CMS_vhbb_puWeight","CMS_vhbb_LHE_weights_scale_muR_Wj2b","CMS_vhbb_LHE_weights_scale_muR_Wj1b","CMS_vhbb_LHE_weights_scale_muR_Wj0b","CMS_vhbb_LHE_weights_scale_muR_TT","CMS_vhbb_LHE_weights_scale_muF_Wj2b","CMS_vhbb_LHE_weights_scale_muF_Wj1b","CMS_vhbb_LHE_weights_scale_muF_Wj0b","CMS_vhbb_LHE_weights_scale_muF_TT","CMS_vhbb_LHE_weights_pdf_TT","CMS_vhbb_LHE_weights_pdf_Wj0b","CMS_vhbb_LHE_weights_pdf_Wj1b","CMS_vhbb_LHE_weights_pdf_Wj2b","CMS_vhbb_LHE_weights_pdf_WH","CMS_vhbb_LHE_weights_pdf_ZH","CMS_vhbb_LHE_weights_pdf_VVHF","CMS_vhbb_LHE_weights_pdf_VVLF","CMS_vhbb_LHE_weights_pdf_Zj0b","CMS_vhbb_LHE_weights_pdf_Zj1b","CMS_vhbb_LHE_weights_pdf_Zj2b"]
#systematics = ["CMS_vhbb_LHE_weights_scale_muF_Wj2b","CMS_vhbb_LHE_weights_scale_muF_Wj1b","CMS_vhbb_LHE_weights_scale_muF_Wj0b","CMS_vhbb_LHE_weights_scale_muF_TT"]
#systematics = ["CMS_vhbb_LHE_weights_scale_muR_Wj2b","CMS_vhbb_LHE_weights_scale_muR_Wj1b","CMS_vhbb_LHE_weights_scale_muR_Wj0b","CMS_vhbb_LHE_weights_scale_muR_TT"]
#systematics = ["CMS_vhbb_scale_high_central_j_13TeV","CMS_vhbb_scale_low_central_j_13TeV","CMS_vhbb_scale_high_forward_j_13TeV","CMS_vhbb_scale_low_forward_j_13TeV","CMS_vhbb_res_high_central_j_13TeV","CMS_vhbb_res_low_central_j_13TeV","CMS_vhbb_res_high_forward_j_13TeV","CMS_vhbb_res_high_central_j_13TeV"]
#systematics = ['CMS_vhbb_bTagWeightHF_high_central', 'CMS_vhbb_bTagWeightHF_high_forward', 'CMS_vhbb_bTagWeightHF_low_central', 'CMS_vhbb_bTagWeightHF_low_forward', 'CMS_vhbb_bTagWeightLF_high_central', 'CMS_vhbb_bTagWeightLF_high_forward', 'CMS_vhbb_bTagWeightLF_low_central', 'CMS_vhbb_bTagWeightLF_low_forward', 'CMS_vhbb_bTagWeightLFStats1_high_central', 'CMS_vhbb_bTagWeightLFStats1_high_forward', 'CMS_vhbb_bTagWeightLFStats1_low_central', 'CMS_vhbb_bTagWeightLFStats1_low_forward', 'CMS_vhbb_bTagWeightHFStats1_high_central', 'CMS_vhbb_bTagWeightHFStats1_high_forward', 'CMS_vhbb_bTagWeightHFStats1_low_central', 'CMS_vhbb_bTagWeightHFStats1_low_forward', 'CMS_vhbb_bTagWeightLFStats2_high_central', 'CMS_vhbb_bTagWeightLFStats2_high_forward', 'CMS_vhbb_bTagWeightLFStats2_low_central', 'CMS_vhbb_bTagWeightLFStats2_low_forward', 'CMS_vhbb_bTagWeightHFStats2_high_central', 'CMS_vhbb_bTagWeightHFStats2_high_forward', 'CMS_vhbb_bTagWeightHFStats2_low_central', 'CMS_vhbb_bTagWeightHFStats2_low_forward', 'CMS_vhbb_bTagWeightcErr1_high_central', 'CMS_vhbb_bTagWeightcErr1_high_forward', 'CMS_vhbb_bTagWeightcErr1_low_central', 'CMS_vhbb_bTagWeightcErr1_low_forward', 'CMS_vhbb_bTagWeightcErr2_high_central', 'CMS_vhbb_bTagWeightcErr2_high_forward', 'CMS_vhbb_bTagWeightcErr2_low_central', 'CMS_vhbb_bTagWeightcErr2_low_forward']
#systematics = ["CMS_vhbb_scale_j_13TeV","CMS_vhbb_res_j_13TeV","CMS_vhbb_bTagWeightHF","CMS_vhbb_bTagWeightLF","CMS_vhbb_bTagWeightLFStats1","CMS_vhbb_bTagWeightHFStats1","CMS_vhbb_bTagWeightLFStats2","CMS_vhbb_bTagWeightHFStats2","CMS_vhbb_btagWeightcErr1","CMS_vhbb_btagWeightcErr2","CMS_vhbb_Wj0bModel_Wln_13TeV","CMS_vhbb_Wj1bModel_Wln_13TeV","CMS_vhbb_Wj2bModel_Wln_13TeV","CMS_vhbb_WHModel_Wln_13TeV"]
channel = args.channel
ifile = ROOT.TFile(args.inputfile.replace("WmnCat0",channel), "r")
#c = ROOT.TCanvas("c","c")
for systematic in systematics:
    for sample in samples:
        print systematic,sample
        #if (systematic == "CMS_vhbb_TTModel_Wln_13TeV" and sample != "TT"): continue
        if ((systematic.find("Model") != -1 or systematic.find("LHE") != -1) and systematic.find(sample) == -1): continue
        hnominal = ifile.Get("BDT_%s_%s" % (channel,sample))
        hUp = ifile.Get("BDT_%s_%s_%sUp" % (channel,sample,systematic))
        hDown = ifile.Get("BDT_%s_%s_%sDown" % (channel,sample,systematic))

        hUp.SetLineColor(ROOT.kRed)
        hDown.SetLineColor(ROOT.kBlue)
        hnominal.SetLineColor(ROOT.kBlack)
 
        hUp.SetTitle("Shape Systematic: %s, Sample: %s" % (systematic,sample))

        hUp.SetMaximum(1.2*max(hUp.GetMaximum(),hDown.GetMaximum(),hnominal.GetMaximum()))
  
        TdrStyles.tdrStyle()
        c = ROOT.TCanvas('','', 600, 600)
        c.SetFillStyle(4000)
        c.SetFrameFillStyle(1000)
        c.SetFrameFillColor(0)
        oben = ROOT.TPad('oben','oben',0,0.3 ,1.0,1.0)
        oben.SetBottomMargin(0)
        oben.SetFillStyle(4000)
        oben.SetFrameFillStyle(1000)
        oben.SetFrameFillColor(0)
        unten = ROOT.TPad('unten','unten',0,0.0,1.0,0.3)
        unten.SetTopMargin(0.)
        unten.SetBottomMargin(0.35)
        unten.SetFillStyle(4000)
        unten.SetFrameFillStyle(1000)
        unten.SetFrameFillColor(0)
        oben.Draw()
        unten.Draw()
        oben.cd()

        hUp.Draw("hist")
        hnominal.Draw("ep0 same")
        hDown.Draw("hist same")
        c.Update()
        print sample,hnominal.Integral()
        l = ROOT.TLegend(0.75, 0.8, 0.95, 0.65)
        l.SetLineWidth(2)
        l.SetBorderSize(0)
        l.SetFillColor(0)
        l.SetFillStyle(4000)
        l.SetTextFont(62)
        l.SetTextSize(0.035)
        l.AddEntry(hnominal,"nominal","PL")
        l.AddEntry(hUp, "up","PL")
        l.AddEntry(hDown, "down","PL")
     
        l.Draw("same")

        unten.cd()
        ROOT.gPad.SetTicks(1,1)
 
        ratioU = hUp.Clone("ratioU")
        ratioU.Divide(hnominal)
        ratioD = hDown.Clone("ratioD")
        ratioD.Divide(hnominal)
 
        ratioU.SetStats(0)
        ratioU.GetYaxis().SetRangeUser(0.9,1.10)
        #ratioU.GetYaxis().SetRangeUser(0.95,1.05)
        #ratioU.GetYaxis().SetRangeUser(0.5,1.5)
        ratioU.GetYaxis().SetNdivisions(502,0)
        ratioD.SetStats(0)
        #ratioD.GetYaxis().SetRangeUser(0.95,1.05)
        ratioD.GetYaxis().SetRangeUser(0.9,1.10)
        #ratioD.GetYaxis().SetRangeUser(0.5,1.5)
        ratioD.GetYaxis().SetNdivisions(502,0)
        ratioD.GetYaxis().SetLabelSize(0.05)
        ratioD.SetLineColor(2)
        ratioD.SetLineStyle(3)
        ratioD.SetLineWidth(2)  
        ratioU.SetLineColor(4)    
        ratioU.SetLineStyle(4)
        ratioU.SetLineWidth(2)

        fitRatioU = ratioU.Fit("pol2","S")
        ratioU.GetFunction("pol2").SetLineColor(4)
        fitRatioD = ratioD.Fit("pol2","S")
        ratioU.Draw("APSAME")
        ratioD.GetXaxis().SetTitle('BDT Output')
        ratioD.GetYaxis().SetTitle('Ratio')
        ratioD.GetYaxis().SetTitleSize(0.1)
        ratioD.GetYaxis().SetTitleOffset(0.2)
        fitRatioU.Draw("SAME")
        fitRatioD.Draw("SAME")
                    
        ratioD.Draw("SAME")

        m_one_line = ROOT.TLine(-1,1,1,1)
        m_one_line.SetLineStyle(7)
        m_one_line.SetLineColor(1)
        m_one_line.Draw("Same")

        c.Update()
        c.SaveAs("%s/%s_%s_%s.png" % (channel,systematic,channel,sample))
        #c.SaveAs("%s_%s_%s.png" % (systematic,channel,sample))
        c.Update()
        #c.Close()
