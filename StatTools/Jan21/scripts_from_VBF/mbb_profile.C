//
// Make nice profile plots of analysis variables. Used for analysis note plots.
// Script adapted from VBF.
//
// Author: Nadezda Chernyavskaya
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "TPaveText.h"
#include "TVirtualPad.h"
#include "TProfile.h"


using namespace RooFit ;
using namespace std;


int main(int argc, char* argv[]){
//void mbb_profile(char* var1, char* var2, int nbinsx, float xlow1, float xhigh1, int nbinsy, float xlow2, float xhigh2, char *xlabel, char* ylabel){

        char* var1 = argv[1];
        char *var2 = argv[2];
        int nbinsx = atoi(argv[3]);
        float xlow1 = atof(argv[4]);
        float xhigh1 = atof(argv[5]);
        float nbinsy = atoi(argv[6]);
        float xlow2 = atof(argv[7]);
        float xhigh2 = atof(argv[8]);
        char* xlabel = argv[9];
        char* ylabel = argv[10];

	//int type_int = atoi(argv[1]);
	int type_int = 100000;
	
//int JECR_unc_mbb(int type_int, int type_jec){
	gROOT->ProcessLine(".x ./setTDRStyle.C");


//	int type_int; //double  =0; single=1
	
	TString type;
	if (type_int==0) type = "_double";
	if (type_int==1) type = "_single";
	TString text_type;
	if (type_int==0) text_type = "DoubleB";
	if (type_int==1) text_type = "SingleB";
	TString jec="";

	int NCATS;
	Float_t cats[10]; 
	Float_t Nevt_up[10];
	Float_t Nevt_nom[10];
	Float_t Nevt_down[10];
	Float_t SysUnc[10];
	TString cats_names[10];
//	int NCATS_single = 0;
//	if (type_int==0) NCATS_single= 4;
	if (type_int==0) NCATS=3+1;
	else if (type_int==1) NCATS=4+1;
	if (type_int==1) {
		Float_t cats_new[10]={-1.,0.28  , 0.72  , 0.87 , 0.93    , 0.9825};
		memcpy(cats, cats_new, sizeof cats_new);
	}
	if (type_int==0) {
		Float_t cats_new[10]={-1.,0.36 ,  0.76  , 0.89  ,0.9725};
		memcpy(cats, cats_new, sizeof cats_new);
	}
	for (int i=0;i<NCATS;i++){
		TString tmp;
		tmp.Form("%1d",i);
		cats_names[i] = "CAT";
		cats_names[i].Append(tmp);
	} 
//		Float_t cats[10] = {-1.,0.12,0.44,0.57,0.6725};

	Float_t lumi = 2320.;
	Float_t xsec[2] =  {2.20,25.69};

	const int num_ss = 2;	

//	TString s_names[num_ss] = {"VBFHToBB_M-125_13TeV_powheg","GF"};
	TString s_names[num_ss] = {"vbf_jecr_v21","gf_jecr_v21"};
//	TString s_names[num_ss] = {"vbf_76"};
	TString tex_s_names[num_ss] = {"VBF, m(H) = 125 GeV", "GF, m(H) = 125 GeV"};
		

//up down variabtion loop


RooWorkspace w("w","workspace");
	char cmd[50];
	char cmd2[50];
	char cmd3[50];








	TCanvas *c3 = new TCanvas();
	//TFile *file_down =  new TFile("BTagCSV_analysis"+type+"_trignom_v21_v21_123jets1_2csv_bdt_trigWt.root");
        //TFile *file_down = new TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/VHbbAnalysis/V21_Wlnu_June28_CR-v2/output_data.root");
        //TFile *file_down = TFile::Open("root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_CR_Nov29/output_data.root");
        TFile *file_down = TFile::Open("root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_April6_data/output_data.root");
        //TFile *file_down = TFile::Open("root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29/output_signal.root");
        //cout<<"here"<<endl;
	//TProfile *hist_mbb_reg= (TProfile*)file_down->Get("hprof");
        TTree *tree_data = (TTree*) file_down->Get("tree");
        TH2F *hdata = new TH2F("hdata","hdata",nbinsx,xlow1,xhigh1,nbinsy,xlow2,xhigh2);
        tree_data->Draw(Form("%s:%s>>hdata",var2,var1),"(sampleIndex==0&&controlSample==1&&(Vtype==2||Vtype==3))*weight"); 
        //tree_data->Draw(Form("%s:%s>>hdata",var2,var1),"(run<=276811&&sampleIndex==0&&controlSample==3&&(Vtype==2||Vtype==3))*weight"); 
        //tree_data->Draw(Form("%s:%s>>hdata",var2,var1),"(selLeptons_relIso_0<0.06&&H_mass>90&&H_mass<150&&cutFlow>=10&&(isWmunu||isWenu)&&V_pt>100&&(Vtype==2||Vtype==3))*weight"); 
        //tree_data->Draw(Form("%s:%s>>hdata",var2,var1),"((sampleIndex!=0||run<=276811)&&selLeptons_relIso_0<0.06&&H_mass>90&&H_mass<150&&cutFlow>=10&&(isWmunu||isWenu)&&V_pt>100&&(Vtype==2||Vtype==3))*weight"); 
        //tree_data->Draw(Form("%s:%s>>hdata",var2,var1),"(run<=276811&&sampleIndex==0&&controlSample==2&&(H_mass<90||H_mass>150)&&(Vtype==2||Vtype==3))*weight"); 
        cout<<hdata->Integral()<<endl;
        TProfile *hist_mbb_reg = hdata->ProfileX("hist_mbb_reg");
        //hist_mbb_reg->SetErrorOption("s");
        cout<<"here"<<endl;
        	
        hist_mbb_reg->SetLineColor(1);
	hist_mbb_reg->SetMarkerColor(1);
	hist_mbb_reg->SetMarkerStyle(20);
	hist_mbb_reg->SetStats(kFALSE);
/////////

	//TFile *file_qcd =  new TFile("QCD"+type+".root");
        //TFile *file_mc = new TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/VHbbAnalysis/V21_Wlnu_June28_CR-v2/output_allmc.root");
        //TFile *file_mc = TFile::Open("root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_CR_Nov29/output_allmc.root");
        //TFile *file_mc = TFile::Open("root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SR_Nov29/output_allmc.root");
	//TProfile *hist_mbb_reg_qcd= (TProfile*)file_qcd->Get("hprof");
        //TTree *tree = (TTree*) file_mc->Get("tree");
        TChain *tree = new TChain("tree");
        tree->Add("/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_April6_v2/haddjobs/sum_*_3.root");
        TH2F *hmc= new TH2F("hmc","hmc",nbinsx,xlow1,xhigh1,nbinsy,xlow2,xhigh2);
        tree->Draw(Form("%s:%s>>hmc",var2,var1),"(sampleIndex>0&&sampleIndex!=6000&&sampleIndex!=6001&&sampleIndex!=6002&&controlSample==1&&(Vtype==2||Vtype==3))*weight*bTagWeightMoriondCMVA*VPtCorrFactorSplit3");    
        //tree->Draw(Form("%s:%s>>hmc",var2,var1),"(sampleIndex>0&&sampleIndex!=6000&&sampleIndex!=6001&&sampleIndex!=6002&&controlSample==3&&(Vtype==2||Vtype==3))*weight");    
        //tree->Draw(Form("%s:%s>>hmc",var2,var1),"((sampleIndex!=0||run<=276811)&&selLeptons_relIso_0<0.06&&H_mass>90&&H_mass<150&&cutFlow>=10&&(isWmunu||isWenu)&&V_pt>100&&(Vtype==2||Vtype==3))*weight");    
        //tree->Draw(Form("%s:%s>>hmc",var2,var1),"(selLeptons_relIso_0<0.06&&H_mass>90&&H_mass<150&&cutFlow>=10&&(isWmunu||isWenu)&&V_pt>100&&(Vtype==2||Vtype==3))*weight*bTagWeightMoriondCMVA*VPtCorrFactorSplit3");    
        //tree->Draw(Form("%s:%s>>hmc",var2,var1),"(sampleIndex>0&&sampleIndex!=6000&&sampleIndex!=6001&&sampleIndex!=6002&&controlSample==2&&(H_mass<90||H_mass>150)&&(Vtype==2||Vtype==3))*weight");    
        cout<<hmc->Integral()<<endl;
        TProfile *hist_mbb_reg_qcd = hmc->ProfileX("hist_mbb_reg_qcd");
	hist_mbb_reg_qcd->SetStats(kFALSE);
        //hist_mbb_reg_qcd->SetErrorOption("s");
        cout<<"here"<<endl;
	//if (type_int==1) hist_mbb_reg_qcd->Rebin(10);
	//if (type_int==0) hist_mbb_reg_qcd->Rebin(4);
        cout<<hist_mbb_reg->Integral()<<endl;
        cout<<hist_mbb_reg_qcd->Integral()<<endl;

 	TCanvas* c_roo = new TCanvas("Fit","",800,800) ;
	TH1F *frame2 = new TH1F("frame2","",nbinsx,xlow1,xhigh1);
//	frame2->SetMinimum(0.);
	frame2->GetXaxis()->SetNdivisions(505);
	frame2->GetYaxis()->SetRangeUser(xlow2,xhigh2);
   frame2->SetStats(0);
	frame2->SetXTitle(xlabel);
	frame2->SetYTitle(ylabel);
	frame2->Draw();

        //hist_mbb_reg->SetFillStyle(1001);
        //hist_mbb_reg->SetFillColorAlpha(kRed,0.7);
        hist_mbb_reg->SetFillStyle(3001);
        hist_mbb_reg->SetMarkerStyle(1);

        //hist_mbb_reg_qcd->SetFillStyle(1001);
        //hist_mbb_reg_qcd->SetFillColorAlpha(kBlue,0.7);
        hist_mbb_reg_qcd->SetFillStyle(3001);
        hist_mbb_reg_qcd->SetFillColor(kBlue);
        hist_mbb_reg_qcd->SetMarkerStyle(1);

	hist_mbb_reg_qcd->Draw("e2same");
	//hist_mbb_reg->Draw("e2same");
	hist_mbb_reg->Draw("Psame");
	



	gPad->RedrawAxis();


	float right   = gStyle->GetPadRightMargin();
	float top = gStyle->GetPadTopMargin();
	float left =  gStyle->GetPadLeftMargin();
	float bottom = gStyle->GetPadBottomMargin();




	TPaveText pave(0.6,0.8,1.-right-top*0.5333,1.-top*1.666,"NDC");
	pave.SetTextAlign(11);
	pave.SetFillStyle(-1);
	pave.SetBorderSize(0);
	pave.SetTextFont(42);
	pave.SetTextSize(top*0.5);
	pave.SetTextColor(kBlue+1);

	//sprintf(cmd,"WHbb TT CS",text_type.Data())	;
	//sprintf(cmd,"WHbb W+HF CS",text_type.Data())	;
	//sprintf(cmd,"WHbb W+LF CS",text_type.Data())	;
	sprintf(cmd,"WHbb SR",text_type.Data())	;
	//sprintf(cmd,"\%s selection",text_type.Data())	;
	pave.AddText(cmd);

	pave.Draw("same");


// CMS info
	float left2 = gStyle->GetPadLeftMargin();
	float right2 = gStyle->GetPadRightMargin();
	float top2 = gStyle->GetPadTopMargin();
	float bottom2 = gStyle->GetPadBottomMargin();




	TPaveText pCMS1(left2,1.-top2,0.4,1.,"NDC");
	pCMS1.SetTextFont(62);
	pCMS1.SetTextSize(top2*0.75);
	pCMS1.SetTextAlign(12);
	pCMS1.SetFillStyle(-1);
	pCMS1.SetBorderSize(0);
	pCMS1.AddText("CMS");
	TPaveText pCMS12(left2+0.1,1.-top2*1.1,0.6,1.,"NDC");
	pCMS12.SetTextFont(52);
	pCMS12.SetTextSize(top2*0.75);
	pCMS12.SetTextAlign(12);
	pCMS12.SetFillStyle(-1);
	pCMS12.SetBorderSize(0);
	pCMS12.AddText("Preliminary");
	TPaveText pCMS2(0.5,1.-top2,1.-right2*0.5,1.,"NDC");
	pCMS2.SetTextFont(42);
	pCMS2.SetTextSize(top2*0.75);
	pCMS2.SetTextAlign(32);
	pCMS2.SetFillStyle(-1);
	pCMS2.SetBorderSize(0);
	pCMS2.AddText("L = 35.9 fb^{-1} (13 TeV)");
	gPad->Update();
	pCMS1.Draw("same");
	pCMS2.Draw("same");
	gPad->Update();
	pave.Draw("same");
	pCMS1.Draw("same");
	pCMS12.Draw("same");
	pCMS2.Draw("same");


	TLegend *leg = new TLegend(0.65,0.75,0.82,0.82);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.025);
	//leg->AddEntry(hist_mbb_reg,"Signal MC","F");
	leg->AddEntry(hist_mbb_reg,"Data","P");
	//leg->AddEntry(hist_mbb_reg_qcd,"Bkg. MC","F");
	leg->AddEntry(hist_mbb_reg_qcd,"MC","F");
	leg->Draw("same");


	c_roo->Print(Form("profile_plots/profile_%s_vs_%s.pdf",var1,var2));
	c_roo->Print(Form("profile_plots/profile_%s_vs_%s.png",var1,var2));
        c_roo->Close();
        file_down->Close();
        //file_mc->Close();

/*	
	TCanvas *c11 = new TCanvas();
	c11->SetBottomMargin(.15);
	c11->SetRightMargin(.25);
	c11->cd();
	TH1F *frame3 = new TH1F("frame3","",15,50.,200.);
	frame3->SetMinimum(0.);
   frame3->SetMaximum(hist_mbb_JEdown->GetMaximum()*1.6);
   frame3->SetStats(0);
	frame3->SetYTitle("Events");
	frame3->SetXTitle(hist_mbb_JEdown->GetXaxis()->GetTitle());
	frame3->GetYaxis()->SetNdivisions(505);
	frame3->GetXaxis()->SetLabelSize(0.0);
  	frame3->GetXaxis()->SetTitleSize(0.05);
  	frame3->GetXaxis()->SetLabelSize(0.04);
	frame3->Draw();
	TLatex* tex = new TLatex(0.75,0.95,"13 TeV");
   tex->SetNDC();
	tex->SetTextAlign(35);
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   TLatex *tex1 = new TLatex(0.17,0.95,"CMS");
   tex1->SetNDC();
   tex1->SetTextAlign(20);
   tex1->SetTextFont(61);
   tex1->SetTextSize(0.04);
   tex1->SetLineWidth(2);
   TLatex* tex2 = new TLatex(0.25,0.89,"Work in progress");
   tex2->SetNDC();
   tex2->SetTextAlign(20);
   tex2->SetTextFont(52);
   tex2->SetTextSize(0.035);
  	tex2->SetLineWidth(2);	
	TLatex* tex_file = new TLatex(0.36,0.95,text_type);
   tex_file->SetNDC();
	tex_file->SetTextAlign(35);
   tex_file->SetTextFont(42);
   tex_file->SetTextSize(0.04);
   tex_file->SetLineWidth(2);	
	tex->Draw();
	tex1->Draw();
	tex2->Draw();
//	tex_file->Draw();
	hist_mbb_JEdown->SetLineColor(8);
	hist_mbb_JEdown->SetLineWidth(2);
	hist_mbb_JEdown->Draw("PEsame");
//	c11->Print("plots/shape/v14_vbf"+type+"_"+jec+"_mbb_"+cats_names[cat_num]+corType[i]+".png");
*/

	//return 0;

}

