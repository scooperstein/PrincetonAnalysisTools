// 
// Fit regressed vs. non-regressed Mjj distributions.
// Used to make analysis note plots. Adpated from VBF script.
//
// Author: Nadezda Chernyavskaya
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
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


using namespace RooFit ;
using namespace std;


int mbb(){

	
//int JECR_unc_mbb(int type_int, int type_jec){
	gROOT->ProcessLine(".x ./setTDRStyle.C");


	int type_int = 0; //double  =0; single=1
	
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
	//TFile *file_down =  new TFile("/afs/cern.ch/work/n/nchernya/Hbb/SysUnc/JECR/Mbb/VBFHToBB_M-125_analysis"+type+"_trignom_v21_v21_123jets1_2csv_bdt_trigWt.root");
        //TFile *file_down = new TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/VHbbAnalysis/V21_Wlnu_June27_SR/output_mc.root");
        //TFile *file_down = new TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/VHbbAnalysis/V21_Wlnu_June28_CR-v2/output_allmc.root");
        TFile *file_down = new TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/VHbbAnalysis/V21_Wlnu_June29_SR_justSignal/output_mc.root");
        TTree *tree = (TTree*) file_down->Get("tree");
        TH1F *hist_mbb_reg = new TH1F("hist_mbb_reg","hist_mbb_reg",100,0,250);
        TH1F *hist_mbb = new TH1F("hist_mbb","hist_mbb",100,0,250);
	
        tree->Draw("H_mass>>hist_mbb_reg","(sampleIndex==-12501&&Pass_nominal==1&&(Vtype==2||Vtype==3))*weight"); 
        tree->Draw("H_mass_noreg>>hist_mbb","(sampleIndex==-12501&&Pass_nominal==1&&(Vtype==2||Vtype==3))*weight");    
        
        //TH1F *hist_mbb_reg = new TH1F("hist_mbb_reg","hist_mbb_reg",100,0,250);
        //TH1F *hist_mbb = new TH1F("hist_mbb","hist_mbb",100,0,250);

        //tree->Draw("Top1_mass_fromLepton_regPT_w4MET>>hist_mbb_reg","(sampleIndex>0&&sampleIndex!=6000&&sampleIndex!=6001&&sampleIndex!=6002&&controlSample==1&&(Vtype==2||Vtype==3))*weight"); 
       // tree->Draw("Top1_mass_fromLepton_w4MET>>hist_mbb","(sampleIndex>0&&sampleIndex!=6000&&sampleIndex!=6001&&sampleIndex!=6002&&controlSample==1&&(Vtype==2||Vtype==3))*weight");    

        TFile *fout = new TFile("hists.root","RECREATE");
        hist_mbb_reg->Write();
        hist_mbb->Write();

        cout<<hist_mbb_reg->Integral()<<endl;
        cout<<hist_mbb->Integral()<<endl;

        //hist_mbb_reg->Scale(1./hist_mbb_reg->Integral());	
	//hist_mbb->Scale(1./hist_mbb->Integral());	
/////////
//workspaces
	

	char nVar[50], xVar[50], Nb[50];
	sprintf(xVar,"mbb_reg");
	RooRealVar x(xVar,xVar,50,250);
	sprintf(xVar,"roofit_hist_mbb_reg");
	RooDataHist rh(xVar,xVar,RooArgList(x),hist_mbb_reg);
	sprintf(xVar,"yield_signal_VBF_reg");
	TString tmp;
	RooRealVar m("mean_reg","mean_reg",125,115,135);
	RooRealVar s("sigma_reg","sigma_reg",12,7,20);
	RooRealVar width("fwhm_reg","fwhm_reg",25,0,50);
	RooRealVar a("alpha_reg","alpha_reg",1,-10,10);
	RooFormulaVar mShift("mean_shifted_reg","mean_shifted_reg","@0",RooArgList(m));
	RooFormulaVar sShift("sigma_shifted_reg","sigma_shifted_reg","@0",RooArgList(s));
	RooRealVar n("exp_reg","exp_reg",1,0,100);
	RooRealVar b0("b0_reg","b0_reg",0.5,0.,1.);
	RooRealVar b1("b1_reg","b1_reg",0.5,0.,1.);
	RooRealVar b2("b2_reg","b2_reg",0.5,0.,1.);
	RooRealVar b3("b3_reg","b3_reg",0.5,0.,1.);
 /// Bkg part: Bernste[i].Data()n     
	RooBernstein bkg("signal_bkg_reg", "signal_bkg_reg" ,x,RooArgList(b0,b1,b2));
 /// Sig part: Crystal Ball   
	RooRealVar fsig("fsig_reg","fsig_reg",0.7,0.,1.);
	RooCBShape sig("signal_gauss_reg",  "signal_gauss_reg",x,mShift,sShift,a,n);
 //// Combined model
	RooAddPdf model("signal_model_reg", "signal_model_reg" ,RooArgList(sig,bkg),RooArgList(fsig));
 //// Fit 
 
	model.fitTo(rh,RooFit::SumW2Error(kFALSE));
//	model.fitTo(rh,RooFit::SumW2Error(kTRUE));


	b0.setConstant(kTRUE);
	b1.setConstant(kTRUE);
	b2.setConstant(kTRUE);
	b3.setConstant(kTRUE);
	a.setConstant(kTRUE);
	m.setConstant(kTRUE);
	s.setConstant(kTRUE);
	n.setConstant(kTRUE);
	fsig.setConstant(kTRUE);
	w.import(rh);
	w.import(model);
	w.import(width);
////////
//
//Mbb
//
char nVar2[50], xVar2[50];
	sprintf(xVar2,"mbb");
	RooRealVar x2(xVar2,xVar2,50,250);
	sprintf(xVar2,"roofit_hist_mbb");
	RooDataHist rh2(xVar2,xVar2,RooArgList(x2),hist_mbb);
	sprintf(xVar2,"yield_signal_VBF");
	RooRealVar m2("mean","mean",125,115,135);
	RooRealVar s2("sigma","sigma",12,7,20);
	RooRealVar width2("fwhm","fwhm",25,0,50);
	RooRealVar a2("alpha","alpha",1,-10,10);
	RooFormulaVar mShift2("mean_shifted","mean_shifted","@0",RooArgList(m2));
	RooFormulaVar sShift2("sigma_shifted","sigma_shifted","@0",RooArgList(s2));
	RooRealVar n2("exp","exp",1,0,100);
	RooRealVar b02("b0","b0",0.5,0.,1.);
	RooRealVar b12("b1","b1",0.5,0.,1.);
	RooRealVar b22("b2","b2",0.5,0.,1.);
	RooRealVar b32("b3","b3",0.5,0.,1.);
 /// Bkg part: Bernste[i].Data()n     
	RooBernstein bkg2("signal_bkg", "signal_bkg" ,x2,RooArgList(b02,b12,b22));
 /// Sig part: Crystal Ball   
	RooRealVar fsig2("fsig","fsig",0.7,0.,1.);
	RooCBShape sig2("signal_gauss",  "signal_gauss",x2,mShift2,sShift2,a2,n2);
 //// Combined model
	RooAddPdf model2("signal_model", "signal_model" ,RooArgList(sig2,bkg2),RooArgList(fsig2));
 //// Fit 
 
	model2.fitTo(rh2,RooFit::SumW2Error(kFALSE));


	b02.setConstant(kTRUE);
	b12.setConstant(kTRUE);
	b22.setConstant(kTRUE);
	b32.setConstant(kTRUE);
	a2.setConstant(kTRUE);
	m2.setConstant(kTRUE);
	s2.setConstant(kTRUE);
	n2.setConstant(kTRUE);
	fsig2.setConstant(kTRUE);

//////////////////////////////////////
 /// Draw 
	 // 			RooDraw(opts,canM,x,rhs[N],hVBF[N],hGF[N],hTOT[N],yVBF[N],yGF[N],model,bkg,fsig,width,mass,m,s,iS,S,C,Cp,archive)
 	TCanvas* c_roo = new TCanvas("RooFit","",800,800) ;
	RooPlot *frame_roo = x.frame(Name("Mbb")) ;
	
	
	RooDataHist hScale("tmp","tmp",RooArgList(x),hist_mbb_reg);
	hScale.plotOn(frame_roo,DrawOption("Psame"),LineWidth(2),LineColor(kBlue+1),FillColor(kBlue-9),MarkerColor(kBlue+1), FillStyle(1001)) ;
	model.plotOn(frame_roo);
	frame_roo->GetXaxis()->SetNdivisions(505);
	//frame_roo->GetXaxis()->SetTitle("M_{bb} (GeV)");
	frame_roo->GetXaxis()->SetTitle("Top Mass (GeV)");
	//frame_roo->GetYaxis()->SetTitle("1/N #times dN/dM_{bb}");
	frame_roo->GetYaxis()->SetTitle("Events / 2.5 GeV");
  	frame_roo->Draw() ;

	hist_mbb_reg->SetFillColor(kBlue-10);
	hist_mbb_reg->SetLineColor(kBlue+1);
	hist_mbb_reg->SetFillStyle(1001);
	hist_mbb_reg->SetMarkerColor(kBlue+1);
	hist_mbb_reg->SetMarkerStyle(20);

	THStack *hS = new THStack("hs","hs");
	hS->Add(hist_mbb_reg);
	hS->Draw("same,hist");
	hist_mbb->SetLineColor(kRed);
	hist_mbb->SetMarkerColor(kRed);
	hist_mbb->SetMarkerStyle(21);
	THStack *hS2 = new THStack("hs2","hs2");
	hS2->Add(hist_mbb);
	hS2->Draw("same,hist");
	frame_roo->Draw("same");
	gPad->RedrawAxis();


	frame_roo=x2.frame(Name("Mbb"));
	RooDataHist hScale2("tmp","tmp",RooArgList(x2),hist_mbb);
	hScale2.plotOn(frame_roo,DrawOption("Psame"),LineWidth(2),LineColor(kRed),MarkerColor(kRed), MarkerStyle(21), FillStyle(1001)) ;
	model2.plotOn(frame_roo,LineColor(kRed));
	frame_roo->Draw("same");




	gPad->RedrawAxis();

	TF1 *ftmp= (TF1*) model.asTF(RooArgList(x),RooArgList(fsig),RooArgSet(x));
	float y_0       = ftmp->GetMaximum();
	float x_0 = ftmp->GetMaximumX();
	float x_1       = ftmp->GetX(y_0/2.,50,x_0);
	float x_2 =ftmp->GetX(y_0/2.,x_0,250);
	float fwhm        = x_2-x_1;

	TF1 *ftmp2= (TF1*) model2.asTF(RooArgList(x2),RooArgList(fsig2),RooArgSet(x2));
	float y_02       = ftmp2->GetMaximum();
	float x_02 = ftmp2->GetMaximumX();
	float x_12       = ftmp2->GetX(y_02/2.,50,x_02);
	float x_22 =ftmp2->GetX(y_02/2.,x_02,250);
	float fwhm2        = x_22-x_12;

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
	pave.SetTextColor(kBlack);
	pave.AddText("M_{H} = 125 GeV");

	sprintf(cmd,"\%s selection",text_type.Data())	;
	//pave.AddText(cmd);

	TLegend *leg = new TLegend(0.6,0.65,1.-right-top*0.3333,0.8);
	leg->AddEntry(hist_mbb_reg,"Regressed","PL");
	leg->AddEntry(hist_mbb,"Raw","PL");
	leg->SetFillStyle(-1);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(top*0.5);
	leg->SetY1(leg->GetY2()-top*leg->GetNRows()*0.87);
	leg->Draw("same");

	TPaveText pave2(0.6,0.59,1.-right-top*0.5333,0.69,"NDC");
	pave2.SetTextAlign(11);
	pave2.SetFillStyle(-1);
	pave2.SetBorderSize(0);
	pave2.SetTextFont(42);
	pave2.SetTextSize(top*0.5);
	pave2.SetTextColor(kBlue+1);
	sprintf(cmd,"PEAK = \%.1f",m.getVal())	;
	pave2.AddText(cmd);
	sprintf(cmd,"FWHM = \%.1f",fwhm)	;
	pave2.AddText(cmd);
	pave2.Draw("same");
	TPaveText pave3(0.6,0.49,1.-right-top*0.5333,0.59,"NDC");
	pave3.SetTextAlign(11);
	pave3.SetFillStyle(-1);
	pave3.SetBorderSize(0);
	pave3.SetTextFont(42);
	pave3.SetTextSize(top*0.5);
	pave3.SetTextColor(kRed);
	sprintf(cmd,"PEAK = \%.1f",m2.getVal())	;
	pave3.AddText(cmd);
	sprintf(cmd,"FWHM = \%.1f",fwhm2)	;
	pave3.AddText(cmd);
	pave3.Draw("same");


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
	pCMS2.AddText("13 TeV");
	gPad->Update();
	pCMS1.Draw("same");
	pCMS2.Draw("same");
	gPad->Update();
	pave.Draw("same");
	pCMS1.Draw("same");
	pCMS12.Draw("same");
	pCMS2.Draw("same");

	c_roo->Print("mbb_reg_v21.pdf");



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
	TLegend *leg = new TLegend(0.77,0.6,0.92,0.9);
	leg->SetFillColor(0);
//	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.025);
	leg->SetHeader(jec+" variation");
//	leg->AddEntry(hist_mbb_JEnom,"central","L");
//	leg->AddEntry(hist_mbb_JEdown,jec+" down","L");
//	leg->AddEntry(hist_mbb_JEup,jec+" up","L");
	
//	leg->Draw("same");
//	c11->Print("plots/shape/v14_vbf"+type+"_"+jec+"_mbb_"+cats_names[cat_num]+corType[i]+".png");
*/

	return 0;

}

