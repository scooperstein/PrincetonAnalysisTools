//
// Perform fits on the pT(W) distribution for backgrounds in CR's.
// used to correct for data/MC discrepancies in pT(W) distribution
// with full 2016 dataset.
//
// Author: Stephane Cooperstein
//

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;


void fitPtWCorrs() {

using namespace RooFit ;

float SF_TT = 0.84;
float SF_Wj0b = 1.01;
float SF_Wj1b = 1.49;
float SF_Wj2b = 1.94; 

TFile *ifile_ttWmn = TFile::Open("hists_ttWmn_to400.root", "r");
TH1F *hist_ttWmn_data = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_data_obs");
TH1F *hist_ttWmn_tt = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_TT");
TH1F *hist_ttWmn_wj2b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Wj2b");
TH1F *hist_ttWmn_wj1b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Wj1b");
TH1F *hist_ttWmn_wj0b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Wj0b");
TH1F *hist_ttWmn_stop = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_s_Top");
TH1F *hist_ttWmn_zj2b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Zj2b");
TH1F *hist_ttWmn_zj1b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Zj1b");
TH1F *hist_ttWmn_zj0b = (TH1F*) ifile_ttWmn->Get("BDT_ttWmn_Zj0b");

TFile *ifile_wlfWmn = TFile::Open("hists_wlfWmn_to400.root", "r");
TH1F *hist_wlfWmn_data = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_data_obs");
TH1F *hist_wlfWmn_tt = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_TT");
TH1F *hist_wlfWmn_wj2b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Wj2b");
TH1F *hist_wlfWmn_wj1b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Wj1b");
TH1F *hist_wlfWmn_wj0b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Wj0b");
TH1F *hist_wlfWmn_stop = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_s_Top");
TH1F *hist_wlfWmn_zj2b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Zj2b");
TH1F *hist_wlfWmn_zj1b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Zj1b");
TH1F *hist_wlfWmn_zj0b = (TH1F*) ifile_wlfWmn->Get("BDT_wlfWmn_Zj0b");

TFile *ifile_whfWmn = TFile::Open("hists_whfWmn_to400.root", "r");
TH1F *hist_whfWmn_data = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_data_obs");
TH1F *hist_whfWmn_tt = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_TT");
TH1F *hist_whfWmn_wj2b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Wj2b");
TH1F *hist_whfWmn_wj1b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Wj1b");
TH1F *hist_whfWmn_wj0b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Wj0b");
TH1F *hist_whfWmn_stop = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_s_Top");
TH1F *hist_whfWmn_zj2b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Zj2b");
TH1F *hist_whfWmn_zj1b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Zj1b");
TH1F *hist_whfWmn_zj0b = (TH1F*) ifile_whfWmn->Get("BDT_whfWmn_Zj0b");

hist_ttWmn_tt->Scale(SF_TT);
hist_ttWmn_wj0b->Scale(SF_Wj0b);
hist_ttWmn_wj1b->Scale(SF_Wj1b);
hist_ttWmn_wj2b->Scale(SF_Wj2b);

hist_wlfWmn_tt->Scale(SF_TT);
hist_wlfWmn_wj0b->Scale(SF_Wj0b);
hist_wlfWmn_wj1b->Scale(SF_Wj1b);
hist_wlfWmn_wj2b->Scale(SF_Wj2b);

hist_whfWmn_tt->Scale(SF_TT);
hist_whfWmn_wj0b->Scale(SF_Wj0b);
hist_whfWmn_wj1b->Scale(SF_Wj1b);
hist_whfWmn_wj2b->Scale(SF_Wj2b);

TH1F *hist_ttWmn_sub = hist_ttWmn_data->Clone();
hist_ttWmn_sub->Reset();
hist_ttWmn_sub->Add(hist_ttWmn_wj2b);
hist_ttWmn_sub->Add(hist_ttWmn_wj1b);
hist_ttWmn_sub->Add(hist_ttWmn_wj0b);
hist_ttWmn_sub->Add(hist_ttWmn_zj2b);
hist_ttWmn_sub->Add(hist_ttWmn_zj1b);
hist_ttWmn_sub->Add(hist_ttWmn_zj0b);
hist_ttWmn_sub->Add(hist_ttWmn_stop);

TH1F *hist_ttWmn_whf = hist_ttWmn_wj2b->Clone();
hist_ttWmn_whf->Add(hist_ttWmn_wj1b);
hist_ttWmn_whf->Add(hist_ttWmn_stop);

TH1F *hist_wlfWmn_sub = hist_wlfWmn_data->Clone();
hist_wlfWmn_sub->Reset();
hist_wlfWmn_sub->Add(hist_wlfWmn_wj2b);
hist_wlfWmn_sub->Add(hist_wlfWmn_wj1b);
hist_wlfWmn_sub->Add(hist_wlfWmn_tt);
hist_wlfWmn_sub->Add(hist_wlfWmn_zj2b);
hist_wlfWmn_sub->Add(hist_wlfWmn_zj1b);
hist_wlfWmn_sub->Add(hist_wlfWmn_zj0b);
hist_wlfWmn_sub->Add(hist_wlfWmn_stop);

TH1F *hist_wlfWmn_whf = hist_wlfWmn_wj2b->Clone();
hist_wlfWmn_whf->Add(hist_wlfWmn_wj1b);
hist_wlfWmn_whf->Add(hist_wlfWmn_stop);

TH1F *hist_whfWmn_sub = hist_whfWmn_data->Clone();
hist_whfWmn_sub->Reset();
hist_whfWmn_sub->Add(hist_whfWmn_wj0b);
hist_whfWmn_sub->Add(hist_whfWmn_tt);

TH1F *hist_whfWmn_whf = hist_whfWmn_wj2b->Clone();
hist_whfWmn_whf->Add(hist_whfWmn_wj1b);
hist_whfWmn_whf->Add(hist_whfWmn_stop);

std::cout<<"Total data: "<<hist_ttWmn_data->Integral()<<std::endl;;
std::cout<<"Total TT: "<<hist_ttWmn_tt->Integral()<<std::endl;;
std::cout<<"Total non-TT: "<<hist_ttWmn_sub->Integral()<<std::endl;;

std::cout<<"Total data: "<<hist_wlfWmn_data->Integral()<<std::endl;;
std::cout<<"Total WLF: "<<hist_wlfWmn_wj0b->Integral()<<std::endl;;
std::cout<<"Total non-WLF: "<<hist_wlfWmn_sub->Integral()<<std::endl;;

std::cout<<"Total data: "<<hist_whfWmn_data->Integral()<<std::endl;;
std::cout<<"Total WHF: "<<hist_whfWmn_whf->Integral()<<std::endl;;
std::cout<<"Total non-WHF: "<<hist_whfWmn_sub->Integral()<<std::endl;;

hist_ttWmn_data->Sumw2();
hist_wlfWmn_data->Sumw2();
hist_whfWmn_data->Sumw2();

//hist_ttWmn_data->Add(hist_ttWmn_sub, -1.0);
//hist_ttWmn_tt->Add(hist_ttWmn_sub);
//hist_wlfWmn_wlf->Add(hist_wlfWmn_sub);

//std::cout<<"Total data after subtraction: "<<hist_ttWmn_data->Integral()<<std::endl;;

//hist_ttWmn_data->Divide(hist_ttWmn_tt);

RooWorkspace w("ws");

RooRealVar vpt("vpt","vpt",100,400);
RooDataHist data_ttWmn("data_ttWmn","data_ttWmn",vpt,hist_ttWmn_data);
RooDataHist tt_ttWmn("tt_ttWmn","tt_ttWmn",vpt,hist_ttWmn_tt);
RooDataHist wj0b_ttWmn("wj0b_ttWmn","wj0b_ttWmn",vpt,hist_ttWmn_wj0b);
RooDataHist wj1b_ttWmn("wj1b_ttWmn","wj1b_ttWmn",vpt,hist_ttWmn_wj1b);
RooDataHist wj2b_ttWmn("wj2b_ttWmn","wj2b_ttWmn",vpt,hist_ttWmn_wj2b);
RooDataHist stop_ttWmn("stop_ttWmn","stop_ttWmn",vpt,hist_ttWmn_stop);
RooDataHist zj0b_ttWmn("zj0b_ttWmn","zj0b_ttWmn",vpt,hist_ttWmn_zj0b);
RooDataHist zj1b_ttWmn("zj1b_ttWmn","zj1b_ttWmn",vpt,hist_ttWmn_zj1b);
RooDataHist zj2b_ttWmn("zj2b_ttWmn","zj2b_ttWmn",vpt,hist_ttWmn_zj2b);
RooDataHist whf_ttWmn("whf_ttWmn","whf_ttWmn",vpt,hist_ttWmn_whf);

RooDataHist data_wlfWmn("data_wlfWmn","data_wlfWmn",vpt,hist_wlfWmn_data);
RooDataHist tt_wlfWmn("tt_wlfWmn","tt_wlfWmn",vpt,hist_wlfWmn_tt);
RooDataHist wj0b_wlfWmn("wj0b_wlfWmn","wj0b_wlfWmn",vpt,hist_wlfWmn_wj0b);
RooDataHist wj1b_wlfWmn("wj1b_wlfWmn","wj1b_wlfWmn",vpt,hist_wlfWmn_wj1b);
RooDataHist wj2b_wlfWmn("wj2b_wlfWmn","wj2b_wlfWmn",vpt,hist_wlfWmn_wj2b);
RooDataHist stop_wlfWmn("stop_wlfWmn","stop_wlfWmn",vpt,hist_wlfWmn_stop);
RooDataHist zj0b_wlfWmn("zj0b_wlfWmn","zj0b_wlfWmn",vpt,hist_wlfWmn_zj0b);
RooDataHist zj1b_wlfWmn("zj1b_wlfWmn","zj1b_wlfWmn",vpt,hist_wlfWmn_zj1b);
RooDataHist zj2b_wlfWmn("zj2b_wlfWmn","zj2b_wlfWmn",vpt,hist_wlfWmn_zj2b);
RooDataHist whf_wlfWmn("whf_wlfWmn","whf_wlfWmn",vpt,hist_wlfWmn_whf);

RooDataHist data_whfWmn("data_whfWmn","data_whfWmn",vpt,hist_whfWmn_data);
RooDataHist tt_whfWmn("tt_whfWmn","tt_whfWmn",vpt,hist_whfWmn_tt);
RooDataHist wj0b_whfWmn("wj0b_whfWmn","wj0b_whfWmn",vpt,hist_whfWmn_wj0b);
RooDataHist wj1b_whfWmn("wj1b_whfWmn","wj1b_whfWmn",vpt,hist_whfWmn_wj1b);
RooDataHist wj2b_whfWmn("wj2b_whfWmn","wj2b_whfWmn",vpt,hist_whfWmn_wj2b);
RooDataHist stop_whfWmn("stop_whfWmn","stop_whfWmn",vpt,hist_whfWmn_stop);
RooDataHist zj0b_whfWmn("zj0b_whfWmn","zj0b_whfWmn",vpt,hist_whfWmn_zj0b);
RooDataHist zj1b_whfWmn("zj1b_whfWmn","zj1b_whfWmn",vpt,hist_whfWmn_zj1b);
RooDataHist zj2b_whfWmn("zj2b_whfWmn","zj2b_whfWmn",vpt,hist_whfWmn_zj2b);
RooDataHist whf_whfWmn("whf_whfWmn","whf_whfWmn",vpt,hist_whfWmn_whf);

//RooRealVar a0_tt("a0_tt","a0_tt",0.);
//RooRealVar a0_wlf("a0_wlf","a0_wlf",0.);
//RooRealVar a0_whf("a0_whf","a0_whf",0.);
RooRealVar a0_tt("a0_tt","a0_tt",-0.00086,-0.01,0.);
RooRealVar a0_wlf("a0_wlf","a0_wlf",-0.00075,-0.01,0.0);
RooRealVar a0_whf("a0_whf","a0_whf",-0.00075,-0.01,0.0);
//RooRealVar a0_tt("a0_tt","a0_tt",-0.00001,-0.01,0.01);
//RooRealVar a1_tt("a1_tt","a1_tt",-0.00001,-0.01,0.01);
//RooRealVar a0_wlf("a0_wlf","a0_wlf",-0.0007,-0.01,0.01);
RooRealVar a1_wlf("a1_wlf","a1_wlf",0.0002,-10000,10000);
//RooRealVar a0_whf("a0_whf","a0_whf",-0.009,-10000,10000);
RooRealVar a1_whf("a1_whf","a1_whf",0.0002,-10000,10000);

RooArgList arglist_tt(a0_tt);
//RooArgList arglist_tt(a0_tt,a1_tt);
RooArgList arglist_wlf(a0_wlf);
//RooArgList arglist_wlf(a0_wlf,a1_wlf);
RooArgList arglist_whf(a0_whf);
//RooArgList arglist(a0,a1,a2,a3,a4,a5,a6,a7,a8);
//arglist.add(a9);
//arglist.add(a10);
//arglist.add(a11);
//arglist.add(a12);
//arglist.add(a13);
RooPolynomial poly_tt("poly_tt","poly_tt",vpt,arglist_tt);
RooPolynomial poly_wlf("poly_wlf","poly_wlf",vpt,arglist_wlf);
RooPolynomial poly_whf("poly_whf","poly_whf",vpt,arglist_whf);


RooHistPdf ttpdf_ttWmn("ttpdf_ttWmn","ttpdf_ttWmn",vpt,tt_ttWmn);
RooProdPdf model_tt_ttWmn("model_tt_ttWmn","model_tt_ttWmn",ttpdf_ttWmn,poly_tt);
RooHistPdf wlfpdf_ttWmn("wlfpdf_ttWmn","wlfpdf_ttWmn",vpt,wj0b_ttWmn);
RooProdPdf model_wlf_ttWmn("model_wlf_ttWmn","model_wlf_ttWmn",wlfpdf_ttWmn,poly_wlf);
RooHistPdf whfpdf_ttWmn("whfpdf_ttWmn","whfpdf_ttWmn",vpt,whf_ttWmn);
RooProdPdf model_whf_ttWmn("model_whf_ttWmn","model_whf_ttWmn",whfpdf_ttWmn,poly_whf);
RooRealVar ttfrac_ttWmn("ttfrac_ttWmn","ttfrac_ttWmn",(hist_ttWmn_tt->Integral()/(hist_ttWmn_tt->Integral()+hist_ttWmn_sub->Integral())));
RooRealVar wlffrac_ttWmn("wlffrac_ttWmn","wlffrac_ttWmn",(hist_ttWmn_wj0b->Integral()/(hist_ttWmn_tt->Integral()+hist_ttWmn_sub->Integral())));
//RooAddPdf model_ttWmn("model_ttWmn","model_ttWmn",RooArgList(model_wlf_ttWmn,model_tt_ttWmn),ttfrac_ttWmn);
RooAddPdf model_ttWmn("model_ttWmn","model_ttWmn",RooArgList(model_tt_ttWmn,model_wlf_ttWmn,model_whf_ttWmn),RooArgList(ttfrac_ttWmn,wlffrac_ttWmn));

RooHistPdf ttpdf_wlfWmn("ttpdf_wlfWmn","ttpdf_wlfWmn",vpt,tt_wlfWmn);
RooProdPdf model_tt_wlfWmn("model_tt_wlfWmn","model_tt_wlfWmn",ttpdf_wlfWmn,poly_tt);
RooHistPdf wlfpdf_wlfWmn("wlfpdf_wlfWmn","wlfpdf_wlfWmn",vpt,wj0b_wlfWmn);
RooProdPdf model_wlf_wlfWmn("model_wlf_wlfWmn","model_wlf_wlfWmn",wlfpdf_wlfWmn,poly_wlf);
RooHistPdf whfpdf_wlfWmn("whfpdf_wlfWmn","whfpdf_wlfWmn",vpt,whf_wlfWmn);
RooProdPdf model_whf_wlfWmn("model_whf_wlfWmn","model_whf_wlfWmn",whfpdf_wlfWmn,poly_whf);
RooRealVar wlffrac_wlfWmn("wlffrac_wlfWmn","wlffrac_wlfWmn",(hist_wlfWmn_wj0b->Integral()/(hist_wlfWmn_wj0b->Integral()+hist_wlfWmn_sub->Integral())));
RooRealVar ttfrac_wlfWmn("ttfrac_wlfWmn","ttfrac_wlfWmn",(hist_wlfWmn_tt->Integral()/(hist_wlfWmn_wj0b->Integral()+hist_wlfWmn_sub->Integral())));
//RooAddPdf model_wlfWmn("model_wlfWmn","model_wlfWmn",RooArgList(model_wlf_wlfWmn,model_tt_wlfWmn),wlffrac_wlfWmn);
RooAddPdf model_wlfWmn("model_wlfWmn","model_wlfWmn",RooArgList(model_wlf_wlfWmn,model_tt_wlfWmn,model_whf_wlfWmn),RooArgList(wlffrac_wlfWmn,ttfrac_wlfWmn));

RooHistPdf ttpdf_whfWmn("ttpdf_whfWmn","ttpdf_whfWmn",vpt,tt_whfWmn);
RooProdPdf model_tt_whfWmn("model_tt_whfWmn","model_tt_whfWmn",ttpdf_whfWmn,poly_tt);
RooHistPdf wlfpdf_whfWmn("wlfpdf_whfWmn","wlfpdf_whfWmn",vpt,wj0b_whfWmn);
RooProdPdf model_wlf_whfWmn("model_wlf_whfWmn","model_wlf_whfWmn",wlfpdf_whfWmn,poly_wlf);
RooHistPdf whfpdf_whfWmn("whfpdf_whfWmn","whfpdf_whfWmn",vpt,whf_whfWmn);
RooProdPdf model_whf_whfWmn("model_whf_whfWmn","model_whf_whfWmn",whfpdf_whfWmn,poly_whf);
RooRealVar wlffrac_whfWmn("wlffrac_whfWmn","wlffrac_whfWmn",(hist_whfWmn_wj0b->Integral()/(hist_whfWmn_whf->Integral()+hist_whfWmn_sub->Integral())));
RooRealVar ttfrac_whfWmn("ttfrac_whfWmn","ttfrac_whfWmn",(hist_whfWmn_tt->Integral()/(hist_whfWmn_whf->Integral()+hist_whfWmn_sub->Integral())));
//RooAddPdf model_whfWmn("model_whfWmn","model_whfWmn",RooArgList(model_wlf_whfWmn,model_tt_whfWmn),wlffrac_whfWmn);
RooAddPdf model_whfWmn("model_whfWmn","model_whfWmn",RooArgList(model_wlf_whfWmn,model_tt_whfWmn,model_whf_whfWmn),RooArgList(wlffrac_whfWmn,ttfrac_whfWmn));
//RooRealVar n1("n1","n1",1,0,1000);
//RooRealVar n2("n2","n2",1,0,1000);
//RooExtendPdf model_ttWmn2("model_ttWmn2","model_ttWmn2",model_ttWmn,n1);
//RooExtendPdf model_wlfWmn2("model_wlfWmn2","model_wlfWmn2",model_wlfWmn,n2);

//RooHistPdf wlfpdf("wlfpdf","wlfpdf",vpt,wj0b_wlfWmn);
//RooProdPdf model_wlf("model_wlf","model_wlf",wlfpdf,poly_wlf);

//prod.chi2FitTo(data);
//prod.fitTo(data);
//model_ttWmn->chi2FitTo(data_ttWmn);
//model_tt_ttWmn->fitTo(data_ttWmn);

RooCategory sample("sample","sample") ;
sample.defineType("ttWmn") ;
sample.defineType("wlfWmn") ;
sample.defineType("whfWmn") ;

RooDataSet *data = model_ttWmn.generate(RooArgSet(vpt),1000); 
RooDataSet *data2 = model_wlfWmn.generate(RooArgSet(vpt),1000); 
RooDataSet *data3 = model_whfWmn.generate(RooArgSet(vpt),1000); 

//RooDataSet combData("combData","combined data",vpt,RooFit::Index(sample),RooFit::Import("ttWmn",*data),RooFit::Import("wlfWmn",*data2),RooFit::Import("whfWmn",*data3)) ;
RooDataSet combData("combData","combined data",RooArgSet(vpt,sample));

//RooRealVar y("y","y",0,100000);
for (int i=1; i<hist_ttWmn_data->GetNbinsX()+1; i++) {
    vpt = hist_ttWmn_data->GetBinLowEdge(i) + 0.5*hist_ttWmn_data->GetBinWidth(i);
    float val = hist_ttWmn_data->GetBinContent(i);
    sample.setLabel("ttWmn");
//    combData.add(RooArgSet(vpt,sample),val);
    cout<<i<<":, "<<vpt<<", "<<val<<std::endl;
    for (int j=0; j<val; j++) {
        combData.add(RooArgSet(vpt,sample));
    }
}
for (int i=1; i<hist_wlfWmn_data->GetNbinsX()+1; i++) {
    vpt = hist_wlfWmn_data->GetBinLowEdge(i) + 0.5*hist_wlfWmn_data->GetBinWidth(i);
    float val = hist_wlfWmn_data->GetBinContent(i);
    sample.setLabel("wlfWmn");
    //combData.add(RooArgSet(vpt,sample));
    //combData.add(RooArgSet(vpt,sample),val);
    cout<<i<<":, "<<vpt<<", "<<val<<std::endl;
    for (int j=0; j<val; j++) {
        combData.add(RooArgSet(vpt,sample));
    }
}
for (int i=1; i<hist_whfWmn_data->GetNbinsX()+1; i++) {
    vpt = hist_whfWmn_data->GetBinLowEdge(i) + 0.5*hist_whfWmn_data->GetBinWidth(i);
    float val = hist_whfWmn_data->GetBinContent(i);
    sample.setLabel("whfWmn");
    //combData.add(RooArgSet(vpt,sample));
    //combData.add(RooArgSet(vpt,sample),val);
    cout<<i<<":, "<<vpt<<", "<<val<<std::endl;
    for (int j=0; j<val; j++) {
        combData.add(RooArgSet(vpt,sample));
    }
}

RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

simPdf.addPdf(model_ttWmn,"ttWmn") ;
simPdf.addPdf(model_wlfWmn,"wlfWmn"); 
simPdf.addPdf(model_whfWmn,"whfWmn"); 
//simPdf.addPdf(model_ttWmn2,"ttWmn") ;
//simPdf.addPdf(model_wlfWmn2,"wlfWmn"); 

RooDataHist dh("dh","dh",*combData.get(),combData) ; 
RooChi2Var chi2("chi2","chi2",simPdf,dh);
RooMinuit m2(chi2);
m2.migrad();
//m2.minos();

//std::cout<<"simPDF.canBeExtend() = "<<simPdf.canBeExtended()<<std::endl;
//simPdf.chi2FitTo(combData) ;
//simPdf.fitTo(combData) ;


//RooChi2Var chi2_var_A("","", model_ttWmn, *data_ttWmn);
//RooMinuit m2_var_A(chi2_var_A);
//m2_var_A.migrad();

//poly_tt.fitTo(data);
//poly_tt.fitTo(data,SumW2Error(kFALSE));

//RooPlot* frame = vpt.frame();
//data_ttWmn.plotOn(frame);
//poly_tt.plotOn(frame);
//poly_wlf.plotOn(frame,LineColor(kGreen));
//model_ttWmn.plotOn(frame,LineColor(kRed));
//model_tt_ttWmn.plotOn(frame,LineColor(kBlue));
//model_wlf_ttWmn.plotOn(frame,LineColor(kGreen));

//combData.plotOn(frame,Cut("sample==sample::ttWmn")) ;
//simPdf.plotOn(frame,Slice(sample,"ttWmn"),ProjWData(sample,combData)) ;
//simPdf.plotOn(frame,Slice(sample,"ttWmn"),Components("model_ttWm"),ProjWData(sample,combData),LineStyle(kDashed)) ;

RooPlot* frame1 = vpt.frame(Bins(60),Title("TT CR")) ;
combData.plotOn(frame1,Cut("sample==sample::ttWmn")) ;
simPdf.plotOn(frame1,Slice(sample,"ttWmn"),ProjWData(sample,combData)) ;
simPdf.plotOn(frame1,Slice(sample,"ttWmn"),Components("model_tt_ttWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kBlue)) ;
simPdf.plotOn(frame1,Slice(sample,"ttWmn"),Components("model_wlf_ttWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen)) ;
simPdf.plotOn(frame1,Slice(sample,"ttWmn"),Components("model_whf_ttWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kViolet)) ;
simPdf.plotOn(frame1,Slice(sample,"ttWmn"),Components("model_ttWmn"),ProjWData(sample,combData),LineStyle(kDashed)) ;

RooPlot* frame2 = vpt.frame(Bins(60),Title("WLF CR")) ;
combData.plotOn(frame2,Cut("sample==sample::wlfWmn")) ;
simPdf.plotOn(frame2,Slice(sample,"wlfWmn"),ProjWData(sample,combData)) ;
simPdf.plotOn(frame2,Slice(sample,"wlfWmn"),Components("model_tt_wlfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kBlue)) ;
simPdf.plotOn(frame2,Slice(sample,"wlfWmn"),Components("model_wlf_wlfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen)) ;
simPdf.plotOn(frame2,Slice(sample,"wlfWmn"),Components("model_whf_wlfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kViolet)) ;
simPdf.plotOn(frame2,Slice(sample,"wlfWmn"),Components("model_wlfWmn"),ProjWData(sample,combData),LineStyle(kDashed));

RooPlot* frame3 = vpt.frame(Bins(60),Title("WHF CR")) ;
combData.plotOn(frame3,Cut("sample==sample::whfWmn")) ;
simPdf.plotOn(frame3,Slice(sample,"whfWmn"),ProjWData(sample,combData)) ;
simPdf.plotOn(frame3,Slice(sample,"whfWmn"),Components("model_tt_whfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kBlue)) ;
simPdf.plotOn(frame3,Slice(sample,"whfWmn"),Components("model_wlf_whfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen)) ;
simPdf.plotOn(frame3,Slice(sample,"whfWmn"),Components("model_whf_whfWmn"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kViolet)) ;
simPdf.plotOn(frame3,Slice(sample,"whfWmn"),Components("model_whfWmn"),ProjWData(sample,combData),LineStyle(kDashed));

//TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",600,600) ;
TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",1600,800) ;
c->Divide(3) ;
c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

//frame3->Draw();

std::cout<<"TT chi2/ndof = "<<frame1.chiSquare()<<std::endl;;
std::cout<<"WHF chi2/ndof = "<<frame3.chiSquare()<<std::endl;;
std::cout<<"WLF chi2/ndof = "<<frame2.chiSquare()<<std::endl;;

frame3.Print();
//ttpdf.plotOn(frame);
//nb.plotOn(frame);
//expo.plotOn(frame);
//bkgmodel.plotOn(frame);
//CBall.plotOn(frame);
//std::cout<<"chi2/ndof = "<<frame.chiSquare()<<std::endl;;
//frame.Draw();

//RooPlot* frame2 = vpt.frame();
//data_wlfWmn.plotOn(frame2);
//poly_tt.plotOn(frame2);
//poly_wlf.plotOn(frame2,LineColor(kGreen));
//model_wlfWmn.plotOn(frame2,LineColor(kRed));
//std::cout<<"chi2/ndof = "<<frame2.chiSquare()<<std::endl;;
//frame2.Draw();

std::cout<<hist_whfWmn_data->Integral()<<std::endl;
std::cout<<ttfrac_ttWmn->getVal()<<std::endl;;
std::cout<<(1-wlffrac_whfWmn->getVal()-ttfrac_whfWmn->getVal())<<std::endl;;
std::cout<<wlffrac_wlfWmn->getVal()<<std::endl;;
//std::cout<<(hist_wlfWmn_wj0b->Integral()/(hist_wlfWmn_wj0b->Integral()+hist_wlfWmn_sub->Integral()))<<std::endl;;

}
