
// Written by Stephane Cooperstein
//
// Creates workspaces used for limits
// Inheriting from Analysis Manager
//

#include "WorkspaceAnalysis.h"

#include "RooRealVar.h"

// initialize parameters
WorkspaceAnalysis::WorkspaceAnalysis(){
    if(debug>10) std::cout<<"Constructing WorkspaceAnalysis"<<std::endl;
}

// remove any 'new'
WorkspaceAnalysis::~WorkspaceAnalysis(){
    if(debug>10) std::cout<<"Deconstructing WorkspaceAnalysis"<<std::endl;
}

void WorkspaceAnalysis::InitAnalysis() {
    if(debug>1000) std::cout<<"Called InitAnalysis()"<<std::endl;
    // initialize histograms
    for(int i=0; i < (int)samples.size(); i++) {
        if(debug>1000) std::cout<<"Initializing histograms for sample: "<<samples[i].sampleName.c_str()<<std::endl;
        const char* sname = samples[i].sampleName.c_str();
        float nBinsX = *f["nBinsX"];
        float varMinX = *f["varMinX"];
        float varMaxX = *f["varMaxX"];
        float nBinsY = *f["nBinsY"];
        float varMinY = *f["varMinY"];
        float varMaxY = *f["varMaxY"];
        hists1D.push_back(new TH1F(sname, sname, nBinsX, varMinX, varMaxX) );
        hists2D.push_back(new TH2F(sname, sname, nBinsX, varMinX, varMaxX, nBinsY, varMinY, varMaxY) );
    }

    // create output file for histograms
    histout = new TFile("hists.root", "recreate");
}

bool WorkspaceAnalysis::Preselection() {
    // for now we won't impose any selections on top of the input tree from VHbbAnalyzer
    return true; 
}

bool WorkspaceAnalysis::Analyze() {
    // no selection imposed on top of VHbbAnalyzer
    if(debug>10000) std::cout<<"Called Analyze()"<<std::endl;
    return true;
}

void WorkspaceAnalysis::FinishEvent() {
    for(int i=0; i < (int)samples.size(); i++) {
        // fill histograms with current values of BDT, Hmass
        float BDT = *f["BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4"];
        float H_mass = *f["H_mass_f"];
        float H_mass_reg = *f["H_mass_reg"];
        int sNum = *in["sampleIndex"];        
        float weight = samples[i].intWeight;
        if(sNum != samples[i].sampleNum) continue; 
        ofile->cd();
        hists1D[i]->Fill(BDT, weight);
        if(debug>20000) std::cout<<"1D hist filled with BDT value: "<<BDT<<" and weight: "<<weight<<std::endl;
        hists2D[i]->Fill(BDT, H_mass_reg, weight);
        if(debug>20000) std::cout<<"2D hist filled with BDT value: "<<BDT<<", Mjj value: "<<H_mass_reg<<" and weight: "<<weight<<std::endl;
    }
    
}

void WorkspaceAnalysis::TermAnalysis() {

    // initialize all the RooHistPDF's
    RooWorkspace *WS = new RooWorkspace("workspace", "workspace");
    RooArgList *obs;
    float varMinX = *f["varMinX"];
    float varMaxX = *f["varMaxX"];
    float varMinY = *f["varMinY"];
    float varMaxY = *f["varMaxY"];
    WS->factory(Form("CMS_vhbb_BDT_Wln_8TeV[%f,%f]",varMinX,varMaxX));
    WS->factory(Form("CMS_vhbb_Mjj_Wln_8TeV[%f,%f]",varMinY,varMaxY));
  
    for(int i=0; i < (int)samples.size(); i++) {
        // N.B.: This will only work correctly if the samples list hasn't changed at all since calling InitAnalysis()
        SampleContainer *csample = &samples[i];
        RooDataHist *tmpRDH = new RooDataHist();
       if(debug>1000) std::cout<<"modelDimension set to: "<<int(*f["modelDimension"])<<::std::endl;
        switch(int(*f["modelDimension"]))
        {
        case 1: 
            obs = new RooArgList(*WS->var("CMS_vhbb_BDT_Wln_8TeV"));
            tmpRDH = new RooDataHist(csample->sampleName.c_str(), csample->sampleName.c_str(), *obs, hists1D[i]);
            break;
        case 2: 
            obs = new RooArgList(*WS->var("CMS_vhbb_BDT_Wln_8TeV"), *WS->var("CMS_vhbb_Mjj_Wln_8TeV"));
            tmpRDH = new RooDataHist(csample->sampleName.c_str(), csample->sampleName.c_str(), *obs, hists2D[i]);
            break;
        default:
            std::cerr<<"modelDimension must be set to 1 or 2! It is currently set to: "<<*f["modelDimension"]<<std::endl;
        }

        RooHistPdf *tmpRHP = new RooHistPdf(samples[i].sampleName.c_str(), samples[i].sampleName.c_str(), *obs, *tmpRDH);        
        if(debug>2000) std::cout<<"created RooHistPDF for sample: "<<samples[i].sampleName.c_str()<<std::endl;
        WS->import(*tmpRHP); 
        if(debug>2000) std::cout<<"imported RooHistPdf into workspace"<<std::endl;
        histout->cd();
        hists1D[i]->Write();
        hists2D[i]->Write();
    }
    // do fits

    // delete variables/etc
    ofile->cd();
    WS->writeToFile(Form("%s.root",outputTreeName.c_str()));
    ofile->Close();
}
