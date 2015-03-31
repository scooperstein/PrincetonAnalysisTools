
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
    // specify categories, for now do this manually
    catTypes.push_back("WenLowPt");    
    catTypes.push_back("WenHighPt");    
    catTypes.push_back("WmnLowPt");    
    catTypes.push_back("WmnMidPt");    
    catTypes.push_back("WmnHighPt");    

    // initialize histograms
    for (int j=0; j < (int)catTypes.size(); j++) {
        for(int i=0; i < (int)samples.size(); i++) {
            if(debug>1000) std::cout<<"Initializing histograms for sample: "<<samples[i].sampleName.c_str()<<std::endl;
            const char* sname = samples[i].sampleName.c_str();
            float nBinsX = *f["nBinsX"];
            float varMinX = *f["varMinX"];
            float varMaxX = *f["varMaxX"];
            float nBinsY = *f["nBinsY"];
            float varMinY = *f["varMinY"];
            float varMaxY = *f["varMaxY"];
            hists1D[catTypes[j]].push_back( new TH1F(sname, sname, nBinsX, varMinX, varMaxX) );
            hists2D[catTypes[j]].push_back(new TH2F(Form("%s_2D",sname), Form("%s_2D",sname), nBinsX, varMinX, varMaxX, nBinsY, varMinY, varMaxY) );
        }
    // create output file for histograms
    histout[catTypes[j]] = new TFile(Form("hists_%s.root",catTypes[j].c_str()), "recreate");
    }
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
    for(int j=0; j < (int)catTypes.size(); j++) {
        for(int i=0; i < (int)samples.size(); i++) {
            // fill histograms with current values of BDT, Hmass
            float BDT = *f["BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4"];
            float H_mass = *f["H_mass_f"];
            float H_mass_reg = *f["H_mass_reg"];
            int sNum = *in["sampleIndex"];        
            //float weight = samples[i].intWeight;
            float weight = *f["weight"];
            //std::cout<<"sNum = "<<sNum<<std::endl;
            //std::cout<<"samples[i].sampleNum = "<<samples[i].sampleNum<<std::endl; 
            if(sNum != samples[i].sampleNum) continue; 
            //std::cout<<sNum<<" == "<<samples[i].sampleNum<<std::endl;
            // match to category
            bool pass = false;
            //std::cout<<catTypes[j].c_str()<<std::endl;
            //std::cout<<"eventClass = "<<*in["eventClass"]<<std::endl;
            if (catTypes[j] == "WenLowPt" && *in["eventClass"]==1002) pass=true;
            if (catTypes[j] == "WenHighPt" && *in["eventClass"]==1001) pass=true;
            if (catTypes[j] == "WmnLowPt" && *in["eventClass"]==3) pass=true;
            if (catTypes[j] == "WmnMidPt" && *in["eventClass"]==2) pass=true;
            if (catTypes[j] == "WmnHighPt" && *in["eventClass"]==1) pass=true;
            if (!pass) continue;
            //std::cout<<"passed!"<<std::endl;
            ofile->cd();
            hists1D[catTypes[j]][i]->Fill(H_mass, weight);
            std::cout<<"1D hist filled with Mjj value: "<<H_mass<<" and weight: "<<weight<<std::endl;
            hists2D[catTypes[j]][i]->Fill(BDT, H_mass_reg, weight);
            if(debug>20000) std::cout<<"2D hist filled with BDT value: "<<BDT<<", Mjj value: "<<H_mass_reg<<" and weight: "<<weight<<std::endl;
        }
    }
}

void WorkspaceAnalysis::TermAnalysis() {

    // initialize all the RooHistPDF's
    float varMinX = *f["varMinX"];
    float varMaxX = *f["varMaxX"];
    float varMinY = *f["varMinY"];
    float varMaxY = *f["varMaxY"];
    for(int j=0; j < (int)catTypes.size(); j++) {
        RooWorkspace *WS = new RooWorkspace(catTypes[j].c_str(), catTypes[j].c_str());
        RooArgList *obs;
        WS->factory(Form("CMS_vhbb_BDT_Wln_8TeV[%f,%f]",varMinX,varMaxX));
        WS->factory(Form("CMS_vhbb_Mjj_Wln_8TeV[%f,%f]",varMinY,varMaxY));
  
        //TDirectory *wsdir = histout[catTypes[j]]->mkdir(catTypes[j].c_str());
        std::cout<<catTypes[j].c_str()<<std::endl; 
        for(int i=0; i < (int)samples.size(); i++) {
            // N.B.: This will only work correctly if the samples list hasn't changed at all since calling InitAnalysis()
            SampleContainer *csample = &samples[i];
            RooDataHist *tmpRDH = new RooDataHist();
            if(debug>1000) std::cout<<"modelDimension set to: "<<int(*f["modelDimension"])<<::std::endl;
            switch(int(*f["modelDimension"]))
            {
            case 1: 
                obs = new RooArgList(*WS->var("CMS_vhbb_Mjj_Wln_8TeV"));
                tmpRDH = new RooDataHist(csample->sampleName.c_str(), csample->sampleName.c_str(), *obs, hists1D[catTypes[j]][i]);
                break;
            case 2: 
                obs = new RooArgList(*WS->var("CMS_vhbb_BDT_Wln_8TeV"), *WS->var("CMS_vhbb_Mjj_Wln_8TeV"));
                tmpRDH = new RooDataHist(csample->sampleName.c_str(), csample->sampleName.c_str(), *obs, hists2D[catTypes[j]][i]);
                break;
            default:
                std::cerr<<"modelDimension must be set to 1 or 2! It is currently set to: "<<*f["modelDimension"]<<std::endl;
            }

            RooHistPdf *tmpRHP = new RooHistPdf(samples[i].sampleName.c_str(), samples[i].sampleName.c_str(), *obs, *tmpRDH);        
            if(debug>2000) std::cout<<"created RooHistPDF for sample: "<<samples[i].sampleName.c_str()<<std::endl;
            WS->import(*tmpRHP); 
            if(debug>2000) std::cout<<"imported RooHistPdf into workspace"<<std::endl; 
            histout[catTypes[j]]->cd();
            //wsdir->cd();
            hists1D[catTypes[j]][i]->Write();
            hists2D[catTypes[j]][i]->Write();
        }
        //histout[catTypes[j]]->Close();
    // do fits

    // delete variables/etc
    ofile->cd();
    WS->writeToFile(Form("%s_%s.root",outputTreeName.c_str(),catTypes[j].c_str()));
    }
    ofile->Close();
}
