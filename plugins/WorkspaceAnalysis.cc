
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
    // also split by reco top mass for study
    if (bool(*f["doTop1MassSplit"])) {
        catTypes.push_back("WenLowPt_0");    
        catTypes.push_back("WenLowPt_1");    
        catTypes.push_back("WenHighPt_0");    
        catTypes.push_back("WenHighPt_1");    
        catTypes.push_back("WmnLowPt_0");    
        catTypes.push_back("WmnLowPt_1");    
        catTypes.push_back("WmnMidPt_0");    
        catTypes.push_back("WmnMidPt_1");    
        catTypes.push_back("WmnHighPt_0"); 
        catTypes.push_back("WmnHighPt_1");
    } 
    else {
        catTypes.push_back("WenLowPt");    
        catTypes.push_back("WenHighPt");    
        catTypes.push_back("WmnLowPt");    
        catTypes.push_back("WmnMidPt");    
        catTypes.push_back("WmnHighPt"); 
    }
    reIndex[-12501] = -12501;
    reIndex[-12502] = -12502;
    reIndex[2200] = 2200;
    reIndex[2201] = 2201;
    reIndex[2202] = 2202;
    reIndex[2300] = 2300;
    reIndex[2301] = 2301;
    reIndex[2302] = 2302;
    reIndex[10] = 101;
    reIndex[11] = 101;
    reIndex[16] = 101;
    reIndex[17] = 101;
    reIndex[20] = 101;
    reIndex[21] = 101;
    reIndex[13] = 100;
    reIndex[14] = 100;
    reIndex[15] = 100;   

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
    if (*in["cutFlow"] < *f["cutFlowLevel"]) return false;
    if (*f["H_mass_f"] < 90 || *f["H_mass_f"]>150) return false; 
    // mass window, use parameterized mean/sigma from jet mass res. study
    double mean = 104.7 + 0.055 * (*f["H_pt_f"]);
    double sigma = 20.9* exp(-0.0024*(*f["H_pt_f"]));

    // parameterized mass window cut
    //if (*f["H_mass_f"] < (mean - 2*sigma) || *f["H_mass_f"] > (mean + 2*sigma)) return false;

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
            float H_pt = *f["H_pt_f"];
            int sNum = reIndex[*in["sampleIndex"]];        
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
            if (bool(*f["doTop1MassSplit"])) {
                if (catTypes[j] == "WenLowPt_0" && *in["eventClass"]==1002 && *f["Top1_mass_bestCSV"] <= *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WenLowPt_1" && *in["eventClass"]==1002 && *f["Top1_mass_bestCSV"] > *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WenHighPt_0" && *in["eventClass"]==1001 && *f["Top1_mass_bestCSV"] <= *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WenHighPt_1" && *in["eventClass"]==1001 && *f["Top1_mass_bestCSV"] > *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnLowPt_0" && *in["eventClass"]==3 && *f["Top1_mass_bestCSV"] <= *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnLowPt_1" && *in["eventClass"]==3 && *f["Top1_mass_bestCSV"] > *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnMidPt_0" && *in["eventClass"]==2 && *f["Top1_mass_bestCSV"] <= *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnMidPt_1" && *in["eventClass"]==2 && *f["Top1_mass_bestCSV"] > *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnHighPt_0" && *in["eventClass"]==1 && *f["Top1_mass_bestCSV"] <= *f["Top1MassSplit"]) pass=true;
                if (catTypes[j] == "WmnHighPt_1" && *in["eventClass"]==1 && *f["Top1_mass_bestCSV"] > *f["Top1MassSplit"]) pass=true;
            }
            else {

                if (catTypes[j] == "WenLowPt" && *in["eventClass"]==1002) pass=true;
                if (catTypes[j] == "WenHighPt" && *in["eventClass"]==1001) pass=true;
                if (catTypes[j] == "WmnLowPt" && *in["eventClass"]==3) pass=true;
                if (catTypes[j] == "WmnMidPt" && *in["eventClass"]==2) pass=true;
                if (catTypes[j] == "WmnHighPt" && *in["eventClass"]==1) pass=true;
            }
            // Let's go back to cutting directly on V_pt for categories so we can optimize
            if (catTypes[j] == "WenLowPt" && *in["isWenu"]&&*f["V_pt"]>100&&*f["V_pt"]<150) pass=true;
            if (catTypes[j] == "WenHighPt" && *in["isWenu"]&&*f["V_pt"]>150) pass=true;
            if (catTypes[j] == "WmnLowPt" && *in["isWmunu"]&&*f["V_pt"]>100&&*f["V_pt"]<130) pass=true;
            if (catTypes[j] == "WmnMidPt" && *in["isWmunu"]&&*f["V_pt"]>130&&*f["V_pt"]<180) pass=true;
            if (catTypes[j] == "WmnHighPt" && *in["isWmunu"]&&*f["V_pt"]>180) pass=true;
                        
            //if (*in["isWenu"] || *in["isWmunu"]) pass=true;
            if (!pass) continue;
            //std::cout<<"passed!"<<std::endl;
            ofile->cd();
            hists1D[catTypes[j]][i]->Fill(H_mass, weight);
            if(debug>20000) std::cout<<"1D hist filled with Mjj value: "<<H_mass<<" and weight: "<<weight<<std::endl;
            //hists2D[catTypes[j]][i]->Fill(H_mass, BDT, weight);
            //if(debug>20000) std::cout<<"2D hist filled with BDT value: "<<BDT<<", Mjj value: "<<H_mass<<" and weight: "<<weight<<std::endl;
            hists2D[catTypes[j]][i]->Fill(H_mass, H_pt, weight);
            if(debug>20000) std::cout<<"2D hist filled with H_pt value: "<<H_pt<<", Mjj value: "<<H_mass<<" and weight: "<<weight<<std::endl;
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
        //WS->factory(Form("CMS_vhbb_BDT_Wln_8TeV[%f,%f]",varMinY,varMaxY));
        WS->factory(Form("CMS_vhbb_Hpt_Wln_8TeV[%f,%f]",varMinY,varMaxY));
        WS->factory(Form("CMS_vhbb_Mjj_Wln_8TeV[%f,%f]",varMinX,varMaxX));
  
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
                //obs = new RooArgList(*WS->var("CMS_vhbb_Mjj_Wln_8TeV"), *WS->var("CMS_vhbb_BDT_Wln_8TeV"));
                obs = new RooArgList(*WS->var("CMS_vhbb_Mjj_Wln_8TeV"), *WS->var("CMS_vhbb_Hpt_Wln_8TeV"));
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
