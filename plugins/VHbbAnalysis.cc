//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "VHbbAnalysis.h"

// initialize parameters
VHbbAnalysis::VHbbAnalysis(){
    if(debug>10) std::cout<<"Constructing VHbbAnalysis"<<std::endl;
}

// remove any 'new'
VHbbAnalysis::~VHbbAnalysis(){
    if(debug>10) std::cout<<"Deconstructing VHbbAnalysis"<<std::endl;
}

//probably want to do bdtSetup here
void VHbbAnalysis::InitAnalysis(){
    SetupBDT();
    return;
}

//get a few entries and then preselect
//if sel, then analyzeevent
//default to false in the future
bool VHbbAnalysis::Preselection(){
    bool sel=true;
    return sel;
}

bool VHbbAnalysis::Analyze(){
    bool sel=false;
    //check if event passes any event class
    if(WmunuHbbSelection()) {
        sel=true;
        *in["eventClass"]=0;
    } else if(WenuHbbSelection()) {
        sel=true;
        *in["eventClass"]=1;
    }

    return sel;
}

void VHbbAnalysis::FinishEvent(){
    *in["sampleIndex"] = cursample->sampleNum;
    
    *f["Vtype_f"] = (float) *in["Vtype"];
    *f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["naLeptonsPassingCuts"] = 0.;
    for(int j=0;j<*in["nalep"]&&j<100;j++){
        if(f["aLepton_pt"][j]>15. && fabs(f["aLepton_eta"][j])<2.5 && f["aLepton_pfCombRelIso"][j]<0.15) *f["naLeptonsPassingCuts"] += 1.;
    }
    *f["naJetsPassingCuts"] = 0.;
    *f["highestCSVaJet"] = 0.0;
    *f["minDeltaRaJet"] = 9999.0;
    for(int j=0;j<*in["naJets"];j++){
        if(f["aJet_pt"][j]>20. && fabs(f["aJet_eta"][j])<4.5) *f["naJetsPassingCuts"] += 1;
        if(f["aJet_pt"][j]>20. && fabs(f["aJet_eta"][j])<2.5 && f["aJet_csvReshapedNew"][j]>*f["highestCSVaJet"]) *f["highestCSVaJet"] = f["aJet_csvReshapedNew"][j];
        if(f["aJet_pt"][j]>20. && fabs(f["aJet_eta"][j])<4.5 && f["aJet_puJetIdL"][j]>0 && EvalDeltaR(f["aJet_eta"][j],f["aJet_phi"][j],H.eta,H.phi)<*f["minDeltaRaJet"]) {
            *f["minDeltaRaJet"] = EvalDeltaR(f["aJet_eta"][j],f["aJet_phi"][j],H.eta,H.phi);
    //*f["BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkgj_newCuts4"] = thereader->EvaluateMVA("BDT method");
        }
    }       
    
    ofile->cd();
    outputTree->Fill();
    return;
}

void VHbbAnalysis::TermAnalysis(){
    ofile->cd();
    outputTree->Write();
    ofile->Close();
    return;
}


bool VHbbAnalysis::WenuHbbSelection(){
    int selectedInd=-1;
    bool selectEvent=false;
    float highestPt=-99;
    if(debug>1000){
        std::cout<<"*in[\"nvlep\"] "<<*in["nvlep"]<<std::endl;
    }
    for(int elInd=0; elInd<*in["nvlep"]; elInd++) {
        if(debug>1000){
            std::cout<<"f[\"vLepton_pt\"][elInd] "<<f["vLepton_pt"][elInd]<<std::endl;
            std::cout<<"in[\"vLepton_wp85\"][elInd] "<<in["vLepton_wp85"][elInd]<<std::endl;
        }
        if(f["vLepton_pt"][elInd] > 20 && in["vLepton_wp85"][elInd]>0 && highestPt<f["vLepton_pt"][elInd]){
            highestPt=f["vLepton_pt"][elInd];
            selectedInd=elInd;
            selectEvent=true;
        }
    }
    
    return selectEvent;
}


bool VHbbAnalysis::WmunuHbbSelection(){
    bool sel=false;
    //if(1) sel=true;
    return sel;
}

