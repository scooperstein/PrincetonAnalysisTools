//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "VHbbAnalysis.h"

#include "TLorentzVector.h"

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
    //SetupBDT();a//FIXME need to update BDT
    return;
}

//get a few entries and then preselect
//if sel, then analyzeevent
//default to false in the future
bool VHbbAnalysis::Preselection(){
    bool sel=false;
    if( *d["Vtype"]==3 ) sel=true;
    return sel;
}

bool VHbbAnalysis::Analyze(){
    bool sel=false;

    if(debug>1000) {
        std::cout<<"selecting bjets"<<std::endl;
    }
    std::pair<int,int> bjets=HighestPtBJets();
   
    // there aren't two acceptable jets
    if(bjets.first==-1 || bjets.second==-1) return sel;
    *in["hJetInd1"]=bjets.first;
    *in["hJetInd2"]=bjets.second;

    if(debug>1000) {
        std::cout<<"found two bjets with pt "
            <<d["Jet_pt"][bjets.first]<<" "
            <<d["Jet_pt"][bjets.second]<<" "
            <<std::endl;
    }
 
    //check if event passes any event class
    if(WmunuHbbSelection()) {
        sel=true;
        *in["eventClass"]=0;
    } else if(WenuHbbSelection()) {
        sel=true;
        *in["eventClass"]=1;
    }

    if(debug>1000) {
        std::cout<<"selecting event"<<std::endl;
    }

    return sel;
}

void VHbbAnalysis::FinishEvent(){
    *in["sampleIndex"] = cursample->sampleNum;
    *f["weight"] = cursample->intWeight;
    
    *f["Vtype_f"] = (float) *d["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["naLeptonsPassingCuts"] = 0.;
    //for(int j=0;j<*in["naLeptons"]&&j<100;j++){
    //    //if(f["aLepton_pt"][j]>15. && fabs(f["aLepton_eta"][j])<2.5 && f["aLepton_pfCombRelIso"][j]<0.15) *f["naLeptonsPassingCuts"] += 1.;
    //    if(f["aLepton_pt"][j]>15. && fabs(f["aLepton_eta"][j])<2.5 && f["aLeptons_eleCutIdCSA14_25ns_v1"][j]) *f["naLeptonsPassingCuts"] += 1.;
    //}
    //*f["naJetsPassingCuts"] = 0.;
    //*f["highestCSVaJet"] = 0.0;
    //*f["minDeltaRaJet"] = 9999.0;
    //for(int j=0;j<*in["naJets"];j++){
    //    if(f["aJet_pt"][j]>20. && fabs(f["aJet_eta"][j])<4.5) *f["naJetsPassingCuts"] += 1;
    //}       
    
    for (unsigned int i=0; i<bdtInfos.size(); i++) {
       BDTInfo tmpBDT = bdtInfos[i];
       if(debug>5000) PrintBDTInfoValues(tmpBDT);
       std::cout<<tmpBDT.reader->EvaluateMVA(tmpBDT.bdtname)<<std::endl;
       *f[tmpBDT.bdtname] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtname);
    }

    if(jetEnergyRegressionIsSet) {
        *f["hJets_pt_0"] = float(d["hJets_pt"][0]);
        *f["hJets_rawPt_0"] = float(d["hJets_rawPt"][0]);
        *f["hJets_eta_0"] = float(d["hJets_eta"][0]);
        *f["hJets_mt_0"] = 0.; // FIXME
        *f["hJets_leadTrackPt_0"] = float(d["hJets_leadTrackPt"][0]);
        *f["hJets_leptonPtRel_0"] = float(d["hJets_leptonPtRel"][0]);
        *f["hJets_leptonPt_0"] = float(d["hJets_leptonPt"][0]);
        *f["hJets_leptonDeltaR_0"] = float(d["hJets_leptonDeltaR"][0]);
        *f["hJets_chEmEF_0"] = float(d["hJets_chEmEF"][0]);
        *f["hJets_chHEF_0"] = float(d["hJets_chHEF"][0]);
        *f["hJets_neHEF_0"] = float(d["hJets_neHEF"][0]);
        *f["hJets_neEmEF_0"] = float(d["hJets_neEmEF"][0]);
        *f["hJets_chMult_0"] = float(d["hJets_chMult"][0]);
        *f["hJets_puId_0"] = float(d["hJets_puId"][0]);
        *f["hJets_vtxMass_0"] = 0.; // FIXME
        *f["hJets_vtx3DVal_0"] = 0.;
        *f["hJets_vtxNtracks_0"] = 0.;
        //*f["hJets_vtxMass_0"] = float(d["hJets_vtxMass"][0]);
        //*f["hJets_vtx3DVal_0"] = float(d["hJets_vtx3DVal"][0]);
        //*f["hJets_vtxNtracks_0"] = float(d["hJets_vtxNtracks"][0]);

        *f["hJets_pt_1"] = float(d["hJets_pt"][1]);
        *f["hJets_rawPt_1"] = float(d["hJets_rawPt"][1]);
        *f["hJets_eta_1"] = float(d["hJets_eta"][1]);
        *f["hJets_mt_1"] = 0.; // FIXME
        *f["hJets_leadTrackPt_1"] = float(d["hJets_leadTrackPt"][1]);
        *f["hJets_leptonPtRel_1"] = float(d["hJets_leptonPtRel"][1]);
        *f["hJets_leptonPt_1"] = float(d["hJets_leptonPt"][1]);
        *f["hJets_leptonDeltaR_1"] = float(d["hJets_leptonDeltaR"][1]);
        *f["hJets_chEmEF_1"] = float(d["hJets_chEmEF"][1]);
        *f["hJets_chHEF_1"] = float(d["hJets_chHEF"][1]);
        *f["hJets_neHEF_1"] = float(d["hJets_neHEF"][1]);
        *f["hJets_neEmEF_1"] = float(d["hJets_neEmEF"][1]);
        *f["hJets_chMult_1"] = float(d["hJets_chMult"][1]);
        *f["hJets_puId_1"] = float(d["hJets_puId"][1]);
        *f["hJets_vtxMass_1"] = 0.; // FIXME
        *f["hJets_vtx3DVal_1"] = 0.;
        *f["hJets_vtxNtracks_1"] = 0.;
        //*f["hJets_vtxMass_1"] = float(d["hJets_vtxMass"][1]);
        //*f["hJets_vtx3DVal_1"] = float(d["hJets_vtx3DVal"][1]);
        //*f["hJets_vtxNtracks_1"] = float(d["hJets_vtxNtracks"][1]);
        
        if(debug>10000) {
            std::cout<<"Evaluating the Jet Energy Regression..."<<std::endl;
            PrintBDTInfoValues(jet1EnergyRegression);
            PrintBDTInfoValues(jet2EnergyRegression);
        }
        double r1Pt = jet1EnergyRegression.reader->EvaluateRegression(jet1EnergyRegression.bdtname)[0];
        double r2Pt = jet2EnergyRegression.reader->EvaluateRegression(jet2EnergyRegression.bdtname)[0];

        *f["Jet1_regWeight"] = r1Pt/(*f["hJets_pt_0"]);
        *f["Jet2_regWeight"] = r2Pt/(*f["hJets_pt_1"]);

        TLorentzVector hJ1 = TLorentzVector();
        TLorentzVector hJ2 = TLorentzVector();
 
        *f["Jet1_pt_reg"] = r1Pt;
        *f["Jet2_pt_reg"] = r2Pt;
        hJ1.SetPtEtaPhiM(r1Pt, *f["hJets_eta_0"], d["hJets_phi"][0], d["hJets_mass"][0]);
        hJ2.SetPtEtaPhiM(r2Pt, *f["hJets_eta_1"], d["hJets_phi"][1], d["hJets_mass"][1]);
           
        TLorentzVector H_vec = hJ1 + hJ2;
        *f["H_mass_reg"] = H_vec.M();
        *f["H_pt_reg"] = H_vec.Pt();
        *f["H_eta_reg"] = H_vec.Eta();
        *f["H_phi_reg"] = H_vec.Phi();
        *f["H_dR_reg"] = hJ1.DeltaR(hJ2);
        *f["H_dPhi_reg"] = hJ1.DeltaPhi(hJ2);
        *f["H_dEta_reg"] = fabs(hJ1.Eta() - hJ2.Eta() );
     }

    ofile->cd();
    outputTree->Fill();
    return;
}

void VHbbAnalysis::TermAnalysis(){
    if(debug>10) std::cout<<"START TermAnalysis()"<<std::endl;
    ofile->cd();
    outputTree->Write();
    ofile->Close();
    if(debug>10) std::cout<<"DONE TermAnalysis()"<<std::endl;
    return;
}


bool VHbbAnalysis::WenuHbbSelection(){
    int selectedInd=-1;
    bool selectEvent=false;
    float highestPt=-99;
    if(debug>1000){
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<d["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<d["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*d[\"met_pt\"] "<<*d["met_pt"]<<std::endl;
    }

    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts
    if(fabs(in["selLeptons_pdgId"][0])==11 
        && d["selLeptons_pt"][0]      > *f["eptcut"] 
        && d["selLeptons_relIso03"][0]< *f["erelisocut"]
        && *d["met_pt"]               > *f["metcut"]){
        *d["elMetDPhi"]=fabs(EvalDeltaPhi(d["selLeptons_phi"][0],*d["met_phi"]));
        //*d["HVDPhi"]   =fabs(EvalDeltaPhi(d["selLeptons_phi"][0],*d["met_phi"]));
        
        if(*d["elMetDPhi"] < *f["elMetDPhiCut"]){
            //&& *d["HVDPhi"]> *f["HVDPhiCut"]  ){
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


std::pair<int,int> VHbbAnalysis::HighestPtBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j1ptCut"]
            && d["Jet_btagCSV"][i]>*f["j1ptCSV"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(d["Jet_pt"][pair.first]<d["Jet_pt"][i]){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j2ptCut"]
            && d["Jet_btagCSV"][i]>*f["j2ptCSV"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(d["Jet_pt"][pair.second]<d["Jet_pt"][i]){
                pair.second = i;
            }
        }
    }

    return pair;
}


std::pair<int,int> VHbbAnalysis::HighestCSVBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j1ptCut"]
            && d["Jet_btagCSV"][i]>*f["j1ptCSV"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(d["Jet_btagCSV"][pair.first]<d["Jet_btagCSV"][i]){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j2ptCut"]
            && d["Jet_btagCSV"][i]>*f["j2ptCSV"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(d["Jet_btagCSV"][pair.second]<d["Jet_btagCSV"][i]){
                pair.second = i;
            }
        }
    }

    return pair;
}

