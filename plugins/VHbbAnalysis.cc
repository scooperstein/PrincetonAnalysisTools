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
    //if( *d["Vtype"]==3 ) sel=true;
    if( *d["Vtype"]==3) sel=true;
    return sel;
}

bool VHbbAnalysis::Analyze(){
    bool sel=false;

    if(debug>1000) {
        std::cout<<"selecting bjets"<<std::endl;
    }
    std::pair<int,int> bjets=HighestPtBJets();
   
    // there aren't two acceptable jets
    //if(bjets.first==-1 || bjets.second==-1) return sel;
    *in["hJetInd1"]=bjets.first;
    *in["hJetInd2"]=bjets.second;

    if(debug>1000) {
        std::cout<<"found two bjets with pt "
            <<d["Jet_pt"][bjets.first]<<" "
            <<d["Jet_pt"][bjets.second]<<" "
            <<std::endl;
    }
    
    // For now use higgs bjet selection from step 2 ntuples
    if(d["hJets_btagCSV"][0] < *f["j1ptCSV"] || d["hJets_btagCSV"][1] < *f["j2ptCSV"]) return false;
 
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
    
    // General use variables
    *in["sampleIndex"] = cursample->sampleNum;
   
    // Split WJets and ZJets samples by jet parton flavor
    if (cursample->doJetFlavorSplit) {
        int nBJets = 0;
        //std::cout<<fabs(d["hJets_mcFlavour"][0])<<std::endl;
        if (fabs(in["hJets_mcFlavour"][0]) == 5) nBJets++;
        if (fabs(in["hJets_mcFlavour"][1]) == 5) nBJets++;
        *in["sampleIndex"] = (*in["sampleIndex"]*1000 + nBJets);
    } 

    // Split by boost category
    if(*d["V_pt"] >= 100 && *d["V_pt"] < 150) *in["eventClass"] = 1;
    else if(*d["V_pt"] >= 150) *in["eventClass"] = 0;
    else *in["eventClass"] = -1;

    *f["weight"] = cursample->intWeight;
    *f["Vtype_f"] = (float) *d["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["selLeptons_pt_0"] = d["selLeptons_pt"][0];
    *f["selLeptons_eta_0"] = d["selLeptons_eta"][0];
    *f["selLeptons_phi_0"] = d["selLeptons_phi"][0];
    *in["selLeptons_pdgId_0"] = in["selLeptons_pdgId"][0];    
    *in["selLeptons_eleCutIdCSA14_25ns_v1_0"] = in["selLeptons_eleCutIdCSA14_25ns_v1"][0];
    *f["selLeptons_relIso03_0"] = d["selLeptons_relIso03"][0];  
    

    *f["naLeptonsPassingCuts"] = 0.;
    for(int j=0;j<*in["naLeptons"]&&j<100;j++){
        //if(f["aLepton_pt"][j]>15. && fabs(f["aLepton_eta"][j])<2.5 && f["aLepton_pfCombRelIso"][j]<0.15) *f["naLeptonsPassingCuts"] += 1.;
        if(d["aLeptons_pt"][j]>15. && fabs(d["aLeptons_eta"][j])<2.5 && in["aLeptons_eleCutIdCSA14_25ns_v1"][j]) *f["naLeptonsPassingCuts"] += 1.;
    }
    *f["naJetsPassingCuts"] = 0.;
    *f["highestCSVaJet"] = 0.0;
    *f["minDeltaRaJet"] = 9999.0; // FIXME add correct initialization
    for(int j=0;j<*in["naJets"];j++){
        if(d["aJets_pt"][j]>20. && fabs(d["aJets_eta"][j])<4.5) *f["naJetsPassingCuts"] += 1;
        if(d["aJets_pt"][j]>20 && fabs(d["aJets_eta"][j])<2.5 && d["aJets_btagCSV"][j]>*f["highestCSVaJet"]) *f["highestCSVaJet"] = d["aJets_btagCSV"][j];
      }       
    
    // Set variables used by the BDT
    *f["H_mass_f"] = (float) *d["H_mass"];
    *f["H_pt_f"] = (float) *d["H_pt"];
    *f["V_pt_f"] = (float) *d["V_pt"];
    *f["hJets_btagCSV_0"] = (float) d["hJets_btagCSV"][0];
    *f["hJets_btagCSV_1"] = (float) d["hJets_btagCSV"][1]; 
    *f["HVdPhi_f"] = (float) *d["HVdPhi"];
    *f["H_dEta"] = fabs(d["hJets_eta"][0] - d["hJets_eta"][1]);
    
    TLorentzVector hJ1 = TLorentzVector();
    TLorentzVector hJ2 = TLorentzVector();
    hJ1.SetPtEtaPhiM(d["hJets_pt"][0], d["hJets_eta"][0], d["hJets_phi"][0], d["hJets_mass"][0]);
    hJ2.SetPtEtaPhiM(d["hJets_pt"][1], d["hJets_eta"][0], d["hJets_phi"][1], d["hJets_mass"][1]);
    *f["hJets_mt_0"] = hJ1.Mt();
    *f["hJets_mt_1"] = hJ2.Mt();
    *f["H_dR"] = (float) hJ1.DeltaR(hJ2);   
    *f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    *f["hJet_ptCorr_0"] = (float) d["hJets_pt"][0];
    *f["hJet_ptCorr_1"] = (float) d["hJets_pt"][1];
    *f["met_sumEt_f"] = (float) *d["met_sumEt"]; // is this the right variable??
 
    for (unsigned int i=0; i<bdtInfos.size(); i++) {

       BDTInfo tmpBDT = bdtInfos[i];
       if(debug>5000) {
           PrintBDTInfoValues(tmpBDT);
           std::cout<<"BDT evaulates to: "<<tmpBDT.reader->EvaluateMVA(tmpBDT.bdtname)<<std::endl;
       }
       *f[tmpBDT.bdtname] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtname);
    }

    if(jet1EnergyRegressionIsSet && jet2EnergyRegressionIsSet) {
        *f["hJets_pt_0"] = float(d["hJets_pt"][0]);
        //*f["hJets_rawPt_0"] = float(d["hJets_rawPt"][0]);
        *f["hJets_eta_0"] = float(d["hJets_eta"][0]);
        //*f["hJets_mt_0"] = 0.; // FIXME
        *f["hJets_leadTrackPt_0"] = float(d["hJets_leadTrackPt"][0]);
        *f["hJets_leptonPtRel_0"] = float(d["hJets_leptonPtRel"][0]);
        *f["hJets_leptonPt_0"] = float(d["hJets_leptonPt"][0]);
        *f["hJets_leptonDeltaR_0"] = float(d["hJets_leptonDeltaR"][0]);
        *f["hJets_chEmEF_0"] = float(d["hJets_chEmEF"][0]);
        *f["hJets_chHEF_0"] = float(d["hJets_chHEF"][0]);
        *f["hJets_neHEF_0"] = float(d["hJets_neHEF"][0]);
        *f["hJets_neEmEF_0"] = float(d["hJets_neEmEF"][0]);
        *f["hJets_chMult_0"] = float(d["hJets_chMult"][0]);
        //*f["hJets_puId_0"] = float(d["hJets_puId"][0]);
        *f["hJets_vtx3DVal_0"] = 0.;
        *f["hJets_vtxNtracks_0"] = 0.;
        *f["hJets_vtxMass_0"] = float(d["hJets_vtxMass"][0]);
        *f["hJets_vtxPt_0"] = float(d["hJets_vtxPt"][0]);
        *f["hJets_vtx3DVal_0"] = float(d["hJets_vtx3DVal"][0]);
        *f["hJets_vtxNtracks_0"] = float(d["hJets_vtxNtracks"][0]);

        *f["hJets_pt_1"] = float(d["hJets_pt"][1]);
        //*f["hJets_rawPt_1"] = float(d["hJets_rawPt"][1]);
        *f["hJets_eta_1"] = float(d["hJets_eta"][1]);
        //*f["hJets_mt_1"] = 0.; // FIXME
        *f["hJets_leadTrackPt_1"] = float(d["hJets_leadTrackPt"][1]);
        *f["hJets_leptonPtRel_1"] = float(d["hJets_leptonPtRel"][1]);
        *f["hJets_leptonPt_1"] = float(d["hJets_leptonPt"][1]);
        *f["hJets_leptonDeltaR_1"] = float(d["hJets_leptonDeltaR"][1]);
        *f["hJets_chEmEF_1"] = float(d["hJets_chEmEF"][1]);
        *f["hJets_chHEF_1"] = float(d["hJets_chHEF"][1]);
        *f["hJets_neHEF_1"] = float(d["hJets_neHEF"][1]);
        *f["hJets_neEmEF_1"] = float(d["hJets_neEmEF"][1]);
        *f["hJets_chMult_1"] = float(d["hJets_chMult"][1]);
        //*f["hJets_puId_1"] = float(d["hJets_puId"][1]);
        *f["hJets_vtx3DVal_1"] = 0.;
        *f["hJets_vtxNtracks_1"] = 0.;
        *f["hJets_vtxMass_1"] = float(d["hJets_vtxMass"][1]);
        *f["hJets_vtx3DVal_1"] = float(d["hJets_vtx3DVal"][1]);
        *f["hJets_vtxPt_1"] = float(d["hJets_vtxPt"][1]);
        *f["hJets_vtxNtracks_1"] = float(d["hJets_vtxNtracks"][1]);
        
        if(debug>10000) {
            std::cout<<"Evaluating the Jet Energy Regression..."<<std::endl;
            PrintBDTInfoValues(jet1EnergyRegression);
            PrintBDTInfoValues(jet2EnergyRegression);
        }
        double r1Pt = jet1EnergyRegression.reader->EvaluateRegression(jet1EnergyRegression.bdtname)[0];
        double r2Pt = jet2EnergyRegression.reader->EvaluateRegression(jet2EnergyRegression.bdtname)[0];

        *f["Jet1_regWeight"] = r1Pt/(*f["hJets_pt_0"]);
        *f["Jet2_regWeight"] = r2Pt/(*f["hJets_pt_1"]);

        TLorentzVector hJ1_reg = TLorentzVector();
        TLorentzVector hJ2_reg = TLorentzVector();
 
        *f["Jet1_pt_reg"] = r1Pt;
        *f["Jet2_pt_reg"] = r2Pt;
        hJ1_reg.SetPtEtaPhiM(r1Pt, *f["hJets_eta_0"], d["hJets_phi"][0], d["hJets_mass"][0]);
        hJ2_reg.SetPtEtaPhiM(r2Pt, *f["hJets_eta_1"], d["hJets_phi"][1], d["hJets_mass"][1]);
       
        //*f["hJets_mt_0"] = hJ1_reg.Mt();
        //*f["hJets_mt_1"] = hJ2_reg.Mt();
    
        TLorentzVector H_vec_reg = hJ1_reg + hJ2_reg;
        *f["H_mass_reg"] = H_vec_reg.M();
        *f["H_pt_reg"] = H_vec_reg.Pt();
        *f["H_eta_reg"] = H_vec_reg.Eta();
        *f["H_phi_reg"] = H_vec_reg.Phi();
        *f["H_dR_reg"] = hJ1_reg.DeltaR(hJ2_reg);
        *f["H_dPhi_reg"] = hJ1_reg.DeltaPhi(hJ2_reg);
        *f["H_dEta_reg"] = fabs(hJ1_reg.Eta() - hJ2_reg.Eta() );
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
        
        TLorentzVector W,El, MET, Hbb, HJ1, HJ2;
        // Reconstruct W
        MET.SetPtEtaPhiM(*d["met_pt"], 0., *d["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
        El.SetPtEtaPhiM(d["selLeptons_pt"][0], d["selLeptons_eta"][0], d["selLeptons_phi"][0], d["selLeptons_mass"][0]); 
        W = MET + El; 
        //*d["V_pt"] = W.Pt() // uncomment this line if we want to recalculate W.pt ourselves
 
        // Reconstruct Higgs
        HJ1.SetPtEtaPhiM(d["hJets_pt"][0], d["hJets_eta"][0], d["hJets_phi"][0], d["hJets_mass"][0]);
        HJ1.SetPtEtaPhiM(d["hJets_pt"][1], d["hJets_eta"][1], d["hJets_phi"][1], d["hJets_mass"][1]);
        Hbb = HJ1 + HJ2;
        
        // Now we can calculate whatever we want (transverse) with W and H four-vectors

        if(*d["elMetDPhi"] < *f["elMetDPhiCut"] && *d["HVdPhi"]> *f["HVDPhiCut"] && *d["V_pt"] > *f["vptcut"]
            && *d["H_pt"] > *f["hptcut"]  ){
            
            if (*in["naJets"] > 0 || *in["naLeptons"] > 0) return false;
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

