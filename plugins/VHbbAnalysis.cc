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
    if( *d["Vtype"]>=0 && *d["Vtype"]<=4) sel=true;
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
    *in["nHJetsMatched"] = 0;
    if(d["hJets_pt"][0] == d["Jet_pt"][*in["hJetInd1"]]) *in["nHJetsMatched"] += 1;
    if(d["hJets_pt"][1] == d["Jet_pt"][*in["hJetInd2"]]) *in["nHJetsMatched"] += 1; 
    
    if(debug>10000) {
        if(*in["nHJetsMatched"]<2) {
            std::cout<<"nHJetsMatched = "<<*in["nHJetsMatched"]<<std::endl;
            std::cout<<"hJets_pt[0] = "<<d["hJets_pt"][0]<<std::endl;
            std::cout<<"hJets_pt[1] = "<<d["hJets_pt"][1]<<std::endl;
            std::cout<<"lead jet1 pt = "<<d["Jet_pt"][*in["hJetInd1"]]<<std::endl;
            std::cout<<"sublead jet2 pt = "<<d["Jet_pt"][*in["hJetInd2"]]<<std::endl;
        }
    }

    /*// For now use higgs bjet selection from step 2 ntuples
    if(d["hJets_btagCSV"][0] < *f["j1ptCSV"] || d["hJets_btagCSV"][1] < *f["j2ptCSV"]) return false;
    if(d["hJets_pt"][0] < *f["j1ptCut"] || d["hJets_pt"][1] < *f["j2ptCut"]) return false;
    */
 
    // Cut on the bjets that we select
    if(d["Jet_btagCSV"][*in["hJetInd1"]] < *f["j1ptCSV"] || d["Jet_btagCSV"][*in["hJetInd2"]] < *f["j2ptCSV"]) return false;
    if(d["Jet_pt"][*in["hJetInd1"]] < *f["j1ptCut"] || d["Jet_pt"][*in["hJetInd2"]] < *f["j2ptCut"]) return false; 
     

    //check if event passes any event class
    *in["lepInd"] = -1;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    if(WmunuHbbSelection()) {
        sel=true;
        *in["isWmunu"] = 1;
        *in["eventClass"]=0;
        *in["lepInd"] = *in["muInd"];
    } else if(WenuHbbSelection()) {
        sel=true;
        *in["isWenu"] = 1;
        if (!(*in["isWmunu"])) {
            *in["eventClass"]=1000;
            *in["lepInd"] = *in["elInd"];
        }
    }
    
    // count the number of additional leptons and jets, then cut on this number
    int nAddJet = 0;
    int nAddLep = 0;          
    for(int i=0; i < *in["naJets"]; i++) {
        if(d["aJets_pt"][i]>20 && fabs(d["aJets_eta"][i])<4.5 && in["aJets_id"][i]>0) {
            nAddJet++;
        }  
    }           
    for(int i=0; i < *in["naLeptons"]; i++) {
        if(d["aLeptons_pt"][i]>15 && fabs(d["aLeptons_eta"][i])<2.5 && d["aLeptons_relIso03"][i]<0.1) {
            nAddLep++;
        }
    }
    *in["nAddJets"] = nAddJet;
    *in["nAddLeptons"] = nAddLep;
    if(nAddJet >= *f["nAddJetsCut"] || nAddLep>= *f["nAddLeptonsCut"]) return false; 
    
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
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd1"]]) == 5) nBJets++;
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd2"]]) == 5) nBJets++;
        *in["sampleIndex"] = (*in["sampleIndex"]*1000 + nBJets);
    } 

    // Split by boost category
    if(*in["isWenu"]) {
        if(*d["V_pt"] >= 100 && *d["V_pt"] < 150) *in["eventClass"] += 2;
        else if(*d["V_pt"] >= 150) *in["eventClass"] += 1;
        else *in["eventClass"] += 3;
    }
    else if(*in["isWmunu"]) {
        if(*d["V_pt"] >= 100 && *d["V_pt"] < 130) *in["eventClass"] += 3;
        else if(*d["V_pt"] >= 130 && *d["V_pt"] < 180) *in["eventClass"] += 2;
        else if(*d["V_pt"] >= 180) *in["eventClass"] += 1;
        else *in["eventClass"] += 4;
    }

    // debugging
    *d["hJets_pt_0_Andrea"] = d["hJets_pt"][0];
    *d["hJets_pt_1_Andrea"] = d["hJets_pt"][1];
    *d["hJets_btagCSV_0_Andrea"] = d["hJets_btagCSV"][0];
    *d["hJets_btagCSV_1_Andrea"] = d["hJets_btagCSV"][1];

    *f["weight"] = cursample->intWeight;
    *f["Vtype_f"] = (float) *d["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["selLeptons_pt_0"] = d["selLeptons_pt"][*in["lepInd"]];
    *f["selLeptons_eta_0"] = d["selLeptons_eta"][*in["lepInd"]];
    *f["selLeptons_phi_0"] = d["selLeptons_phi"][*in["lepInd"]];
    *in["selLeptons_pdgId_0"] = in["selLeptons_pdgId"][*in["lepInd"]];    
    *in["selLeptons_eleCutIdCSA14_25ns_v1_0"] = in["selLeptons_eleCutIdCSA14_25ns_v1"][*in["lepInd"]];
    *in["selLeptons_tightId_0"] = in["selLeptons_tightId"][*in["lepInd"]];
    *f["selLeptons_relIso03_0"] = d["selLeptons_relIso03"][*in["lepInd"]];  
   
    if (*in["nGenLep"] == 0) {
        // gen lep is originally a tau?
        *f["selLeptons_genEleDR_0"] = -1;
    }
    else {
        TLorentzVector GenLep, El;
        GenLep.SetPtEtaPhiM(d["GenLep_pt"][0],d["GenLep_eta"][0],d["GenLep_phi"][0],d["GenLep_mass"][0]);
        El.SetPtEtaPhiM(d["selLeptons_pt"][*in["lepInd"]], d["selLeptons_eta"][*in["lepInd"]], d["selLeptons_phi"][*in["lepInd"]], d["selLeptons_mass"][*in["lepInd"]]);
        *f["selLeptons_genEleDR_0"] = El.DeltaR(GenLep);
    }
 

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
    *f["hJets_btagCSV_0"] = (float) d["Jet_btagCSV"][*in["hJetInd1"]];
    *f["hJets_btagCSV_1"] = (float) d["Jet_btagCSV"][*in["hJetInd2"]]; 
    *f["HVdPhi_f"] = (float) *d["HVdPhi"];
    *f["H_dEta"] = fabs(d["Jet_eta"][*in["hJetInd1"]] - d["Jet_eta"][*in["hJetInd2"]]);
   
    TLorentzVector hJ1 = TLorentzVector();
    TLorentzVector hJ2 = TLorentzVector();
    hJ1.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd1"]], d["Jet_eta"][*in["hJetInd1"]], d["Jet_phi"][*in["hJetInd1"]], d["Jet_mass"][*in["hJetInd1"]]);
    hJ2.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd2"]], d["Jet_eta"][*in["hJetInd2"]], d["Jet_phi"][*in["hJetInd2"]], d["Jet_mass"][*in["hJetInd2"]]);
    *f["hJets_mt_0"] = hJ1.Mt();
    *f["hJets_mt_1"] = hJ2.Mt();
    *f["H_dR"] = (float) hJ1.DeltaR(hJ2);   
    *f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    *f["hJet_ptCorr_0"] = (float) d["Jet_pt"][*in["hJetInd1"]];
    *f["hJet_ptCorr_1"] = (float) d["Jet_pt"][*in["hJetInd2"]];
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
        *f["hJets_pt_0"] = float(d["Jet_pt"][*in["hJetInd1"]]);
        //*f["hJets_rawPt_0"] = float(d["hJets_rawPt"][0]);
        *f["hJets_eta_0"] = float(d["Jet_eta"][*in["hJetInd1"]]);
        //*f["hJets_mt_0"] = 0.; // FIXME
        *f["hJets_leadTrackPt_0"] = float(d["Jet_leadTrackPt"][*in["hJetInd1"]]);
        *f["hJets_leptonPtRel_0"] = float(d["Jet_leptonPtRel"][*in["hJetInd1"]]);
        *f["hJets_leptonPt_0"] = float(d["Jet_leptonPt"][*in["hJetInd1"]]);
        *f["hJets_leptonDeltaR_0"] = float(d["Jet_leptonDeltaR"][*in["hJetInd1"]]);
        *f["hJets_chEmEF_0"] = float(d["Jet_chEmEF"][*in["hJetInd1"]]);
        *f["hJets_chHEF_0"] = float(d["Jet_chHEF"][*in["hJetInd1"]]);
        *f["hJets_neHEF_0"] = float(d["Jet_neHEF"][*in["hJetInd1"]]);
        *f["hJets_neEmEF_0"] = float(d["Jet_neEmEF"][*in["hJetInd1"]]);
        *f["hJets_chMult_0"] = float(d["Jet_chMult"][*in["hJetInd1"]]);
        //*f["hJets_puId_0"] = float(d["hJets_puId"][0]);
        *f["hJets_vtx3DVal_0"] = 0.;
        *f["hJets_vtxNtracks_0"] = 0.;
        *f["hJets_vtxMass_0"] = float(d["Jet_vtxMass"][*in["hJetInd1"]]);
        *f["hJets_vtxPt_0"] = float(d["Jet_vtxPt"][*in["hJetInd1"]]);
        *f["hJets_vtx3DVal_0"] = float(d["Jet_vtx3DVal"][*in["hJetInd1"]]);
        *f["hJets_vtxNtracks_0"] = float(d["Jet_vtxNtracks"][*in["hJetInd1"]]);

        *f["hJets_pt_1"] = float(d["Jet_pt"][*in["hJetInd2"]]);
        //*f["hJets_rawPt_1"] = float(d["hJets_rawPt"][1]);
        *f["hJets_eta_1"] = float(d["Jet_eta"][*in["hJetInd2"]]);
        //*f["hJets_mt_1"] = 0.; // FIXME
        *f["hJets_leadTrackPt_1"] = float(d["Jet_leadTrackPt"][*in["hJetInd2"]]);
        *f["hJets_leptonPtRel_1"] = float(d["Jet_leptonPtRel"][*in["hJetInd2"]]);
        *f["hJets_leptonPt_1"] = float(d["Jet_leptonPt"][*in["hJetInd2"]]);
        *f["hJets_leptonDeltaR_1"] = float(d["Jet_leptonDeltaR"][*in["hJetInd2"]]);
        *f["hJets_chEmEF_1"] = float(d["Jet_chEmEF"][*in["hJetInd2"]]);
        *f["hJets_chHEF_1"] = float(d["Jet_chHEF"][*in["hJetInd2"]]);
        *f["hJets_neHEF_1"] = float(d["Jet_neHEF"][*in["hJetInd2"]]);
        *f["hJets_neEmEF_1"] = float(d["Jet_neEmEF"][*in["hJetInd2"]]);
        *f["hJets_chMult_1"] = float(d["Jet_chMult"][*in["hJetInd2"]]);
        //*f["hJets_puId_1"] = float(d["hJets_puId"][1]);
        *f["hJets_vtx3DVal_1"] = 0.;
        *f["hJets_vtxNtracks_1"] = 0.;
        *f["hJets_vtxMass_1"] = float(d["Jet_vtxMass"][*in["hJetInd2"]]);
        *f["hJets_vtx3DVal_1"] = float(d["Jet_vtx3DVal"][*in["hJetInd2"]]);
        *f["hJets_vtxPt_1"] = float(d["Jet_vtxPt"][1]);
        *f["hJets_vtxNtracks_1"] = float(d["Jet_vtxNtracks"][*in["hJetInd2"]]);
        
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
        std::cout<<"Running Wenu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<d["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<d["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*d[\"met_pt\"] "<<*d["met_pt"]<<std::endl;
    }
    
    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts
    *in["elInd"] = -1;
    float elMaxPt = 0; // max pt of the electrons we select
    for (int i =0; i<*in["nselLeptons"]; i++) {
        if(fabs(in["selLeptons_pdgId"][i])==11 
            && d["selLeptons_pt"][i]      > *f["eptcut"] 
            && fabs(d["selLeptons_eta"][i]) < *f["eletacut"]
            && d["selLeptons_relIso03"][i]< *f["erelisocut"]
            && in["selLeptons_eleCutIdCSA14_25ns_v1"][i] >= *f["elidcut"]
            && *d["met_pt"]               > *f["metcut"]){
            *d["lepMetDPhi"]=fabs(EvalDeltaPhi(d["selLeptons_phi"][i],*d["met_phi"]));
            if (*d["lepMetDPhi"] < *f["elMetDPhiCut"] && d["selLeptons_pt"][i] > elMaxPt) {
                elMaxPt = d["selLeptons_pt"][i];
                *in["elInd"] = i;
            }
        }
    }
    if (*in["elInd"] == -1) return false;
        
    TLorentzVector W,El, MET, Hbb, HJ1, HJ2;
    // Reconstruct W
    MET.SetPtEtaPhiM(*d["met_pt"], 0., *d["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    El.SetPtEtaPhiM(d["selLeptons_pt"][*in["elInd"]], d["selLeptons_eta"][*in["elInd"]], d["selLeptons_phi"][*in["elInd"]], d["selLeptons_mass"][*in["elInd"]]); 
    W = MET + El; 
    *d["V_pt"] = W.Pt(); // uncomment this line if we want to recalculate W.pt ourselves
 
    // Reconstruct Higgs
    HJ1.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd1"]], d["Jet_eta"][*in["hJetInd1"]], d["Jet_phi"][*in["hJetInd1"]], d["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd2"]], d["Jet_eta"][*in["hJetInd2"]], d["Jet_phi"][*in["hJetInd2"]], d["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;
        
    // Now we can calculate whatever we want (transverse) with W and H four-vectors
    *d["HVdPhi"] = Hbb.DeltaPhi(W);

    if(*d["HVdPhi"]> *f["HVDPhiCut"] && *d["V_pt"] > *f["vptcut"]
        && *d["H_pt"] > *f["hptcut"]  ){
        selectEvent=true;
    }

    return selectEvent;
}

bool VHbbAnalysis::WmunuHbbSelection(){
    
    bool selectEvent=false;
    if(debug>1000){
        std::cout<<"Running Wmunu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<d["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<d["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*d[\"met_pt\"] "<<*d["met_pt"]<<std::endl;
    }
    
    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts
    *in["muInd"] = -1;
    float muMaxPt = 0; // max pt of the muons we select 
    for (int i=0; i<*in["nselLeptons"]; i++) {
        if(fabs(in["selLeptons_pdgId"][i])==13 
            && d["selLeptons_pt"][i]      > *f["muptcut"]
            && fabs(d["selLeptons_eta"][i]) < *f["muetacut"]
            && d["selLeptons_relIso03"][i]< *f["murelisocut"]
            && in["selLeptons_tightId"][i] >= *f["muidcut"]
            && *d["met_pt"]               > *f["metcut"]){
            *d["lepMetDPhi"]=fabs(EvalDeltaPhi(d["selLeptons_phi"][i],*d["met_phi"]));
            if (*d["lepMetDPhi"] < *f["muMetDPhiCut"] && d["selLeptons_pt"][i] > muMaxPt) {
                muMaxPt = d["selLeptons_pt"][i];
                *in["muInd"] = i;
            }
        }
    }

    if (*in["muInd"] == -1) return false;
    TLorentzVector W,Mu, MET, Hbb, HJ1, HJ2;
    // Reconstruct W
    MET.SetPtEtaPhiM(*d["met_pt"], 0., *d["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Mu.SetPtEtaPhiM(d["selLeptons_pt"][*in["muInd"]], d["selLeptons_eta"][*in["muInd"]], d["selLeptons_phi"][*in["muInd"]], d["selLeptons_mass"][*in["muInd"]]); 
    W = MET + Mu; 
    *d["V_pt"] = W.Pt(); // uncomment this line if we want to recalculate W.pt ourselves
 
    // Reconstruct Higgs
    HJ1.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd1"]], d["Jet_eta"][*in["hJetInd1"]], d["Jet_phi"][*in["hJetInd1"]], d["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd2"]], d["Jet_eta"][*in["hJetInd2"]], d["Jet_phi"][*in["hJetInd2"]], d["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;
        
    // Now we can calculate whatever we want (transverse) with W and H four-vectors
    *d["HVdPhi"] = Hbb.DeltaPhi(W);

    if(*d["HVdPhi"]> *f["HVDPhiCut"] && *d["V_pt"] > *f["vptcut"]
        && *d["H_pt"] > *f["hptcut"]  ){  
        selectEvent=true;
    }
    return selectEvent;
}

std::pair<int,int> VHbbAnalysis::HighestPtBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        /*if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j1ptCut"]
            && d["Jet_btagCSV"][i]>*f["j1ptCSV"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(d["Jet_pt"][pair.first]<d["Jet_pt"][i]){
                pair.first = i;
            }
        }*/
        if( pair.first == -1) {
            pair.first = i;
        }
        else if(d["Jet_pt"][pair.first]<d["Jet_pt"][i]){
            pair.first = i;
        }
    }
   
    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        /*if(in["Jet_puId"][i] == 1
            && d["Jet_pt"][i]>*f["j2ptCut"]
            && d["Jet_btagCSV"][i]>*f["j2ptCSV"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(d["Jet_pt"][pair.second]<d["Jet_pt"][i]){
                pair.second = i;
            }
        }*/
        if (pair.second == -1) {
            pair.second = i;
        }
        else if(d["Jet_pt"][pair.second]<d["Jet_pt"][i]){
            pair.second = i;
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

