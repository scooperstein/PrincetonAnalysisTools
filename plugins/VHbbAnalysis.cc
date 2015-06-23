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
    //return true; // for the moment don't impose any preselection
    bool sel=true;
    //if( *d["Vtype"]==3 ) sel=true;
    //if( *d["Vtype"]>=0 && *d["Vtype"]<=4) sel=true;
    // Preselect for two jets and one lepton which pass some minimum pt threshold
    int nPreselJets = 0;
    for (int i=0; i < *in["nJet"]; i++) {
        if (d["Jet_pt"][i] > *f["JetPtPresel"]) nPreselJets++;
    }
    int nPreselLep = 0;
    for (int i=0; i < *in["nselLeptons"]; i++) {
        if (d["selLeptons_pt"][i] > *f["LepPtPresel"]) nPreselLep++;
    }
    //for (int i=0; i < *in["naLeptons"]; i++) {
    //    if (d["aLeptons_pt"][i] > *f["LepPtPresel"]) nPreselLep++;
    //}
    if (nPreselJets < 2 || nPreselLep < 1) sel = false;
    return sel;
}

bool VHbbAnalysis::Analyze(){
    bool sel=true;
    bool doCutFlow = bool(*f["doCutFlow"]);
    *in["cutFlow"] = 0;

    if(debug>1000 && doCutFlow) {
        std::cout<<"Running cutflow"<<std::endl;
    }

    if(debug>1000) {
        std::cout<<"cutting on V and H pt"<<std::endl;
    }
    if(*d["V_pt"] < *f["vptcut"] || *d["H_pt"] < *f["hptcut"]) sel = false;
    if (sel) *in["cutFlow"] += 1;

    if(debug>1000) {
        std::cout<<"selecting bjets"<<std::endl;
    }
    //std::pair<int,int> bjets=HighestPtBJets();
    std::pair<int,int> bjets=HighestCSVBJets();
 
    if (bjets.first == -1) {
        sel = false;
        *in["hJetInd1"] = 0;
    }
    else *in["hJetInd1"] = bjets.first;
    if (bjets.second == -1) {
         sel = false;
        *in["hJetInd2"] = 1;
    }
    else *in["hJetInd2"] = bjets.second;

   /* *in["nHJetsMatched"] = 0;
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
    }*/

    /*// For now use higgs bjet selection from step 2 ntuples
    if(d["hJets_btagCSV"][0] < *f["j1ptCSV"] || d["hJets_btagCSV"][1] < *f["j2ptCSV"]) return false;
    if(d["hJets_pt"][0] < *f["j1ptCut"] || d["hJets_pt"][1] < *f["j2ptCut"]) return false;
    */

    /*// Cut on the bjets that we select
    if(d["Jet_btagCSV"][*in["hJetInd1"]] < *f["j1ptCSV"] || d["Jet_btagCSV"][*in["hJetInd2"]] < *f["j2ptCSV"]) sel = false;
    if(d["Jet_pt"][*in["hJetInd1"]] < *f["j1ptCut"] || d["Jet_pt"][*in["hJetInd2"]] < *f["j2ptCut"]) sel = false; 
    */
    if (sel) *in["cutFlow"] += 1; 

    if(debug>1000) {
        std::cout<<"nJet = "<<*in["nJet"]<<std::endl;
        std::cout<<"hJetInd1 = "<<*in["hJetInd1"]<<std::endl;
        std::cout<<"hJetInd2 = "<<*in["hJetInd2"]<<std::endl;
        std::cout<<"found two bjets with pt and CSV "
            <<d["Jet_pt"][*in["hJetInd1"]]<<" "
            <<d["Jet_btagCSV"][*in["hJetInd1"]]<<" "
            <<d["Jet_pt"][*in["hJetInd2"]]<<" "
            <<d["Jet_btagCSV"][*in["hJetInd2"]]<<" "
            <<std::endl;
    }
    
    if (*d["met_pt"] < *f["metcut"]) sel = false;
    if (sel) *in["cutFlow"] += 1;
    
    //check if event passes any event class
    *in["lepInd"] = -1;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    if(WmunuHbbSelection()) {
        *in["isWmunu"] = 1;
        *in["eventClass"]=0;
        *in["lepInd"] = *in["muInd"];
    } 
    if(WenuHbbSelection()) {
        *in["isWenu"] = 1;
        if (!(*in["isWmunu"])) {
            *in["eventClass"]=1000;
            *in["lepInd"] = *in["elInd"];
        }
    }
    if (*in["isWmunu"] == 0 && *in["isWenu"] == 0) sel = false;
    if (sel) *in["cutFlow"] += 1;
    
    if (*in["lepInd"] == -1) {
        // not Wenu or Wmunu, use preselected lepton
        *in["lepInd"] = 0;
    }

    if(debug>1000) {
        std::cout<<"cutting on dphi(lep, met)"<<std::endl;
    }
    *d["lepMetDPhi"]=fabs(EvalDeltaPhi(d["selLeptons_phi"][*in["lepInd"]],*d["met_phi"]));
    if (*in["isWmunu"] && *d["lepMetDPhi"] > *f["muMetDPhiCut"]) sel = false;
    else if (*in["isWenu"] && *d["lepMetDPhi"] > *f["elMetDPhiCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1;   


    TLorentzVector W,Lep, MET, Hbb, HJ1, HJ2, GenBJ1, GenBJ2, GenBJJ;
    // Reconstruct W
    MET.SetPtEtaPhiM(*d["met_pt"], 0., *d["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(d["selLeptons_pt"][*in["lepInd"]], d["selLeptons_eta"][*in["lepInd"]], d["selLeptons_phi"][*in["lepInd"]], d["selLeptons_mass"][*in["lepInd"]]); 
    W = MET + Lep; 
    *d["V_pt"] = W.Pt(); // uncomment this line if we want to recalculate W.pt ourselves
 
    // Reconstruct Higgs
    HJ1.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd1"]], d["Jet_eta"][*in["hJetInd1"]], d["Jet_phi"][*in["hJetInd1"]], d["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd2"]], d["Jet_eta"][*in["hJetInd2"]], d["Jet_phi"][*in["hJetInd2"]], d["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;

    // Compare gen kinematics for b jets for signal vs. ttbar
    if (*in["nGenBQuarkFromHafterISR"] > 1) {
        // signal event
        GenBJ1.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][0], d["GenBQuarkFromHafterISR_eta"][0], d["GenBQuarkFromHafterISR_phi"][0], d["GenBQuarkFromHafterISR_mass"][0]);
        GenBJ2.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][1], d["GenBQuarkFromHafterISR_eta"][1], d["GenBQuarkFromHafterISR_phi"][1], d["GenBQuarkFromHafterISR_mass"][1]);
    }        
    else if (*in["nGenBQuarkFromTop"] == 2){
        // ttbar event 
        GenBJ1.SetPtEtaPhiM(d["GenBQuarkFromTop_pt"][0], d["GenBQuarkFromTop_eta"][0], d["GenBQuarkFromTop_phi"][0], d["GenBQuarkFromTop_mass"][0]);
        GenBJ2.SetPtEtaPhiM(d["GenBQuarkFromTop_pt"][1], d["GenBQuarkFromTop_eta"][1], d["GenBQuarkFromTop_phi"][1], d["GenBQuarkFromTop_mass"][1]);
    }
    else {
       // let's not worry about these for now since ttbar is so dominant
       GenBJ1.SetPtEtaPhiM(0, 0, 0, 0);
       GenBJ2.SetPtEtaPhiM(0, 0, 0, 0);
    }

    *f["GenBJ1_pt"] = GenBJ1.Pt();
    *f["GenBJ1_eta"] = GenBJ1.Eta();
    *f["GenBJ1_phi"] = GenBJ1.Phi();
    *f["GenBJ1_mass"] = GenBJ1.M();
    *f["GenBJ2_pt"] = GenBJ2.Pt();
    *f["GenBJ2_eta"] = GenBJ2.Eta();
    *f["GenBJ2_phi"] = GenBJ2.Phi();
    *f["GenBJ2_mass"] = GenBJ2.M();

    GenBJJ = GenBJ1 + GenBJ2;
    *f["GenBJJ_pt"] = GenBJJ.Pt();
    *f["GenBJJ_eta"] = GenBJJ.Eta();
    *f["GenBJJ_phi"] = GenBJJ.Phi();
    *f["GenBJJ_mass"] = GenBJJ.M();
    *f["GenBJJ_dPhi"] = GenBJ2.DeltaPhi(GenBJ1);
    *f["GenBJJ_dR"] = GenBJ2.DeltaR(GenBJ1);
    double dEta = GenBJ2.Eta() - GenBJ1.Eta();
    // dEta should be on interval [-Pi,Pi]
    if (dEta > PI) {
        dEta -= 2*PI;
    }
    else if (dEta < -1*PI) {
        dEta += 2*PI;
    }
    *f["GenBJJ_dEta"] = dEta;
 
    // Now we can calculate whatever we want (transverse) with W and H four-vectors
    *d["HVdPhi"] = Hbb.DeltaPhi(W);
    *f["H_mass_step2"] = *d["H_mass"];
    *d["H_mass"] = Hbb.M(); // mass window cut? regression applied in FinishEvent
    *d["H_pt"] = Hbb.Pt();
 
    // Match Jets with Gen Higgs Jets
    for (int i=0; i<*in["nJet"]; i++) {
        TLorentzVector GenHJ1, GenHJ2, Jet;
        GenHJ1.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][0],d["GenBQuarkFromHafterISR_eta"][0],d["GenBQuarkFromHafterISR_phi"][0],d["GenBQuarkFromHafterISR_mass"][0]);
        GenHJ2.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][1],d["GenBQuarkFromHafterISR_eta"][1],d["GenBQuarkFromHafterISR_phi"][1],d["GenBQuarkFromHafterISR_mass"][1]);
        Jet.SetPtEtaPhiM(d["Jet_pt"][i], d["Jet_eta"][i], d["Jet_phi"][i], d["Jet_mass"][i]);
      
        double dR1 = Jet.DeltaR(GenHJ1);
        double dR2 = Jet.DeltaR(GenHJ2);
        d["Jet_genHJetMinDR"][i] = min(dR1, dR2);
        if(dR1 <= dR2) d["Jet_genHJetIndex"][i] = 1;
        else d["Jet_genHJetIndex"][i] = 2;
    }
    
    // count the number of additional leptons and jets, then cut on this number
    int nAddJet = 0;
    int nAddJet20 = 0;
    int nAddJet25 = 0;
    int nAddJet30 = 0;
    int nAddJet35 = 0;
    int nAddJet40 = 0;
    int nAddJet45 = 0;
    int nAddJet50 = 0;
    int nAddJet203p5 = 0;
    int nAddJet202p5 = 0;
    int nAddJet253p5 = 0;
    int nAddJet252p5 = 0;
    int nAddJet303p5 = 0;
    int nAddJet302p5 = 0;
    int nAddJet353p5 = 0;
    int nAddJet352p5 = 0;
    int nAddJet403p5 = 0;
    int nAddJet402p5 =0;
    int nAddJet453p5 = 0;
    int nAddJet452p5 =0;
    int nAddJet503p5 = 0;
    int nAddJet502p5 =0;
    int nAddJet_puid = 0;
    int nAddJet20_puid = 0;
    int nAddJet25_puid = 0;
    int nAddJet30_puid = 0;
    int nAddJet35_puid = 0;
    int nAddJet40_puid = 0;
    int nAddJet45_puid = 0;
    int nAddJet50_puid = 0;
    int nAddJet203p5_puid = 0;
    int nAddJet202p5_puid = 0;
    int nAddJet253p5_puid = 0;
    int nAddJet252p5_puid = 0;
    int nAddJet303p5_puid = 0;
    int nAddJet302p5_puid = 0;
    int nAddJet353p5_puid = 0;
    int nAddJet352p5_puid = 0;
    int nAddJet403p5_puid = 0;
    int nAddJet402p5_puid =0;
    int nAddJet453p5_puid = 0;
    int nAddJet452p5_puid =0;
    int nAddJet503p5_puid = 0;
    int nAddJet502p5_puid = 0;
    int nAddLep = 0;          
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet20++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet25++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet30++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet35++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet40++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet45++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_id"][i]>0) {
            nAddJet50++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet203p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet202p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet253p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet252p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet303p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet302p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet353p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet352p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet403p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet402p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet453p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet452p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_id"][i]>0) {
            nAddJet503p5++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_id"][i]>0) {
            nAddJet502p5++;
        }  
    }          

    
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet_puid++;
        }  
    }
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet20_puid++;
        }
    }
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet25_puid++;
        }
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet30_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet35_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet40_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet45_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<4.5 && in["Jet_puId"][i]>0) {
            nAddJet50_puid++;
        }  
    }       
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet203p5_puid++;
        }
    }
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>20 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet202p5_puid++;
        }
    }
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet253p5_puid++;
        }
    }
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>25 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet252p5_puid++;
        }
    }    
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet303p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>30 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet302p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet353p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>35 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet352p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet403p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>40 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet402p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet453p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>45 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet452p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<3.5 && in["Jet_puId"][i]>0) {
            nAddJet503p5_puid++;
        }  
    }           
    for(int i=0; i < *in["nJet"]; i++) {
        if(i == *in["hJetInd1"] || i == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][i]>50 && fabs(d["Jet_eta"][i])<2.5 && in["Jet_puId"][i]>0) {
            nAddJet502p5_puid++;
        }  
    }           
 
    *in["nAddJets"] = nAddJet;
    *in["nAddJets20"] = nAddJet20;
    *in["nAddJets25"] = nAddJet25;
    *in["nAddJets30"] = nAddJet30;
    *in["nAddJets35"] = nAddJet35;
    *in["nAddJets40"] = nAddJet40;
    *in["nAddJets45"] = nAddJet45;
    *in["nAddJets50"] = nAddJet50;
    *in["nAddJets203p5"] = nAddJet203p5;
    *in["nAddJets202p5"] = nAddJet202p5; 
    *in["nAddJets253p5"] = nAddJet253p5;
    *in["nAddJets252p5"] = nAddJet252p5; 
    *in["nAddJets303p5"] = nAddJet303p5;
    *in["nAddJets302p5"] = nAddJet302p5; 
    *in["nAddJets353p5"] = nAddJet353p5;
    *in["nAddJets352p5"] = nAddJet352p5;
    *in["nAddJets403p5"] = nAddJet403p5;
    *in["nAddJets402p5"] = nAddJet402p5;
    *in["nAddJets453p5"] = nAddJet453p5;
    *in["nAddJets452p5"] = nAddJet452p5;
    *in["nAddJets503p5"] = nAddJet503p5;
    *in["nAddJets502p5"] = nAddJet502p5;
    
    *in["nAddJets_puid"] = nAddJet_puid;
    *in["nAddJets20_puid"] = nAddJet20_puid;
    *in["nAddJets25_puid"] = nAddJet25_puid;
    *in["nAddJets30_puid"] = nAddJet30_puid;
    *in["nAddJets35_puid"] = nAddJet35_puid;
    *in["nAddJets40_puid"] = nAddJet40_puid;
    *in["nAddJets45_puid"] = nAddJet45_puid;
    *in["nAddJets50_puid"] = nAddJet50_puid;
    *in["nAddJets203p5_puid"] = nAddJet203p5_puid;
    *in["nAddJets202p5_puid"] = nAddJet202p5_puid; 
    *in["nAddJets253p5_puid"] = nAddJet253p5_puid;
    *in["nAddJets252p5_puid"] = nAddJet252p5_puid; 
    *in["nAddJets303p5_puid"] = nAddJet303p5_puid;
    *in["nAddJets302p5_puid"] = nAddJet302p5_puid; 
    *in["nAddJets353p5_puid"] = nAddJet353p5_puid;
    *in["nAddJets352p5_puid"] = nAddJet352p5_puid;
    *in["nAddJets403p5_puid"] = nAddJet403p5_puid;
    *in["nAddJets402p5_puid"] = nAddJet402p5_puid;
    *in["nAddJets453p5_puid"] = nAddJet453p5_puid;
    *in["nAddJets452p5_puid"] = nAddJet452p5_puid;
    *in["nAddJets503p5_puid"] = nAddJet503p5_puid;
    *in["nAddJets502p5_puid"] = nAddJet502p5_puid;

    // count additional leptons (check both collections, which are exclusive)
    for (int i=0; i<*in["nselLeptons"]; i++) {
        if (i == *in["lepInd"]) continue; // don't look at the lepton we've selected from the W
        if (d["selLeptons_pt"][i]>15 && fabs(d["selLeptons_eta"][i])<2.5 && d["selLeptons_relIso03"][i]<0.1) {
            nAddLep++;
        }
    }
    for (int i=0; i<*in["naLeptons"]; i++) {
        if (d["aLeptons_pt"][i]>15 && fabs(d["aLeptons_eta"][i])<2.5 && d["aLeptons_relIso03"][i]<0.1) {
            nAddLep++;
        }
    } 

    *in["nAddLeptons"] = nAddLep;
    if (nAddLep>= *f["nAddLeptonsCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1;

    if(nAddJet252p5_puid >= *f["nAddJetsCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1; 
    
    if(fabs(*d["HVdPhi"]) < *f["HVDPhiCut"]) sel = false;
     if (sel) *in["cutFlow"] += 1;
 
    if(debug>1000) {
        std::cout<<"selecting event"<<std::endl;
    }

    if (doCutFlow) return true; // keep all preselected events for cutflow
    else return sel;
}

void VHbbAnalysis::FinishEvent(){
   
    //if (bool(*f["doCutFlow"])) {
    //    ofile->cd();
    //    outputTree->Fill();
    //    return;
    //} 
    // General use variables
    *in["sampleIndex"] = cursample->sampleNum;
   
    // Split WJets and ZJets samples by jet parton flavor
    if (cursample->doJetFlavorSplit) {
        int nBJets = 0;
        //std::cout<<fabs(d["hJets_mcFlavour"][0])<<std::endl;
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd1"]]) == 5) nBJets++;
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd2"]]) == 5) nBJets++;
        *in["sampleIndex"] = (*in["sampleIndex"]*100 + nBJets);
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

    // compare our Higgs jet selection with step 2 
    *d["hJets_pt_0_step2"] = d["hJets_pt"][0];
    *d["hJets_pt_1_step2"] = d["hJets_pt"][1];
    *d["hJets_btagCSV_0_step2"] = d["hJets_btagCSV"][0];
    *d["hJets_btagCSV_1_step2"] = d["hJets_btagCSV"][1];

    *f["weight"] = cursample->intWeight;
    //*f["weight"] = 1.0; // HACK FIXME
    *f["Vtype_f"] = (float) *d["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["selLeptons_pt_0"] = d["selLeptons_pt"][*in["lepInd"]];
    *f["selLeptons_eta_0"] = d["selLeptons_eta"][*in["lepInd"]];
    *f["selLeptons_phi_0"] = d["selLeptons_phi"][*in["lepInd"]];
    *in["selLeptons_pdgId_0"] = in["selLeptons_pdgId"][*in["lepInd"]];    
    *in["selLeptons_eleCutIdCSA14_25ns_v1_0"] = in["selLeptons_eleCutIdCSA14_25ns_v1"][*in["lepInd"]];
    *in["selLeptons_tightId_0"] = in["selLeptons_tightId"][*in["lepInd"]];
    *f["selLeptons_relIso03_0"] = d["selLeptons_relIso03"][*in["lepInd"]];  
  
    // electron ID variables
    *f["selLeptons_eleSieie_0"] = d["selLeptons_eleSieie"][*in["lepInd"]];
    *f["selLeptons_eleDEta_0"] = d["selLeptons_eleDEta"][*in["lepInd"]];
    *f["selLeptons_eleDPhi_0"] = d["selLeptons_eleDPhi"][*in["lepInd"]];
    *f["selLeptons_eleHoE_0"] = d["selLeptons_eleHoE"][*in["lepInd"]];
    *f["selLeptons_eleMissingHits_0"] = d["selLeptons_eleMissingHits"][*in["lepInd"]];
    *f["selLeptons_eleChi2_0"] = d["selLeptons_eleChi2"][*in["lepInd"]];
     
 
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
 

    *f["naLeptonsPassingCuts"] = *in["nAddLeptons"]; // no need to calculate this twice
    *f["naJetsPassingCuts"] = *in["nAddJets"];
    *f["highestCSVaJet"] = 0.0;
    *f["minDeltaRaJet"] = 9999.0; // FIXME add correct initialization
    for(int j=0;j<*in["nJet"];j++){
        if (j == *in["hJetInd1"] || j == *in["hJetInd2"]) continue;
        if(d["Jet_pt"][j]>20 && fabs(d["Jet_eta"][j])<2.5 && d["Jet_btagCSV"][j]>*f["highestCSVaJet"]) *f["highestCSVaJet"] = d["Jet_btagCSV"][j];
      }       
    
    // Reconstruct Higgs and W and recalculate variables ourselves
    TLorentzVector MET,Lep,W,HJ1,HJ2,Hbb;
    MET.SetPtEtaPhiM(*d["met_pt"], 0., *d["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(d["selLeptons_pt"][*in["lepInd"]], d["selLeptons_eta"][*in["lepInd"]], d["selLeptons_phi"][*in["lepInd"]], d["selLeptons_mass"][*in["lepInd"]]);
    W = MET + Lep;

    HJ1.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd1"]], d["Jet_eta"][*in["hJetInd1"]], d["Jet_phi"][*in["hJetInd1"]], d["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(d["Jet_pt"][*in["hJetInd2"]], d["Jet_eta"][*in["hJetInd2"]], d["Jet_phi"][*in["hJetInd2"]], d["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;
    
    // We already calculate these in Analyze()
    //*d["H_mass"] = Hbb.M();
    //*d["H_pt"] = Hbb.Pt();
    //*d["V_pt"] = W.Pt();
    //*d["HVdPhi"] = Hbb.DeltaPhi(W);
    
    // Set variables used by the BDT
    *f["H_mass_f"] = (float) *d["H_mass"];
    *f["H_pt_f"] = (float) *d["H_pt"];
    *f["V_pt_f"] = (float) *d["V_pt"];
    *f["hJets_btagCSV_0"] = (float) d["Jet_btagCSV"][*in["hJetInd1"]];
    *f["hJets_btagCSV_1"] = (float) d["Jet_btagCSV"][*in["hJetInd2"]]; 
    *f["HVdPhi_f"] = (float) *d["HVdPhi"];
    *f["H_dEta"] = fabs(d["Jet_eta"][*in["hJetInd1"]] - d["Jet_eta"][*in["hJetInd2"]]);
   
    *f["hJets_mt_0"] = HJ1.Mt();
    *f["hJets_mt_1"] = HJ2.Mt();
    *f["H_dR"] = (float) HJ1.DeltaR(HJ2);   
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
    bool selectEvent=true;
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
            ){
            if (d["selLeptons_pt"][i] > elMaxPt) {
                elMaxPt = d["selLeptons_pt"][i];
                *in["elInd"] = i;
            }
        }
    }
    if (*in["elInd"] == -1) selectEvent = false;

    return selectEvent;
}

bool VHbbAnalysis::WmunuHbbSelection(){
    
    bool selectEvent=true;
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
            ){
            if (d["selLeptons_pt"][i] > muMaxPt) {
                muMaxPt = d["selLeptons_pt"][i];
                *in["muInd"] = i;
            }
        }
    }
    if (*in["muInd"] == -1) selectEvent = false;

    return selectEvent;
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
        /*if( pair.first == -1) {
            pair.first = i;
        }
        else if(d["Jet_pt"][pair.first]<d["Jet_pt"][i]){
            pair.first = i;
        }*/
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
        /*if (pair.second == -1) {
            pair.second = i;
        }
        else if(d["Jet_pt"][pair.second]<d["Jet_pt"][i]){
            pair.second = i;
        }*/
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

