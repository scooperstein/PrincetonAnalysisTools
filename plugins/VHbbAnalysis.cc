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
    if (*f["onlyEvenEvents"] && *f["onlyOddEvents"]) {
        std::cout<<"Cannot set both onlyEvenEvents and onlyOddEvents to true!!"<<std::endl;
        return false;
    }
    if (*f["onlyEvenEvents"]) {
        if (*in["evt"]%2 == 1) sel=false;
    }    
    else if (*f["onlyOddEvents"]) {
        if (*in["evt"]%2 == 0) sel=false;
    }    
    // stitch WJets inclusive sample to HT-binned samples
    if (cursample->sampleNum == 22 && *f["lheHT"] > 100) sel=false; 

    // Preselect for two jets and one lepton which pass some minimum pt threshold
    int nPreselJets = 0;
    for (int i=0; i < *in["nJet"]; i++) {
        if (f["Jet_pt"][i] > *f["JetPtPresel"]) nPreselJets++;
        //std::cout<<"Jet_pt["<<i<<"] = "<<f["Jet_pt"][i]<<std::endl;
    }
    int nPreselLep = 0;
    for (int i=0; i < *in["nselLeptons"]; i++) {
        if (f["selLeptons_pt"][i] > *f["LepPtPresel"]) nPreselLep++;
    }
    //for (int i=0; i < *in["naLeptons"]; i++) {
    //    if (f["aLeptons_pt"][i] > *f["LepPtPresel"]) nPreselLep++;
    //}
    if (nPreselJets < 2 || nPreselLep < 1) sel = false;
    return sel;
}

bool VHbbAnalysis::Analyze(){
    *in["sampleIndex"] = cursample->sampleNum;
    bool sel=true;
    bool doCutFlow = bool(*f["doCutFlow"]);
    *in["cutFlow"] = 0;
    

    if(debug>1000) {
         std::cout<<"Imposing trigger and json requirements"<<std::endl;
    }
    // Impose trigger requirements
    //if (*in["HLT_WenHbbLowLumi"]!=1 && *in["HLT_WmnHbbLowLumi"]!=1) sel=false;
    if (*in["HLT_BIT_HLT_Ele23_WPLoose_Gsf_v"]!=1 && *in["HLT_BIT_HLT_Ele27_WPLoose_Gsf_v"] && *in["HLT_BIT_HLT_IsoMu20_v"]!=1 && *in["HLT_BIT_HLT_IsoTkMu20_v"]!=1 ) sel=false;
    if (*in["sampleIndex"]==0) {
        if (*f["json"]!=1) sel=false;
    }
    if (sel) *in["cutFlow"]+=1;   
 

    // let's also keep track of each jet selection method separately, for science!
    std::pair<int,int> bjets_bestCSV = HighestCSVBJets();
    *in["hJetInd1_bestCSV"] = bjets_bestCSV.first;
    *in["hJetInd2_bestCSV"] = bjets_bestCSV.second; 
    std::pair<int,int> bjets_highestPt = HighestPtBJets();
    *in["hJetInd1_highestPt"] = bjets_highestPt.first;
    *in["hJetInd2_highestPt"] = bjets_highestPt.second;
    std::pair<int,int> bjets_highestPtJJ = HighestPtJJBJets();
    *in["hJetInd1_highestPtJJ"] = bjets_highestPtJJ.first;
    *in["hJetInd2_highestPtJJ"] = bjets_highestPtJJ.second; 

    // the jet selection algorithm we actually use for the rest of the analysis chain
    //std::pair<int,int> bjets=HighestPtBJets();
    std::pair<int,int> bjets=HighestCSVBJets();
    //std::pair<int,int> bjets=HighestPtJJBJets(); 

    // put CSV cuts out of selection functions
    if(bjets.first != -1 && bjets.second != -1){
        if(f["Jet_btagCSV"][bjets.first]<*f["j1ptCSV"] || f["Jet_btagCSV"][bjets.second]<*f["j2ptCSV"]){
            sel=false;
        }
    }

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
  
    if (sel) *in["cutFlow"] += 1; // selected jets

    if(debug>1000) {
        std::cout<<"nJet = "<<*in["nJet"]<<std::endl;
        std::cout<<"hJetInd1 = "<<*in["hJetInd1"]<<std::endl;
        std::cout<<"hJetInd2 = "<<*in["hJetInd2"]<<std::endl;
        std::cout<<"found two bjets with pt and CSV "
            <<f["Jet_pt_reg"][*in["hJetInd1"]]<<" "
            <<f["Jet_btagCSV"][*in["hJetInd1"]]<<" "
            <<f["Jet_pt_reg"][*in["hJetInd2"]]<<" "
            <<f["Jet_btagCSV"][*in["hJetInd2"]]<<" "
            <<std::endl;
    }
    
    // Reconstruct Higgs
    TLorentzVector HJ1,HJ2,Hbb;
    HJ1.SetPtEtaPhiM(f["Jet_pt_reg"][*in["hJetInd1"]], f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(f["Jet_pt_reg"][*in["hJetInd2"]], f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;

    *f["H_pt"] = Hbb.Pt();
    if (*f["H_pt"] < *f["hptcut"]) sel=false;
    if (sel) *in["cutFlow"] += 1; // pT(jj) cut
   
    if (*f["met_pt"] < *f["metcut"]) sel = false;
    if (sel) *in["cutFlow"] += 1; // met cut
    
    //check if event passes any event class
    *in["lepInd"] = -1;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    if(MuonSelection()) {
        *in["isWmunu"] = 1;
        *in["eventClass"]=0;
        *in["lepInd"] = *in["muInd"];
    } 
    if(ElectronSelection()) {
        *in["isWenu"] = 1;
        if (!(*in["isWmunu"])) {
            *in["eventClass"]=1000;
            *in["lepInd"] = *in["elInd"];
        }
    }
    if (*in["isWmunu"] == 0 && *in["isWenu"] == 0) sel = false;
    if (sel) *in["cutFlow"] += 1; // lepton selection
    
    
    if (*in["lepInd"] == -1) {
        // not Wenu or Wmunu, use preselected lepton
        *in["lepInd"] = 0;
    }

    if(debug>1000) {
        std::cout<<"cutting on dphi(lep, met)"<<std::endl;
        std::cout<<*in["lepInd"]<<" "<<f["selLeptons_phi"][*in["lepInd"]]<<" "<<f["selLeptons_eta"][*in["lepInd"]]<<" "<<f["selLeptons_mass"][*in["lepInd"]]<<std::endl;
    }
    *f["lepMetDPhi"]=fabs(EvalDeltaPhi(f["selLeptons_phi"][*in["lepInd"]],*f["met_phi"]));
    if (*in["isWmunu"] && *f["lepMetDPhi"] > *f["muMetDPhiCut"]) sel = false;
    else if (*in["isWenu"] && *f["lepMetDPhi"] > *f["elMetDPhiCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1;   


    TLorentzVector W,Lep, MET, GenBJ1, GenBJ2, GenBJJ;
    // Reconstruct W
    if(debug>1000) {
        std::cout<<"met "<<*f["met_pt"]<<" "<<*f["met_phi"]<<std::endl;
    }
    MET.SetPtEtaPhiM(*f["met_pt"], 0., *f["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], f["selLeptons_mass"][*in["lepInd"]]); 
    W = MET + Lep; 
    *f["V_pt"] = W.Pt(); // uncomment this line if we want to recalculate W.pt ourselves


    *f["Lep_HJ1_dPhi"] = Lep.DeltaPhi(HJ1);
    *f["Lep_HJ2_dPhi"] = Lep.DeltaPhi(HJ2);
    *f["HJ1_HJ2_dPhi"] = HJ1.DeltaPhi(HJ2);

    *f["Top1_mass_fromLepton"] = GetRecoTopMass(Lep, false, 0, false); // construct top mass from closest jet to lepton
    *f["Top1_mass_fromLepton_wMET"] = GetRecoTopMass(Lep, false, 1, false); // construct top mass from closest jet to lepton
    
    *f["Top1_mass_fromLepton_regPT"] = GetRecoTopMass(Lep, false, 0, true); // construct top mass from closest jet to lepton
    *f["Top1_mass_fromLepton_regPT_wMET"] = GetRecoTopMass(Lep, false, 1, true); // construct top mass from closest jet to lepton

    // Let's try to reconstruct the second top from a hadronically decaying W
    double min1 = 999.;
    double min2 = 999.;
    int jetInd1 = -1;
    int jetInd2 = -1;
    for (int i=0; i<*in["nJet"]; i++) {
        if (i==*in["hJetInd1"] || i==*in["hJetInd2"]) continue;
        // make some basic jet id/kinematic cut
        if (f["Jet_pt_reg"][i] < 30 || in["Jet_puId"][i]==0) continue;
        //if (in["Jet_mcMatchId"][i] == 0) continue;
        TLorentzVector j;
        j.SetPtEtaPhiM(f["Jet_pt_reg"][i], f["Jet_eta"][i], f["Jet_phi"][i], f["Jet_mass"][i] );
        double DR = j.DeltaR(HJ2);
        if (DR <= min1) {
            min2 = min1;
            jetInd2 = jetInd1;
            min1 = DR;
            jetInd1 = i;
        }
        else if (DR <= min2) {
            min2 = DR;
            jetInd2 = i;
        }
    }
    *f["HJ2_minJetDR1"] = min1;
    *f["HJ2_minJetDR2"] = min2;
    if (jetInd1!=-1 && jetInd2!=-1) {
        TLorentzVector WJet1, WJet2;
        WJet1.SetPtEtaPhiM(f["Jet_pt_reg"][jetInd1], f["Jet_eta"][jetInd1], f["Jet_phi"][jetInd1], f["Jet_mass"][jetInd1]);
        WJet2.SetPtEtaPhiM(f["Jet_pt_reg"][jetInd2], f["Jet_eta"][jetInd2], f["Jet_phi"][jetInd2], f["Jet_mass"][jetInd2]);
        TLorentzVector Top_hadW = WJet1 + WJet2 + HJ2;
        *f["Top2_mass_hadW"] = Top_hadW.M();
        *f["HJ2_WJet1_dPhi"] = HJ2.DeltaPhi(WJet1);
        *f["HJ2_WJet2_dPhi"] = HJ2.DeltaPhi(WJet2);
        *f["HJ2_WJet1_dEta"] = fabs(HJ2.Eta() - WJet1.Eta()); 
        *f["HJ2_WJet2_dEta"] = fabs(HJ2.Eta() - WJet2.Eta()); 
    }
    else {
        *f["Top2_mass_hadW"] = -999;
        *f["HJ2_WJet1_dPhi"] = -999;
        *f["HJ2_WJet2_dPhi"] = -999;
        *f["HJ2_WJet1_dEta"] = -999;
        *f["HJ2_WJet2_dEta"] = -999;
    }

 
    // Now we can calculate whatever we want (transverse) with W and H four-vectors
    *f["HVdPhi"] = Hbb.DeltaPhi(W);
    *f["H_mass_step2"] = *f["H_mass"];
        
    *f["H_mass"] = Hbb.M(); // mass window cut? regression applied in FinishEvent
    if (cursyst->name != "nominal") {
        *f[Form("H_mass_%s",cursyst->name.c_str())] = Hbb.M();
    }
    //*f["H_pt"] = Hbb.Pt(); // we already do this


    // Compare gen kinematics for b jets for signal vs. ttbar
    if(*in["sampleIndex"]!=0) {
        if (*in["nGenBQuarkFromH"] > 1) {
            // signal event
            GenBJ1.SetPtEtaPhiM(f["GenBQuarkFromH_pt"][0], f["GenBQuarkFromH_eta"][0], f["GenBQuarkFromH_phi"][0], f["GenBQuarkFromH_mass"][0]);
            GenBJ2.SetPtEtaPhiM(f["GenBQuarkFromH_pt"][1], f["GenBQuarkFromH_eta"][1], f["GenBQuarkFromH_phi"][1], f["GenBQuarkFromH_mass"][1]);
        }        
        else if (*in["nGenBQuarkFromTop"] > 0){
            GenBJ1.SetPtEtaPhiM(f["GenBQuarkFromTop_pt"][0], f["GenBQuarkFromTop_eta"][0], f["GenBQuarkFromTop_phi"][0], f["GenBQuarkFromTop_mass"][0]);
            if (*in["nGenBQuarkFromTop"] > 1) {
                GenBJ2.SetPtEtaPhiM(f["GenBQuarkFromTop_pt"][1], f["GenBQuarkFromTop_eta"][1], f["GenBQuarkFromTop_phi"][1], f["GenBQuarkFromTop_mass"][1]);
            }
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
        *f["GenBJJ_dEta"] = fabs(GenBJ1.Eta() - GenBJ2.Eta());

        TLorentzVector GenLep1, GenLep2; // closest gen lep to either jet1 or jet2. Sometimes these could be the same lepton.
        double minDR1 = 999;
        double minDR2 = 999;
        double minSecDR2 = 999; //second closest gen lepton to second jet
        int GenLepIndex1 = -1; // index of the lepton closest to jet 1
        int GenLepIndex2 = -1; // index of the lepton closest to jet 2
        int GenLepSecIndex2 = -1; // index of the lepton second closest to jet 2
        for (int i=0; i<*in["nGenLep"]; i++) {
           TLorentzVector gl;
           gl.SetPtEtaPhiM(f["GenLep_pt"][i], f["GenLep_eta"][i], f["GenLep_phi"][i], f["GenLep_mass"][i] );
           double DR1 = gl.DeltaR(GenBJ1);
           double DR2 = gl.DeltaR(GenBJ2);
           if (DR1 <= minDR1) {
               minDR1 = DR1;
               GenLepIndex1 = i;
           }
           if (DR2 <= minDR2) {
               minSecDR2 = minDR2;
               GenLepSecIndex2 = GenLepIndex2;
               minDR2 = DR2;
               GenLepIndex2 = i;
           } 
           else if (DR2 < minSecDR2) {
               minSecDR2 = DR2;
               GenLepSecIndex2 = i;
           }
        }  
        if (GenLepIndex1 == GenLepIndex2) {
            // don't allow us to use the same lepton for each jet
            GenLepIndex2 = GenLepSecIndex2;
        }
            
        *in["GenLepIndex1"] = GenLepIndex1;
        *in["GenLepIndex2"] = GenLepIndex2;
        
        if (GenLepIndex1!=-1) {
            GenLep1.SetPtEtaPhiM(f["GenLep_pt"][GenLepIndex1], f["GenLep_eta"][GenLepIndex1], f["GenLep_phi"][GenLepIndex1], f["GenLep_mass"][GenLepIndex1] );
            *f["GenLep_GenBJ1_dR"] = GenLep1.DeltaR(GenBJ1);
            *f["GenLep_GenBJ1_dEta"] = fabs(GenLep1.Eta() - GenBJ1.Eta());
            *f["GenLep_GenBJ1_dPhi"] = GenLep1.DeltaPhi(GenBJ1);
            
            // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
            TLorentzVector GenTop1 = GenLep1 + GenBJ1;
            *f["GenTop1_mass"] = GenTop1.M();
        }
        else { 
            *f["GenLep_GenBJ1_dR"] = -999;
            *f["GenLep_GenBJ1_dEta"] = -999;
            *f["GenLep_GenBJ1_dPhi"] = -999;
            *f["GenTop1_mass"] = -999;
        }
        if (GenLepIndex2!=-1) {
            GenLep2.SetPtEtaPhiM(f["GenLep_pt"][GenLepIndex2], f["GenLep_eta"][GenLepIndex2], f["GenLep_phi"][GenLepIndex2], f["GenLep_mass"][GenLepIndex2] );
            *f["GenLep_GenBJ2_dR"] = GenLep2.DeltaR(GenBJ2);
            *f["GenLep_GenBJ2_dEta"] = (GenLep2.Eta(), GenBJ2.Eta());
            *f["GenLep_GenBJ2_dPhi"] = GenLep2.DeltaPhi(GenBJ2);

            // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
            TLorentzVector GenTop2 = GenLep2 + GenBJ2; 
            *f["GenTop2_mass"] = GenTop2.M();
        } 
        else {
            *f["GenLep_GenBJ2_dR"] = -999;
            *f["GenLep_GenBJ2_dEta"] = -999;
            *f["GenLep_GenBJ2_dPhi"] = -999;
            *f["GenTop2_mass"] = -999;
        }
    
        // construct Gen W
        TLorentzVector GenW1, GenW2;
        if (*in["nGenVbosons"]>0 && fabs(in["GenVbosons_pdgId"][0])==24) {
            GenW1.SetPtEtaPhiM(f["GenVbosons_pt"][0], f["GenVbosons_eta"][0], f["GenVbosons_phi"][0], f["GenVbosons_mass"][0]);
            *f["GenW_GenBJJ_dPhi"] = GenW1.DeltaPhi(GenBJJ);
            *f["GenW_GenBJJ_dEta"] = fabs(GenW1.Eta() - GenBJJ.Eta());
            // grab both W's in ttbar events
            if (*in["nGenVbosons"]>1 && fabs(in["GenVbosons_pdgId"][1])==24) {
                GenW2.SetPtEtaPhiM(f["GenVbosons_pt"][1], f["GenVbosons_eta"][1], f["GenVbosons_phi"][1], f["GenVbosons_mass"][1]);
            }
        }
        else {
        
            *f["GenW_GenBJJ_dPhi"] = -999;
            *f["GenW_GenBJJ_dEta"] = -999;
        }
  
        std::vector<TLorentzVector> genWQuarks; // gen quarks from hadronic gen W decay
        for (int i=0; i<*in["nGenWZQuark"]; i++) {
            TLorentzVector v;
            v.SetPtEtaPhiM(f["GenWZQuark_pt"][i], f["GenWZQuark_eta"][i], f["GenWZQuark_phi"][i], f["GenWZQuark_mass"][i]);
            genWQuarks.push_back(v);
        }
 
        //int nSelectedJetsMatched = 0; // count the number (0, 1, 2) of selected jets matched to the real bottom quarks
        // Match Jets with Gen B Jets from Higgs/Tops
        for (int i=0; i<*in["nJet"]; i++) {
            in["Jet_genJetMatchId"][i] = 0; // 0 if no gen match, 1 for pt-leading b-jet, 2 for pt sub-leading b-jet, 3 if matched to jet from hadronic W decay
            TLorentzVector Jet;
            //GenHJ1.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][0],d["GenBQuarkFromHafterISR_eta"][0],d["GenBQuarkFromHafterISR_phi"][0],d["GenBQuarkFromHafterISR_mass"][0]);
            //GenHJ2.SetPtEtaPhiM(d["GenBQuarkFromHafterISR_pt"][1],d["GenBQuarkFromHafterISR_eta"][1],d["GenBQuarkFromHafterISR_phi"][1],d["GenBQuarkFromHafterISR_mass"][1]);
            Jet.SetPtEtaPhiM(f["Jet_pt_reg"][i], f["Jet_eta"][i], f["Jet_phi"][i], f["Jet_mass"][i]);
          
            //double dR1 = Jet.DeltaR(GenHJ1);
            //double dR2 = Jet.DeltaR(GenHJ2);
            double dR1 = 999;
            if (GenBJ1.Pt() > 0) dR1 = Jet.DeltaR(GenBJ1);
            double dR2 = 999;
            if (GenBJ2.Pt() > 0) dR2 = Jet.DeltaR(GenBJ2);
            
            // try to match the jet to one of the jets from hadronic W decay
            double dR3 = 999;
            for (int j=0; j<(int)genWQuarks.size(); j++) {
                double Jet_genWQuarkDR = Jet.DeltaR(genWQuarks[j]);
                if (Jet_genWQuarkDR < dR3) {
                    dR3 = Jet_genWQuarkDR;
                }
            }
            f["Jet_genWQuarkDR"][i] = dR3;
            if (dR3 < min(dR1, dR2) && dR3 < 0.5) {
                in["Jet_genJetMatchId"][i] = 3; 
            }
            else if(dR1 <= dR2 && dR1<0.5) in["Jet_genJetMatchId"][i] = 1;
            else if (dR2 < 0.5) in["Jet_genJetMatchId"][i] = 2;
            
            f["Jet_genHJetMinDR"][i] = min(dR1, dR2);

            if (i == *in["hJetInd1"]) {
                *f["hJet1_matchedMinDR"] = f["Jet_genHJetMinDR"][i];
            }
            else if (i == *in["hJetInd2"]) {
                *f["hJet2_matchedMinDR"] = f["Jet_genHJetMinDR"][i];
            }
        }
    }

    if(debug>1000) std::cout<<"counting additional leptons"<<std::endl;
    // count the number of additional leptons and jets, then cut on this number
    int nAddLep = 0;       
    // 15 to 30 by 1 GeV, 1.5 to 3 w/ 0.1 in eta 
    //std::vector<int> ptCuts = {15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
    //std::vector<double> etaCuts = {1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5};
    std::vector<int> ptCuts = {25};
    std::vector<double> etaCuts = {2.9};

    for (int i=0; i < (int)ptCuts.size(); i++) {
        for (int j=0; j < (int)etaCuts.size(); j++) {
            int nAddJet_tmp = 0;
            for (int k=0; k < *in["nJet"]; k++) {
                if (k == *in["hJetInd1"] || k == *in["hJetInd2"]) continue;
                if (f["Jet_pt"][k]>ptCuts[i] && fabs(f["Jet_eta"][k])<etaCuts[j] && in["Jet_puId"][k]>0) {
                    nAddJet_tmp++;
                }
            }
            std::string eta_cut = Form("%.1f", etaCuts[j]); // convert, say, 2.5 to '2p5' 
            std::replace(eta_cut.begin(), eta_cut.end(), '.', 'p');
            //std::cout<<Form("nAddJets%i%s_puid",ptCuts[i],eta_cut.c_str())<<std::endl;
            *in[Form("nAddJets%i%s_puid",ptCuts[i],eta_cut.c_str())] = nAddJet_tmp;
        }
    }
    
    // count additional leptons (check both collections, which are exclusive)
    for (int i=0; i<*in["nselLeptons"]; i++) {
        if (i == *in["lepInd"]) continue; // don't look at the lepton we've selected from the W
        if (f["selLeptons_pt"][i]>15 && fabs(f["selLeptons_eta"][i])<2.5 && f["selLeptons_relIso03"][i]<0.1) {
            nAddLep++;
        }
    }
    for (int i=0; i<*in["naLeptons"]; i++) {
        if (f["aLeptons_pt"][i]>15 && fabs(f["aLeptons_eta"][i])<2.5 && f["aLeptons_relIso03"][i]<0.1) {
            nAddLep++;
        }
    } 

    *in["nAddLeptons"] = nAddLep;
    if (nAddLep>= *f["nAddLeptonsCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1; // additional lepton veto

    if(*in["nAddJets252p9_puid"] >= *f["nAddJetsCut"]) sel = false;
    if (sel) *in["cutFlow"] += 1; // additional jet veto
    
    if(fabs(*f["HVdPhi"]) < *f["HVDPhiCut"]) sel = false;
     if (sel) *in["cutFlow"] += 1; // dPhi(jj,W) cut

    // let's do this for now in WorkspaceAnalysis since it's fast to rerun there
    // mass window, use parameterized mean/sigma from jet mass res. study
    double mean = 104.7 + 0.055 * (*f["H_pt"]);
    double sigma = 20.9* exp(-0.0024*(*f["H_pt"]));

    /*// mass window cut
    if (*f["H_mass"] < (mean - 2*sigma) || *f["H_mass"] > (mean + 2*sigma)) sel=false;
    //if (*d["H_mass"] < 100 || *d["H_mass"] > 150) sel=false;
    if (sel) *in["cutFlow"] += 1; 
    */

    // Control samples

    // if cutflow==0 then jets are bogus
    // jet pt
    // lepton pt, iso, id
    // met_pt > threshold
    // deltaPhi(met,lep)
    // V_pt > 100 and bb_mass<250
    bool baseCSSelection= 
            (*in["cutFlow"]!=0
            && f["Jet_pt_reg"][*in["hJetInd1"]]>*f["j1ptCut"]
            && f["Jet_pt_reg"][*in["hJetInd2"]]>*f["j2ptCut"]
            && (*in["isWmunu"] != 0 || *in["isWenu"] != 0)
            && *f["met_pt"] > *f["metcut"]
            && ((*in["isWmunu"] && *f["lepMetDPhi"] < *f["muMetDPhiCut"])
              || (*in["isWenu"] && *f["lepMetDPhi"] < *f["elMetDPhiCut"]))
            && *f["V_pt"]> *f["vptcut"] && *f["H_mass"]<250 && *f["H_pt"] > *f["hptcut"]);
    
    *in["controlSample"]=-1; // maybe this should go much early (before any returns)
    float maxCSV=std::max(f["Jet_btagCSV"][*in["hJetInd1"]],f["Jet_btagCSV"][*in["hJetInd2"]]);
    if(baseCSSelection){
        *in["controlSample"]=0;
        if (maxCSV > 0.97){ //ttbar or W+HF
            if(*in["nAddJets252p9_puid"]>1.5) { //ttbar
                *in["controlSample"]=1;
            //} else if (*in["nAddJets252p9_puid"]<1.5 && *f["met_pt"]/ sqrt(*f["met_sumEt"])> 2. && (*f["H_mass"]>150 || *f["H_mass"]<90) ){ //W+HF
            } else if (*in["nAddJets252p9_puid"]<1.5 && (*f["H_mass"]>150 || *f["H_mass"]<90) ){ //W+HF
                *in["controlSample"]=2;
            }
        //} else if (maxCSV > 0.3 && *f["met_pt"]/ sqrt(*f["met_sumEt"]) > 2.){ //W+LF
        } else if (maxCSV > 0.3 && maxCSV<0.8){ //W+LF
            *in["controlSample"]=3;
        }
        if(*in["sampleIndex"]==0){
            std::cout<<"data CS event "<<*in["controlSample"]<<" maxCSV "<<maxCSV<<" nAddJets252p9_puid "<<*in["nAddJets252p9_puid"]<<" met_pt "<<*f["met_pt"]<<" met_sumEt "<<*f["met_sumEt"]<<" H_mass "<<*f["H_mass"]<<std::endl;
        }
    }
    
    if(*f["doControlSamples"]==1){
        //if(!sel) sel=( *in["controlSample"]>-1 ); 
        sel=( *in["controlSample"]>-1 ); 
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
    *f["weight_ptQCD"] = 1.0;
    *f["weight_ptEWK"] = 1.0;
    if(*in["sampleIndex"]!=0){
        *f["weight_PU"]=ReWeightMC(*in["nTrueInt"]);
        if (*in["nGenTop"]==0 && *in["nGenVbosons"]>0) {
            // only apply to Z/W+jet samples
            *f["weight_ptQCD"]=ptWeightQCD(*in["nGenVbosons"], *f["lheHT"], in["GenVbosons_pdgId"][0]);
            *f["weight_ptEWK"]=ptWeightEWK(*in["nGenVbosons"], f["GenVbosons_pt"][0], *f["VtypeSim"], in["GenVbosons_pdgId"][0]);
        }
    } else {
        *f["weight_PU"]=1;
        *f["puWeight"]=1;
    }

    //*f["weight"]= *f["weight"] * *f["weight_PU"];

    // we need to just save the bTagWeight since we only want to apply it
    // for the nominal shape
    //if (cursyst->name == "nominal") {
    //    *f["weight"] = *f["weight"] * *f["bTagWeight"];
    //}

    // Split WJets and ZJets samples by jet parton flavor
    *in["bMCFlavorSum"] = 0;
    *in["bMCFlavorSumSelected"] = 0;
    *in["bGenJetSum"] = 0;
    *in["bGenJetBSum"] = 0;
    
    *in["sampleIndex_sel"] = 
    *in["sampleIndex_GenJetSum"] = 
    *in["sampleIndex_GenJetSumNB"] = 
    *in["sampleIndex"];
    
    if (cursample->doJetFlavorSplit) {
        *in["sampleIndex"] = *in["sampleIndex"]*100;
        *in["sampleIndex_sel"] = *in["sampleIndex_sel"]*100;
        *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"]*100;
        *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"]*100;
       
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd1"]]) == 5)  *in["bMCFlavorSumSelected"]=*in["bMCFlavorSumSelected"]+1;
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd2"]]) == 5)  *in["bMCFlavorSumSelected"]=*in["bMCFlavorSumSelected"]+1;
    
        for(int iJet=0;iJet<*in["nJet"];iJet++){
            if(fabs(in["Jet_mcFlavour"][iJet])==5) *in["bMCFlavorSum"]=*in["bMCFlavorSum"]+1;
        }


        for(int indGJ=0; indGJ<*in["nGenJet"]; indGJ++){
            *in["bGenJetSum"]=*in["bGenJetSum"]+1;
        }

        for(int indGJ=0; indGJ<*in["nGenJet"]; indGJ++){
            *in["bGenJetBSum"]=*in["bGenJetBSum"]+in["GenJet_numBHadrons"][indGJ];
        }
        
        //if(*in["bMCFlavorSum"]==1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 1;
        //}else if(*in["bMCFlavorSum"]>1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 2;
        //}

        if(*in["bMCFlavorSumSelected"]==1){
            *in["sampleIndex_sel"] = *in["sampleIndex_sel"] + 1;
        }else if(*in["bMCFlavorSumSelected"]>1){
            *in["sampleIndex_sel"] = *in["sampleIndex_sel"] + 2;
        }
        
        if(*in["bGenJetBSum"]==1){
            *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"] + 1;
        }else if(*in["bGenJetBSum"]>1){
            *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"] + 2;
        }
        
        if(*in["bGenJetSum"]==1){
            *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"] + 1;
        }else if(*in["bGenJetSum"]>1){
            *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"] + 2;
        }

        // default now the number of b's counted
        *in["sampleIndex"]=*in["sampleIndex_GenJetSum"];
    } 

    // Split by boost category
    if(*in["isWenu"]) {
        if(*f["V_pt"] >= 100 && *f["V_pt"] < 150) *in["eventClass"] += 2;
        else if(*f["V_pt"] >= 150) *in["eventClass"] += 1;
        else *in["eventClass"] += 3;
    }
    else if(*in["isWmunu"]) {
        if(*f["V_pt"] >= 100 && *f["V_pt"] < 130) *in["eventClass"] += 3;
        else if(*f["V_pt"] >= 130 && *f["V_pt"] < 180) *in["eventClass"] += 2;
        else if(*f["V_pt"] >= 180) *in["eventClass"] += 1;
        else *in["eventClass"] += 4;
    }

    // these variables don't exist anymore
    /*// compare our Higgs jet selection with step 2 
    *d["hJets_pt_0_step2"] = d["hJets_pt"][0];
    *d["hJets_pt_1_step2"] = d["hJets_pt"][1];
    *d["hJets_btagCSV_0_step2"] = d["hJets_btagCSV"][0];
    *d["hJets_btagCSV_1_step2"] = d["hJets_btagCSV"][1];*/

    if(*in["sampleIndex"]!=0){
        if (*f["genWeight"] > 0) {
            *f["weight"] = cursample->intWeight;
        } else {
            *f["weight"] = cursample->intWeight * -1;
        }
    }
    else {
        *f["weight"] = 1.0;
    } 

    //*f["Vtype_f"] = (float) *f["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    *f["selLeptons_pt_0"] = f["selLeptons_pt"][*in["lepInd"]];
    *f["selLeptons_eta_0"] = f["selLeptons_eta"][*in["lepInd"]];
    //*f["selLeptons_phi_0"] = f["selLeptons_phi"][*in["lepInd"]];
    //*in["selLeptons_pdgId_0"] = in["selLeptons_pdgId"][*in["lepInd"]];    
    //*in["selLeptons_eleCutIdCSA14_25ns_v1_0"] = in["selLeptons_eleCutIdCSA14_25ns_v1"][*in["lepInd"]];
    //*in["selLeptons_tightId_0"] = in["selLeptons_tightId"][*in["lepInd"]];
    *f["selLeptons_relIso03_0"] = f["selLeptons_relIso03"][*in["lepInd"]];  
  
    // electron ID variables
    //*f["selLeptons_eleSieie_0"] = f["selLeptons_eleSieie"][*in["lepInd"]];
    //*f["selLeptons_eleDEta_0"] = f["selLeptons_eleDEta"][*in["lepInd"]];
    //*f["selLeptons_eleDPhi_0"] = f["selLeptons_eleDPhi"][*in["lepInd"]];
    //*f["selLeptons_eleHoE_0"] = f["selLeptons_eleHoE"][*in["lepInd"]];
    //*f["selLeptons_eleMissingHits_0"] = f["selLeptons_eleMissingHits"][*in["lepInd"]];
    //*f["selLeptons_eleChi2_0"] = f["selLeptons_eleChi2"][*in["lepInd"]];
     

    if(*in["sampleIndex"]!=0){
        if (*in["nGenLep"] == 0) {
            // gen lep is originally a tau?
            *f["selLeptons_genLepDR_0"] = -1;
        } else {
            TLorentzVector GenLep, El;
            GenLep.SetPtEtaPhiM(f["GenLep_pt"][0],f["GenLep_eta"][0],f["GenLep_phi"][0],f["GenLep_mass"][0]);
            El.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], f["selLeptons_mass"][*in["lepInd"]]);
            *f["selLeptons_genLepDR_0"] = El.DeltaR(GenLep);
        }       
    }
 

    //*f["naLeptonsPassingCuts"] = *in["nAddLeptons"]; // no need to calculate this twice
    //*f["naJetsPassingCuts"] = *in["nAddJets252p9_puid"];
    //*f["highestCSVaJet"] = 0.0;
    //*f["minDeltaRaJet"] = 9999.0; // FIXME add correct initialization
    //for(int j=0;j<*in["nJet"];j++){
    //    if (j == *in["hJetInd1"] || j == *in["hJetInd2"]) continue;
    //    if(f["Jet_pt_reg"][j]>20 && fabs(f["Jet_eta"][j])<2.5 && f["Jet_btagCSV"][j]>*f["highestCSVaJet"]) *f["highestCSVaJet"] = f["Jet_btagCSV"][j];
    //  }       
    
    // Reconstruct Higgs and W and recalculate variables ourselves
    if(debug>1000) std::cout<<"Making composite candidates"<<std::endl;
    TLorentzVector MET,Lep,W,HJ1,HJ2,Hbb;
    MET.SetPtEtaPhiM(*f["met_pt"], 0., *f["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], f["selLeptons_mass"][*in["lepInd"]]);
    W = MET + Lep;

    HJ1.SetPtEtaPhiM(f["Jet_pt_reg"][*in["hJetInd1"]], f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]]);
    HJ2.SetPtEtaPhiM(f["Jet_pt_reg"][*in["hJetInd2"]], f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]]);
    Hbb = HJ1 + HJ2;
    
    // We already calculate these in Analyze()
    //*d["H_mass"] = Hbb.M();
    //*d["H_pt"] = Hbb.Pt();
    //*d["V_pt"] = W.Pt();
    //*d["HVdPhi"] = Hbb.DeltaPhi(W);
    
    // Set variables used by the BDT
    //*f["H_mass_f"] = (float) *f["H_mass"];
    //*f["H_pt_f"] = (float) *f["H_pt"];
    //*f["V_pt_f"] = (float) *f["V_pt"];
    *f["hJets_btagCSV_0"] = (float) f["Jet_btagCSV"][*in["hJetInd1"]];
    *f["hJets_btagCSV_1"] = (float) f["Jet_btagCSV"][*in["hJetInd2"]]; 
    //*f["HVdPhi_f"] = (float) *f["HVdPhi"];
    //*f["H_dEta"] = fabs(f["Jet_eta"][*in["hJetInd1"]] - f["Jet_eta"][*in["hJetInd2"]]);
  
    //*f["hJets_mt_0"] = HJ1.Mt();
    //*f["hJets_mt_1"] = HJ2.Mt();
    //*f["H_dR"] = (float) HJ1.DeltaR(HJ2);   
    //*f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    *f["hJets_pt_0"] = (float) f["Jet_pt_reg"][*in["hJetInd1"]];
    *f["hJets_pt_1"] = (float) f["Jet_pt_reg"][*in["hJetInd2"]];
    //*f["met_sumEt_f"] = (float) *f["met_sumEt"]; // is this the right variable??
    *f["nAddJet_f"] = (float) *in["nAddJets252p9_puid"];
    *f["nAddLep_f"] = (float) *in["nAddLeptons"];
    *f["isWenu_f"] = (float) *in["isWenu"];
    *f["isWmunu_f"] = (float) *in["isWmunu"];

    if (BDTisSet) {
        if(debug>5000) {std::cout<<"Evaluating BDT..."<<std::endl; }
        for (unsigned int i=0; i<bdtInfos.size(); i++) {

            BDTInfo tmpBDT = bdtInfos[i];
            if(debug>5000) {
                PrintBDTInfoValues(tmpBDT);
                std::cout<<"BDT evaluates to: "<<tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod)<<std::endl;
            }
            std::string bdtname(tmpBDT.bdtname);
            if (cursyst->name != "nominal") {
                bdtname.append("_");
                bdtname.append(cursyst->name);
            }
            *f[bdtname] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod);
        }
    }   

    if(jet1EnergyRegressionIsSet && jet2EnergyRegressionIsSet) {
        *f["hJets_pt_0"] = float(f["Jet_pt"][*in["hJetInd1"]]);
        *f["hJets_corr_0"] = f["Jet_corr"][*in["hJetInd1"]];
        //*f["hJets_rawPt_0"] = float(d["hJets_rawPt"][0]);
        *f["hJets_eta_0"] = float(f["Jet_eta"][*in["hJetInd1"]]);
        *f["hJets_mt_0"] = HJ1.Mt();
        *f["hJets_leadTrackPt_0"] = float(f["Jet_leadTrackPt"][*in["hJetInd1"]]);
        *f["hJets_leptonPtRel_0"] = float(f["Jet_leptonPtRel"][*in["hJetInd1"]]);
        *f["hJets_leptonPt_0"] = float(f["Jet_leptonPt"][*in["hJetInd1"]]);
        *f["hJets_leptonDeltaR_0"] = float(f["Jet_leptonDeltaR"][*in["hJetInd1"]]);
        //*f["hJets_chEmEF_0"] = float(f["Jet_chEmEF"][*in["hJetInd1"]]);
        //*f["hJets_chHEF_0"] = float(f["Jet_chHEF"][*in["hJetInd1"]]);
        *f["hJets_neHEF_0"] = float(f["Jet_neHEF"][*in["hJetInd1"]]);
        *f["hJets_neEmEF_0"] = float(f["Jet_neEmEF"][*in["hJetInd1"]]);
        *f["hJets_chMult_0"] = float(f["Jet_chMult"][*in["hJetInd1"]]);
        //*f["hJets_puId_0"] = float(d["hJets_puId"][0]);
        *f["hJets_vtxMass_0"] = float(f["Jet_vtxMass"][*in["hJetInd1"]]);
        *f["hJets_vtxPt_0"] = float(f["Jet_vtxPt"][*in["hJetInd1"]]);
        //*f["hJets_vtx3DVal_0"] = float(f["Jet_vtx3DVal"][*in["hJetInd1"]]);
        *f["hJets_vtx3dL_0"] = f["Jet_vtx3DVal"][*in["hJetInd1"]];
        *f["hJets_vtxNtracks_0"] = float(f["Jet_vtxNtracks"][*in["hJetInd1"]]);
        *f["hJets_vtx3deL_0"] = f["Jet_vtx3DSig"][*in["hJetInd1"]];

        *f["hJets_pt_1"] = float(f["Jet_pt"][*in["hJetInd2"]]);
        *f["hJets_corr_1"] = f["Jet_corr"][*in["hJetInd2"]];
        //*f["hJets_rawPt_1"] = float(d["hJets_rawPt"][1]);
        *f["hJets_eta_1"] = float(f["Jet_eta"][*in["hJetInd2"]]);
        *f["hJets_mt_1"] = HJ2.Mt();
        *f["hJets_leadTrackPt_1"] = float(f["Jet_leadTrackPt"][*in["hJetInd2"]]);
        *f["hJets_leptonPtRel_1"] = float(f["Jet_leptonPtRel"][*in["hJetInd2"]]);
        *f["hJets_leptonPt_1"] = float(f["Jet_leptonPt"][*in["hJetInd2"]]);
        *f["hJets_leptonDeltaR_1"] = float(f["Jet_leptonDeltaR"][*in["hJetInd2"]]);
        //*f["hJets_chEmEF_1"] = float(f["Jet_chEmEF"][*in["hJetInd2"]]);
        //*f["hJets_chHEF_1"] = float(f["Jet_chHEF"][*in["hJetInd2"]]);
        *f["hJets_neHEF_1"] = float(f["Jet_neHEF"][*in["hJetInd2"]]);
        *f["hJets_neEmEF_1"] = float(f["Jet_neEmEF"][*in["hJetInd2"]]);
        *f["hJets_chMult_1"] = float(f["Jet_chMult"][*in["hJetInd2"]]);
        //*f["hJets_puId_1"] = float(f["hJets_puId"][1]);
        *f["hJets_vtxMass_1"] = float(f["Jet_vtxMass"][*in["hJetInd2"]]);
        //*f["hJets_vtx3DVal_1"] = float(f["Jet_vtx3DVal"][*in["hJetInd2"]]); 
        *f["hJets_vtx3dL_1"] = f["Jet_vtx3DVal"][*in["hJetInd2"]];
        *f["hJets_vtxPt_1"] = float(f["Jet_vtxPt"][1]);
        *f["hJets_vtxNtracks_1"] = float(f["Jet_vtxNtracks"][*in["hJetInd2"]]);
        *f["hJets_vtx3deL_1"] = f["Jet_vtx3DSig"][*in["hJetInd2"]];
        
        if(debug>10000) {
            std::cout<<"Evaluating the Jet Energy Regression..."<<std::endl;
            PrintBDTInfoValues(jet1EnergyRegression);
            PrintBDTInfoValues(jet2EnergyRegression);
        }
        double r1Pt = jet1EnergyRegression.reader->EvaluateRegression(jet1EnergyRegression.bdtmethod)[0];
        double r2Pt = jet2EnergyRegression.reader->EvaluateRegression(jet2EnergyRegression.bdtmethod)[0];

        *f["Jet1_regWeight"] = r1Pt/(*f["hJets_pt_0"]);
        *f["Jet2_regWeight"] = r2Pt/(*f["hJets_pt_1"]);

        TLorentzVector hJ1_reg = TLorentzVector();
        TLorentzVector hJ2_reg = TLorentzVector();
 
        *f["Jet1_pt_reg"] = r1Pt;
        *f["Jet2_pt_reg"] = r2Pt;
        hJ1_reg.SetPtEtaPhiM(r1Pt, f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]]);
        hJ2_reg.SetPtEtaPhiM(r2Pt, f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]]);
       
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


    // FIXME nominal must be last
    if(cursyst->name=="nominal"){
        ofile->cd();
        outputTree->Fill();
    }
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

bool VHbbAnalysis::ElectronSelection(){
    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"Running Wenu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<f["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*d[\"met_pt\"] "<<*f["met_pt"]<<std::endl;
        std::cout<<"*f[selLeptons_eleSieie_0] = "<<*f["selLeptons_eleSieie_0"]<<std::endl;
        std::cout<<"*f[hJets_pt_1] = "<<*f["hJets_pt_1"]<<std::endl;
    }
    
    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts
    *in["elInd"] = -1;
    float elMaxPt = 0; // max pt of the electrons we select
    for (int i =0; i<*in["nselLeptons"]; i++) {
        if(fabs(in["selLeptons_pdgId"][i])==11 
            && f["selLeptons_pt"][i]      > *f["eptcut"] 
            && fabs(f["selLeptons_eta"][i]) < *f["eletacut"]
            && f["selLeptons_relIso03"][i]< *f["erelisocut"]
            && in["selLeptons_eleCutIdSpring15_25ns_v1"][i] >= *f["elidcut"]
            ){
            if (f["selLeptons_pt"][i] > elMaxPt) {
                elMaxPt = f["selLeptons_pt"][i];
                *in["elInd"] = i;
            }
        }
    }
    if (*in["elInd"] == -1) selectEvent = false;

    return selectEvent;
}

bool VHbbAnalysis::MuonSelection(){
    
    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"Running Wmunu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso04\"] "<<f["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*d[\"met_pt\"] "<<*f["met_pt"]<<std::endl;
    }
    
    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts

    *in["muInd"] = -1;
    float muMaxPt = 0; // max pt of the muons we select 
    for (int i=0; i<*in["nselLeptons"]; i++) {
        if(fabs(in["selLeptons_pdgId"][i])==13 
            && f["selLeptons_pt"][i]      > *f["muptcut"]
            && fabs(f["selLeptons_eta"][i]) < *f["muetacut"]
            && f["selLeptons_relIso03"][i]< *f["murelisocut"]
            && in["selLeptons_tightId"][i] >= *f["muidcut"]
            ){
            if (f["selLeptons_pt"][i] > muMaxPt) {
                muMaxPt = f["selLeptons_pt"][i];
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
            && f["Jet_pt_reg"][i]>*f["j1ptCut"]
            && f["Jet_btagCSV"][i]>*f["j1ptCSV"]&&fabs(f["Jet_eta"][i])<=*f["j1etaCut"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(f["Jet_pt_reg"][pair.first]<f["Jet_pt"][i]){
                pair.first = i;
            }
        }
        /*if( pair.first == -1) {
            pair.first = i;
        }
        else if(d["Jet_pt_reg"][pair.first]<d["Jet_pt"][i]){
            pair.first = i;
        }*/
    }
   
    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] == 1
            && f["Jet_pt_reg"][i]>*f["j2ptCut"]
            && f["Jet_btagCSV"][i]>*f["j2ptCSV"]&&fabs(f["Jet_eta"][i])<*f["j2etaCut"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(f["Jet_pt_reg"][pair.second]<f["Jet_pt"][i]){
                pair.second = i;
            }
        }
        /*if (pair.second == -1) {
            pair.second = i;
        }
        else if(d["Jet_pt_reg"][pair.second]<d["Jet_pt"][i]){
            pair.second = i;
        }*/
    }

    return pair;
}


std::pair<int,int> VHbbAnalysis::HighestCSVBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        if(in["Jet_puId"][i] == 1
            && f["Jet_pt_reg"][i]>*f["j1ptCut"]
            &&fabs(f["Jet_eta"][i])<=*f["j1etaCut"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(f["Jet_btagCSV"][pair.first]<f["Jet_btagCSV"][i]){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] == 1
            && f["Jet_pt_reg"][i]>*f["j2ptCut"]
            &&fabs(f["Jet_eta"][i])<*f["j2etaCut"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(f["Jet_btagCSV"][pair.second]<f["Jet_btagCSV"][i]){
                pair.second = i;
            }
        }
    }

    return pair;
}

std::pair<int,int> VHbbAnalysis::HighestPtJJBJets(){
    std::pair<int,int> pair(-1,-1);

    // dont implement the btag CSV cuts until we've already selected the highest pT(jj) jets
    double maxPtJJ = 0.;
    for (int i=0; i<*in["nJet"]; i++) {
        if(in["Jet_puId"][i] == 1
            && f["Jet_pt_reg"][i]>*f["j1ptCut"]
            && fabs(f["Jet_eta"][i])<*f["j1etaCut"]) {
            TLorentzVector jet1;
            jet1.SetPtEtaPhiM(f["Jet_pt_reg"][i],f["Jet_eta"][i],f["Jet_phi"][i],f["Jet_mass"][i]);
            for (int j=0; j<*in["nJet"]; j++) {
                if (i == j) continue;
                if(in["Jet_puId"][j] == 1
                    && f["Jet_pt_reg"][j]>*f["j2ptCut"]
                    && fabs(f["Jet_eta"][j])<*f["j2etaCut"]) {
                    TLorentzVector jet2;
                    jet2.SetPtEtaPhiM(f["Jet_pt_reg"][j],f["Jet_eta"][j],f["Jet_phi"][j],f["Jet_mass"][j]);
                    TLorentzVector jj = jet1 + jet2;
                    double ptJJ = jj.Pt();
                    if (ptJJ >= maxPtJJ) {
                        if (j == pair.first && pair.second == i) {
                            // we've already picked this pair in the other order, order by CSV
                            // quick sanity check, make sure it's the same maxPtJJ
                            if (ptJJ != maxPtJJ) {
                                std::cout<<"Picked both orderings of highest pT(jj) jets, but the orderings don't have the same pT(jj)!!"<<std::endl;
                                std::cout<<"ptJJ = "<<ptJJ<<std::endl;
                                std::cout<<"maxPtJJ = "<<maxPtJJ<<std::endl;
                            }
                            if (f["Jet_btagCSV"][j] > f["Jet_btagCSV"][i]) {
                            //if (f["Jet_pt_reg"][j] > f["Jet_pt"][i]) {
                                pair.first = j;
                                pair.second = i; 
                            }
                            else {
                                pair.first = i;
                                pair.second = j;
                            }
                        }
                        else {
                            pair.first = i;
                            pair.second = j;
                            maxPtJJ = ptJJ;
                        }
                    }
                }
            } 
        }      
    }
    // important to cut on CSV here to kill TTbar
    if (f["Jet_btagCSV"][pair.first] < *f["j1ptCSV"]) pair.first = -1;
    if (f["Jet_btagCSV"][pair.second] < *f["j2ptCSV"]) pair.second = -1;
    return pair;    
}

double VHbbAnalysis::GetRecoTopMass(TLorentzVector Obj, bool isJet, int useMET, bool regPT) {

    // Try to reconstruct the top in ttbar with leptonic W decay
    // if isJet is true, construct top as given jet + closest lepton
    // if isJet is false, construct top as given lepton + closest lepton
    //
    // if useMET is 1, construct the top using the jet + lepton + met and take the transverse mass
    // if useMET is 2, construct the top using the jet + lepton + met, calculate met pZ by assuming mW
    double Top_mass = -999; // return value

    double minDR = 999;
    int ObjClosestIndex = -1; // index of closest lepton if isJet, closet jet otherwise

    TLorentzVector Obj2; // closest lepton if isJet, closest jet otherwise
    if (isJet) {
        // look in aleptons too FIXME
        // find the closest lepton to the given jet
        for (int i=0; i<*in["nselLeptons"]; i++) {
            TLorentzVector l;
            l.SetPtEtaPhiM(f["selLeptons_pt"][i], f["selLeptons_eta"][i], f["selLeptons_phi"][i], f["selLeptons_mass"][i] );
            double d1 = l.DeltaR(Obj);
            if (d1 <= minDR) {
                minDR = d1;
                ObjClosestIndex = i;
             }
        }
        if (ObjClosestIndex!=-1) {
            Obj2.SetPtEtaPhiM(f["selLeptons_pt"][ObjClosestIndex], f["selLeptons_eta"][ObjClosestIndex], f["selLeptons_phi"][ObjClosestIndex], f["selLeptons_mass"][ObjClosestIndex]);
        }
        else return -999;
    }
    else {
        // find closest jet to the given lepton
        float thisPT=0;
        for (int i=0; i<*in["nJet"]; i++) {
            if(regPT){
                thisPT=f["Jet_pt_reg"][i];
            } else {
                thisPT=f["Jet_pt"][i];
            }
            if (thisPT< 30 || f["Jet_btagCSV"][i] < 0.5) continue; // only consider jets with some minimal preselection
            TLorentzVector j;
            j.SetPtEtaPhiM(thisPT, f["Jet_eta"][i], f["Jet_phi"][i], f["Jet_mass"][i] );
            double d1 = j.DeltaR(Obj);
            if (d1 <= minDR) {
                minDR = d1;
                ObjClosestIndex = i;
             }
        }
        if (ObjClosestIndex!=-1) {
            if(regPT){
                thisPT=f["Jet_pt_reg"][ObjClosestIndex];
            } else {
                thisPT=f["Jet_pt"][ObjClosestIndex];
            }
            Obj2.SetPtEtaPhiM(thisPT, f["Jet_eta"][ObjClosestIndex], f["Jet_phi"][ObjClosestIndex], f["Jet_mass"][ObjClosestIndex]);
        }
        else return -999;
    }

    if (useMET==1) {
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*f["met_pt"],0.,*f["met_phi"],0.);
        TLorentzVector Obj_transverse, Obj2_transverse; // can only consider transverse (x-y plane) 4-vector components if using MET
        Obj_transverse.SetPxPyPzE(Obj.Px(),Obj.Py(),0,TMath::Sqrt(TMath::Power(Obj.M(),2) + TMath::Power(Obj.Pt(),2)));
        Obj2_transverse.SetPxPyPzE(Obj2.Px(),Obj2.Py(),0,TMath::Sqrt(TMath::Power(Obj2.M(),2) + TMath::Power(Obj2.Pt(),2)));
        TLorentzVector Top_transverse = Obj_transverse + Obj2_transverse + MET;
        /*if (Top_transverse.M() != Top_transverse.Mt()) {
            std::cout<<"Top_transverse.M() = "<<Top_transverse.M()<<std::endl; 
            std::cout<<"Top_transverse.Mt() = "<<Top_transverse.Mt()<<std::endl;
            std::cout<<"Top_transverse.Px() = "<<Top_transverse.Px()<<std::endl;
            std::cout<<"Top_transverse.Py() = "<<Top_transverse.Py()<<std::endl;
            std::cout<<"Top_transverse.Pz() = "<<Top_transverse.Pz()<<std::endl;
            std::cout<<"Top_transverse.E() = "<<Top_transverse.E()<<std::endl;
        }*/
        return Top_transverse.M(); // particle physics two-particle Mt = sqrt(Et**2 - pt**2)
    }
    else if (useMET==2) {
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*f["met_pt"],0.,*f["met_phi"],0.);
        TLorentzVector lep, jet;
        if (isJet) {
            lep = Obj2; 
            jet = Obj;
        }
        else {
            lep = Obj;
            jet = Obj2;
        }  
        TLorentzVector W_trans = lep + MET; // W w/ transverse neutrino info
        TLorentzVector W; // full W four-vector, to be constructed assuming PDG W mass
        // solve analytically for neutrino pZ
        double mW = 80.376;
        double alpha = TMath::Power(mW,2) + TMath::Power(W_trans.Pt(),2) + TMath::Power(lep.Pz(),2) - TMath::Power(lep.E(),2) - TMath::Power(MET.Pt(), 2);
        double a = 4*TMath::Power(lep.E(),2) - 4*TMath::Power(lep.Pz(),2);
        double b = -4 * alpha * lep.Pz();
        double c = 4*TMath::Power(lep.E(),2)*TMath::Power(MET.Pt(),2) - TMath::Power(alpha,2);
        double pzm_p = (-1*b + TMath::Sqrt(TMath::Power(b,2) - 4*a*c) ) / (2 * a);
        double pzm_m = (-1*b - TMath::Sqrt(TMath::Power(b,2) - 4*a*c) ) / (2 * a);
        //MET.SetPz(pzm_p);
        W.SetPxPyPzE(W_trans.Px(), W_trans.Py(), lep.Pz() + pzm_p, TMath::Sqrt(TMath::Power(lep.E(),2) + TMath::Power(MET.Pt(),2) + TMath::Power(pzm_p,2) ) );  
        std::cout<<"met_pt = "<<MET.Pt()<<std::endl;
        std::cout<<"E_lep = "<<lep.E()<<std::endl;
        std::cout<<"lep_pz = "<<lep.Pz()<<std::endl;;
        std::cout<<"lep_pt = "<<lep.Pt()<<std::endl;;
        std::cout<<"pt_W = "<<W_trans.Pt()<<std::endl;
        std::cout<<"alpha = "<<alpha<<std::endl;
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"c = "<<c<<std::endl;
        std::cout<<"pzm_p = "<<pzm_p<<std::endl;
        std::cout<<"pzm_m = "<<pzm_m<<std::endl;
        std::cout<<"W.M() = "<<W.M()<<std::endl;
        std::cout<<"W.Pt() = "<<W.Pt()<<std::endl;
        std::cout<<"W.Pz() = "<<W.Pz()<<std::endl;
        std::cout<<"W.E() = "<<W.E()<<std::endl;
        TLorentzVector Top = W + jet;
        return Top.M(); 
            
    } 

    TLorentzVector Top = Obj + Obj2;
    return Top.M();
}


float VHbbAnalysis::ReWeightMC(int nPU){

double data2[52]={
4.49873379421785E-005 ,
0.0002734329  ,
0.0002829032  ,
0.0003713701  ,
0.0005505258  ,
0.0008813604  ,
0.0026826432  ,
0.0166379415  ,
0.0618837545  ,
0.1203749801  ,
0.1715793027  ,
0.1959297551  ,
0.1788795228  ,
0.1285242708  ,
0.0723801445  ,
0.0322818448  ,
0.0116512156  ,
0.0035181198  ,
0.0009365339  ,
0.0002378211  ,
6.3484922942407E-005  ,
1.9640165964275E-005  ,
7.51884233725539E-006 ,
3.47269001094576E-006 ,
1.75552059572915E-006 ,
8.94169705326947E-007 ,
4.40644428362066E-007 ,
0.000000207 ,
9.23073857687134E-008 ,
0.000000039 ,
1.56317481250602E-008 ,
5.93482507474814E-009 ,
2.1353501155199E-009  ,
7.28092380422255E-010 ,
2.35266751336713E-010 ,
7.20423931836899E-011 ,
2.09059068601469E-011 ,
5.74917354326211E-012 ,
1.49829687744647E-012 ,
3.70036607069327E-013 ,
8.66064573134044E-014 ,
1.92095530821743E-014 ,
4.03783143396351E-015 ,
8.04352330952033E-016 ,
1.51850475885846E-016 ,
2.71662388084133E-017 ,
4.60882697565656E-018 ,
7.37768210770165E-019 ,
1.15691685161167E-019 ,
1.71283468612116E-020 ,
0 ,
0 ,
};

double mc2[52]={4.8551E-07,
		1.74806E-06,
		3.30868E-06,
		1.62972E-05,
		4.95667E-05,
		0.000606966,
		0.003307249,
		0.010340741,
		0.022852296,
		0.041948781,
		0.058609363,
		0.067475755,
		0.072817826,
		0.075931405,
		0.076782504,
		0.076202319,
		0.074502547,
		0.072355135,
		0.069642102,
		0.064920999,
		0.05725576,
		0.047289348,
		0.036528446,
		0.026376131,
		0.017806872,
		0.011249422,
		0.006643385,
		0.003662904,
		0.001899681,
		0.00095614,
		0.00050028,
		0.000297353,
		0.000208717,
		0.000165856,
		0.000139974,
		0.000120481,
		0.000103826,
		8.88868E-05,
		7.53323E-05,
		6.30863E-05,
		5.21356E-05,
		4.24754E-05,
		3.40876E-05,
		2.69282E-05,
		2.09267E-05,
		1.5989E-05,
		4.8551E-06,
		2.42755E-06,
		4.8551E-07,
		2.42755E-07,
		1.21378E-07,
		4.8551E-08};

  if(nPU<0) return 1;
  if(nPU>51) return 1;
  return data2[nPU]/mc2[nPU];
}

// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float VHbbAnalysis::ptWeightQCD(int nGenVbosons, float lheHT, int GenVbosons_pdgId){
  float SF = 1.;
  if (lheHT>100 && nGenVbosons==1){
    if (GenVbosons_pdgId == 23){ // Z
    SF =   ((lheHT>100 && lheHT<200)*1.588 * ( 280.35 / (409.860000) ) + (lheHT>200 && lheHT<400)*1.438 * ( 77.67 / ( 110.880000 )) + (lheHT>400 && lheHT<600)*1.494 * (10.73 / (13.189 )) + (lheHT>600)*1.139 * ( 4.116 / (4.524300) ));
    }
    if (abs(GenVbosons_pdgId) == 24){
      SF =   ((lheHT>100 && lheHT<200)*1.588 * ( 1345 / (1.23 *  1.29e3) ) + (lheHT>200 && lheHT<400)*1.438 * ( 359.7 / ( 1.23 *  3.86e2)) + (lheHT>400 && lheHT<600)*1.494 * (48.91 / (1.23 * 47.9 )) + (lheHT>600)*1.139 * ( 18.77 / (1.23 * 19.9) ));
    }
  }
  return SF>0?SF:0;
}

// weights correction for EWK NLO correction
// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float VHbbAnalysis::ptWeightEWK(int nGenVbosons,float GenVbosons_pt,int VtypeSim,int GenVbosons_pdgId){
  float SF = 1.;
  if (nGenVbosons ==1)
    {
      if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5)
    {
      if (GenVbosons_pdgId == 23)
        {
          //for Z options
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
        }
    }
      else if (VtypeSim == 2 || VtypeSim == 3)
    {
      //for W options
      if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24)
        {
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041+7.93714*(TMath::Power((GenVbosons_pt+877.978),-0.213831));
        }
    }
    }
  return SF>0?SF:0;
}
