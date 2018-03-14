//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "EWKAnalysis.h"
#include "HelperClasses/EquationSolver.h"
#include "TRandom3.h"

// initialize parameters
EWKAnalysis::EWKAnalysis(){
    if(debug>10) std::cout<<"Constructing EWKAnalysis"<<std::endl;
}

// remove any 'new'
EWKAnalysis::~EWKAnalysis(){
    if(debug>10) std::cout<<"Deconstructing EWKAnalysis"<<std::endl;
}

//probably want to do bdtSetup here
void EWKAnalysis::InitAnalysis(){
    //SetupBDT();a//FIXME need to update BDT
    return;
}

//get a few entries and then preselect
//if sel, then analyzeevent
//default to false in the future
bool EWKAnalysis::Preselection(){
    //return true; // for the moment don't impose any preselection
    bool sel=true;
    //if (*in["evt"]!=378187 && *in["evt"]!=378183 && *in["evt"]!=378179 && *in["evt"]!=378181) { return false; }
    //std::cout<<*in["evt"]<<std::endl;
    //std::cout<<"found it!"<<std::endl;
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
    //if (cursample->sampleNum == 22 && *f["lheHT"] > 100) sel=false; 
    // stitch ZJets inclusive sample to HT-binned samples
    if (cursample->sampleNum == 23 && *f["lheHT"] > 100) return false; 
    
    // Impose trigger requirements
    if (*in["HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v"]!=1 && *in["HLT_BIT_HLT_IsoMu24_v"]!=1 && *in["HLT_BIT_HLT_IsoTkMu24_v"]!=1) return false;
    //if (*in["HLT_BIT_HLT_Ele27_WPTight_Gsf_v"]!=1 && *in["HLT_BIT_HLT_IsoMu24_v"]!=1 && *in["HLT_BIT_HLT_IsoTkMu24_v"]!=1) sel=false;
    
    *f["met_pt_JECUp_ratio"] = *f["met_shifted_JetEnUp_pt"] / *f["met_pt"] ;
    *f["met_pt_JECDown_ratio"] = *f["met_shifted_JetEnDown_pt"] / *f["met_pt"] ;

    if (*f["met_pt"] < *f["metcut"]) return false;

    if (cursample->sampleNum==0) {
        if (*f["json"]!=1) return false;
    }

    // Heppy jet corrections for JER/JEC are full correction, it's easier to just use the
    // correction on top of the nominal
    for (int i=0; i<*in["nJet"]; i++) {
        //if (int(*f["doReg"]) == 0) {
        //    // don't apply the regression. The easiest way to do this is to just reset the value of
        //    // the regressed jet pt to the nominal jet pt, since we build everything later from the 
        //    // jet pt's.
        //    f["Jet_pt_reg"][i] = f["Jet_pt"][i];
       // }
        //else if (int(*f["reReg"]) != 0) {
        //    //std::cout<<"Jet_pt_reg was "<<f["Jet_pt_reg"][i]<<std::endl;
        //    //std::cout<<"regression evaluates to: "<<evaluateRegression(i)<<std::endl;
        //    f["Jet_pt_reg"][i] = evaluateRegression(i);
        //    //std::cout<<"Now jet_pt_reg is "<<f["Jet_pt_reg"][i]<<std::endl;
       // }
        //float JERScale = *f["JERScale"]; // apply JER smearing x times the nominal smearing amount
        ////float JERScale = 1.0; // apply JER smearing x times the nominal smearing amount
        //if (JERScale != 1.0) {
        //    smearJets(JERScale);
       // }

        //std::cout<<*in["nJet"]<<": "<<f["Jet_pt_reg_corrJECUp"][i]<<" / "<<f["Jet_pt"][i]<<std::endl;
        f["Jet_corr_JERUp_ratio"][i] = f["Jet_corr_JERUp"][i] / f["Jet_corr_JER"][i] ; 
        f["Jet_corr_JERDown_ratio"][i] = f["Jet_corr_JERDown"][i] / f["Jet_corr_JER"][i] ; 
        f["Jet_corr_JECUp_ratio"][i] = f["Jet_corr_JECUp"][i] / f["Jet_corr"][i] ; 
        f["Jet_corr_JECDown_ratio"][i] = f["Jet_corr_JECDown"][i] / f["Jet_corr"][i] ; 
        //f["Jet_pt_reg_corrJECUp_ratio"][i] = f["Jet_pt_reg_corrJECUp"][i] / f["Jet_pt"][i] ; 
        //f["Jet_pt_reg_corrJECDown_ratio"][i] = f["Jet_pt_reg_corrJECDown"][i] / f["Jet_pt"][i] ; 
        //f["Jet_pt_reg_corrJERUp_ratio"][i] = f["Jet_pt_reg_corrJERUp"][i] / f["Jet_pt"][i] ; 
        //f["Jet_pt_reg_corrJERDown_ratio"][i] = f["Jet_pt_reg_corrJERDown"][i] / f["Jet_pt"][i] ; 
    }
  
    //SetupFactorizedJECs(cursyst->name);

    //// for some reason sometimes in Heppy trees met < 0, which leads to annoying error messages
    //if (*f["met_pt"] <= 0.) {
    //    *f["met_pt"] == 0.0001;
    //}    

    // Preselect for two jets and one lepton which pass some minimum pt threshold
    int nPreselJets = 0;
    int indJet1 = -1;
    int indJet2 = -1;
    for (int i=0; i < *in["nJet"]; i++) {
        if (fabs(f["Jet_eta"][i]) < *f["jetEtaCut"] && in["Jet_puId"][i] > 0 && in["Jet_id"][i] > 0 ) {
            if (f["Jet_pt"][i] > *f["Jet1PtPresel"] && nPreselJets==0) {
                nPreselJets++;
                indJet1 = i;
            }
            else if (nPreselJets==1 && f["Jet_pt"][i] > *f["Jet2PtPresel"]) {
                nPreselJets++;
                indJet2 = i;
                break;
            }
        }
    }
    int nPreselLep = 0;
    for (int i=0; i < *in["nselLeptons"]; i++) {
        if (f["selLeptons_pt"][i] > *f["LepPtPresel"] && fabs(f["selLeptons_eta"][i]) < 2.5 ) nPreselLep++;
    }
    //for (int i=0; i < *in["naLeptons"]; i++) {
    //    if (f["aLeptons_pt"][i] > *f["LepPtPresel"]) nPreselLep++;
    //}
   
 
    if ((*f["Vtype"]!=2 && *f["Vtype"]!=3) || nPreselJets < 2 || nPreselLep < 1) return false;

    TLorentzVector J1,J2,Hjj;
    J1.SetPtEtaPhiM(f["Jet_pt"][indJet1], f["Jet_eta"][indJet1], f["Jet_phi"][indJet1], f["Jet_mass"][indJet1] * f["Jet_pt"][indJet1] );
    J2.SetPtEtaPhiM(f["Jet_pt"][indJet2], f["Jet_eta"][indJet2], f["Jet_phi"][indJet2], f["Jet_mass"][indJet2] * f["Jet_pt"][indJet2] );
    Hjj = J1 + J2;
    float Mjj = Hjj.M();
    if (Mjj < *f["mjjcut"]) return false;

    return sel;
}

bool EWKAnalysis::Analyze(){
    *in["sampleIndex"] = cursample->sampleNum;
    bool sel=true;

    //for (int i=0; i<*in["nJet"]; i++) {
        //if (int(*f["doReg"]) == 0) {
        //    // don't apply the regression. The easiest way to do this is to just reset the value of
        //    // the regressed jet pt to the nominal jet pt, since we build everything later from the 
            // jet pt's.
        //    f["Jet_pt"][i] = f["Jet_pt"][i];
        //}
        //else if (int(*f["reReg"]) != 0) {
        //    f["Jet_pt"][i] = evaluateRegression(i);
        //}
        //float JERScale = *f["JERScale"]; // apply JER smearing x times the nominal smearing amount
        //if (JERScale != 1.0) {
        //    smearJets(JERScale);
        //}
    //}

    if(debug>1000) {
         std::cout<<"Imposing trigger and json requirements"<<std::endl;
    }
    // Impose trigger requirements
    if (*in["HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v"]!=1 && *in["HLT_BIT_HLT_IsoMu24_v"]!=1 && *in["HLT_BIT_HLT_IsoTkMu24_v"]!=1) return false;
    //if (*in["HLT_BIT_HLT_Ele27_WPTight_Gsf_v"]!=1 && *in["HLT_BIT_HLT_IsoMu24_v"]!=1 && *in["HLT_BIT_HLT_IsoTkMu24_v"]!=1) sel=false;


    if (*in["sampleIndex"]==0) {
        if (*f["json"]!=1) return false;
    }
 
    *in["JetInd1"] = -1;
    *in["JetInd2"] = -1;
    bool leadJetIsFound = false;
    for (int i=0; i < *in["nJet"]; i++) {
        if (fabs(f["Jet_eta"][i]) < *f["jetEtaCut"] && in["Jet_puId"][i] > 0 && in["Jet_id"][i] > 0 ) {
            if (f["Jet_pt"][i] > *f["j1ptCut"]) {
                if (leadJetIsFound) {
                    *in["JetInd2"] = i;
                     break;
                }   
                else {
                    *in["JetInd1"] = i;
                    leadJetIsFound = true;
                }
            }
            else if (leadJetIsFound && f["Jet_pt"][i] > *f["j2ptCut"]) {
                *in["JetInd2"] = i;
                break;
            }
        }
    }

    if (*in["JetInd1"]==-1 || *in["JetInd2"] == -1) return false;

    *in["JetInd3"] = -1;
    *in["JetInd4"] = -1;

    for (int i=0; i < *in["nJet"]; i++) {
        if (i == *in["JetInd1"] || i == *in["JetInd2"]) {
            continue;
        }
        if (*in["JetInd3"] == -1) {
            *in["JetInd3"] = i;
        }
        else if (*in["JetInd4"] == -1) {
            *in["JetInd4"] = i;
        }
    } 
    
 

    if(debug>1000) {
        std::cout<<"nJet = "<<*in["nJet"]<<std::endl;
        std::cout<<"JetInd1 = "<<*in["JetInd1"]<<std::endl;
        std::cout<<"JetInd2 = "<<*in["JetInd2"]<<std::endl;
        std::cout<<"found two jets with pt "
            <<f["Jet_pt"][*in["JetInd1"]]<<" "
            <<f["Jet_pt"][*in["JetInd2"]]<<" "
            <<std::endl;
    }
    
    TLorentzVector J1,J2,Hjj;
    J1.SetPtEtaPhiM(f["Jet_pt"][*in["JetInd1"]], f["Jet_eta"][*in["JetInd1"]], f["Jet_phi"][*in["JetInd1"]], f["Jet_mass"][*in["JetInd1"]] * f["Jet_pt"][*in["JetInd1"]] );
    J2.SetPtEtaPhiM(f["Jet_pt"][*in["JetInd2"]], f["Jet_eta"][*in["JetInd2"]], f["Jet_phi"][*in["JetInd2"]], f["Jet_mass"][*in["JetInd2"]] * f["Jet_pt"][*in["JetInd2"]] );
    Hjj = J1 + J2;
    *f["Mjj"] = Hjj.M();
    if (debug > 10000) {
        std::cout<<"Cutting on M(jj): "<<*f["Mjj"]<<std::endl;
    }
    if (*f["Mjj"] < *f["mjjcut"]) return false;
    
    TLorentzVector J1_noreg,J2_noreg,Hjj_noreg;
    J1_noreg.SetPtEtaPhiM(f["Jet_pt"][*in["JetInd1"]], f["Jet_eta"][*in["JetInd1"]], f["Jet_phi"][*in["JetInd1"]], f["Jet_mass"][*in["JetInd1"]]);
    J2_noreg.SetPtEtaPhiM(f["Jet_pt"][*in["JetInd2"]], f["Jet_eta"][*in["JetInd2"]], f["Jet_phi"][*in["JetInd2"]], f["Jet_mass"][*in["JetInd2"]]);
    Hjj_noreg = J1_noreg + J2_noreg;
    *f["Mjj_noreg"] = Hjj_noreg.M();

    *f["PtJJ"] = Hjj.Pt();
   
    if (*f["met_pt"] < *f["metcut"]) return false;

    //check if event passes any event class
    *in["lepInd"] = -1;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    if(MuonSelection() && (cursample->lepFlav==-1 || cursample->lepFlav==0)) {
        *in["isWmunu"] = 1;
        *in["lepInd"] = *in["muInd"];
    }
    if(ElectronSelection() && (cursample->lepFlav==-1 || cursample->lepFlav==1)) {
        *in["isWenu"] = 1;
        if (!(*in["isWmunu"])) {
            *in["lepInd"] = *in["elInd"];
        }
    }
    if (*in["isWmunu"] == 0 && *in["isWenu"] == 0) return false;
    
    // calculate and apply Rochester corrections for muon momentum
    float rochesterCorr = 1.0;
    if (*in["isWmunu"] == 1) {
        //RoccoR  *rc = new RoccoR("/uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/EWKAnalysis/aux/roccor.2016.v3/rcdata.2016.v3/"); 
        if (*in["sampleIndex"] == 0) {
            // data event
            rochesterCorr = rc->kScaleDT(in["selLeptons_charge"][*in["lepInd"]], f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], 0, 0);
        }
        else {
            TRandom3 *rand = new TRandom3();
            double u1 = rand->Rndm();
            double u2 = rand->Rndm();
            rochesterCorr = rc->kScaleAndSmearMC(in["selLeptons_charge"][*in["lepInd"]], f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]],  in["selLeptons_trackerLayers"][*in["lepInd"]], u1, u2, 0, 0);
        }
        f["selLeptons_pt"][*in["lepInd"]] = f["selLeptons_pt"][*in["lepInd"]] * rochesterCorr;
    }
    *f["rochesterCorr"] = rochesterCorr;
     
    *f["selLeptons_pt_0"] = f["selLeptons_pt"][*in["lepInd"]];
    *f["selLeptons_eta_0"] = f["selLeptons_eta"][*in["lepInd"]];
    
    if (*in["lepInd"] == -1) {
        // not Wenu or Wmunu, use preselected lepton
        *in["lepInd"] = 0;
    }

    if(debug>1000) {
        std::cout<<"cutting on dphi(lep, met)"<<std::endl;
        std::cout<<*in["lepInd"]<<" "<<f["selLeptons_phi"][*in["lepInd"]]<<" "<<f["selLeptons_eta"][*in["lepInd"]]<<" "<<f["selLeptons_mass"][*in["lepInd"]]<<std::endl;
    }
    *f["lepMetDPhi"]=fabs(EvalDeltaPhi(f["selLeptons_phi"][*in["lepInd"]],*f["met_phi"]));
    if (*in["isWmunu"] && *f["lepMetDPhi"] > *f["muMetDPhiCut"])  return false;
    else if (*in["isWenu"] && *f["lepMetDPhi"] > *f["elMetDPhiCut"]) return false;

    TLorentzVector W,Lep, MET, GenBJ1;
    // Reconstruct W
    if(debug>1000) {
        std::cout<<"met "<<*f["met_pt"]<<" "<<*f["met_phi"]<<std::endl;
    }
    MET.SetPtEtaPhiM(*f["met_pt"], 0., *f["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], f["selLeptons_mass"][*in["lepInd"]]); 
    W = MET + Lep; 
    *f["V_pt"] = W.Pt(); // uncomment this line if we want to recalculate W.pt ourselves
    *f["V_eta"] = W.Eta();   
 
    double cosPhi12 = ( Lep.Px()*MET.Px() + Lep.Py()*MET.Py() ) / ( Lep.Pt()*MET.Pt() ); // cos of the angle between the lepton and the missing energy
    *f["V_mt"] = TMath::Sqrt( 2*Lep.Pt()*MET.Pt() * (1 - cosPhi12) );
    if (*f["V_mt"] < *f["wmtcut"]) return false;
    *f["jjWPtBalance"] = *f["PtJJ"] / *f["V_pt"];

    *f["Lep_J1_dPhi"] = Lep.DeltaPhi(J1);
    *f["Lep_J2_dPhi"] = Lep.DeltaPhi(J2);
    *f["Lep_J1_dEta"] = fabs( Lep.Eta() - J1.Eta() );
    *f["Lep_J2_dEta"] = fabs( Lep.Eta() - J1.Eta() );
    *f["Lep_J1_dR"] = Lep.DeltaR(J1);
    *f["Lep_J2_dR"] = Lep.DeltaR(J2);
    *f["J1_J2_dPhi"] = J1.DeltaPhi(J2);
    *f["J1_J2_dEta"] = fabs( J1.Eta() - J2.Eta());
    *f["J1_J2_dR"] = J1.DeltaR(J2);

    float minDPhi = 100; // minimum deltaphi between missing energy and jet (only selected jets)
    float minDPhi2 = 100; // minimum deltaphi between missing energy and jet (all jets)
    float minDPhi3 = 100; // minimum deltaphi between missing energy and two pT-leading jets
    for (int i=0; i<*in["nJet"]; i++) {
        TLorentzVector jet_tmp;
        jet_tmp.SetPtEtaPhiM(f["Jet_pt"][i],f["Jet_eta"][i],f["Jet_phi"][i],f["Jet_mass"][i]);
        float dphi = fabs(jet_tmp.DeltaPhi(MET));
        if (dphi < minDPhi2) {
            minDPhi2 = dphi;
        }
        if (dphi < minDPhi && (i == *in["JetInd1"] || i == *in["JetInd2"])) {
                minDPhi = dphi;
        }
        if (dphi < minDPhi3 && i < 2) {
                minDPhi3 = dphi;
        }
        
    }
    *f["minJMetDPhi"] = minDPhi;
    *f["minJMetDPhiAll"] = minDPhi2;
    *f["minJMetDPhiPtL"] = minDPhi3;
    


    // Now we can calculate whatever we want (transverse) with W and H four-vectors
    *f["JJVdPhi"] = Hjj.DeltaPhi(W);
    *f["JJVdEta"] = fabs(Hjj.Eta() - W.Eta());
    *f["yW_noMET"] = f["selLeptons_eta"][*in["lepInd"]] - (f["Jet_eta"][*in["JetInd1"]] + f["Jet_eta"][*in["JetInd2"]] )/2.;
    *f["zW_noMET"] = fabs(*f["yW_noMET"]) / (fabs(f["Jet_eta"][*in["JetInd1"]] - f["Jet_eta"][*in["JetInd2"]]) );
    *f["yW"] = *f["V_eta"] - (f["Jet_eta"][*in["JetInd1"]] + f["Jet_eta"][*in["JetInd2"]] )/2.;
    *f["zW"] = fabs(*f["yW"]) / (fabs(f["Jet_eta"][*in["JetInd1"]] - f["Jet_eta"][*in["JetInd2"]]) );
    TLorentzVector neutrino = getNu4Momentum(Lep, MET);
    TLorentzVector W_4MET = neutrino + Lep;
    *f["JJVdEta_4MET"] = fabs(Hjj.Eta() - W_4MET.Eta());
    *f["V_eta_4MET"] = W_4MET.Eta();
    *f["yW_4MET"] = *f["V_eta_4MET"] - (f["Jet_eta"][*in["JetInd1"]] + f["Jet_eta"][*in["JetInd2"]] )/2.;
    *f["zW_4MET"] = fabs(*f["yW_4MET"]) / (fabs(f["Jet_eta"][*in["JetInd1"]] - f["Jet_eta"][*in["JetInd2"]]) );

    *f["JJEtaBal"] = ( fabs(f["Jet_eta"][*in["JetInd1"]] + f["Jet_eta"][*in["JetInd2"]]) ) / (fabs(f["Jet_eta"][*in["JetInd1"]] - f["Jet_eta"][*in["JetInd2"]]) ) ;
    *f["JJPtBal"] = *f["PtJJ"] / (f["Jet_pt"][*in["JetInd1"]] + f["Jet_pt"][*in["JetInd2"]]);
    TLorentzVector JJL = Hjj + Lep; 
    *f["JJL_pt"] = JJL.Pt();
    TLorentzVector JJW = Hjj + W;
     *f["JJW_pt"] = JJW.Pt();
    //*f["RPt"] = *f["JJL_pt"] / (f["Jet_pt"][*in["JetInd1"]] + f["Jet_pt"][*in["JetInd2"]] + f["selLeptons_pt"][*in["lepInd"]]);
    *f["RPt"] = *f["JJW_pt"] / (f["Jet_pt"][*in["JetInd1"]] + f["Jet_pt"][*in["JetInd2"]] + *f["V_pt"]);
 
    //std::cout<<cursyst->name.c_str()<<": Hbb.M() = "<<Hbb.M()<<", Jet_pt_reg[jetInd1] = "<<f["Jet_pt"][jetInd1]<<", Jet_pt_reg[jetInd2] = "<<f["Jet_pt"][jetInd2]<<", "<<(f["Jet_corr_JECUp"][jetInd1] / f["Jet_corr"][jetInd1])<<", "<<(f["Jet_corr_JECDown"][jetInd1] / f["Jet_corr"][jetInd1])<<std::endl;
    if (cursyst->name != "nominal") {
        *f[Form("Mjj_%s",cursyst->name.c_str())] = Hjj.M();
        *f[Form("PtJJ_%s",cursyst->name.c_str())] = Hjj.Pt();
    }


    if(debug>1000) std::cout<<"counting additional leptons"<<std::endl;
    // count the number of additional leptons and jets, then cut on this number
    int nAddLep = 0;       
    // 15 to 30 by 1 GeV, 1.5 to 3 w/ 0.1 in eta 
    //std::vector<int> ptCuts = {15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
    //std::vector<double> etaCuts = {1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5};
    std::vector<int> ptCuts = {25};
    std::vector<double> etaCuts = {2.5};

    for (int i=0; i < (int)ptCuts.size(); i++) {
        for (int j=0; j < (int)etaCuts.size(); j++) {
            float maxPt = 0; // max. pT of additional jets
            int nAddJet_tmp = 0;
            std::string eta_cut = Form("%.1f", etaCuts[j]); // convert, say, 2.5 to '2p5' 
            std::replace(eta_cut.begin(), eta_cut.end(), '.', 'p');
            for (int k=0; k < *in["nJet"]; k++) {
                if (k == *in["JetInd1"] || k == *in["JetInd2"]) continue;
                if (f["Jet_pt"][k]>ptCuts[i] && fabs(f["Jet_eta"][k])<etaCuts[j] && in["Jet_puId"][k]>0 && in["Jet_id"][k]>0) {
                    nAddJet_tmp++;
                    if (f["Jet_pt"][k] > maxPt) {
                        maxPt = f["Jet_pt"][k];
                        *f[Form("AddJets%i%s_leadJet_pt",ptCuts[i],eta_cut.c_str())] = maxPt;
                        *f[Form("AddJets%i%s_leadJet_eta",ptCuts[i],eta_cut.c_str())] = f["Jet_eta"][k];
                        *f[Form("AddJets%i%s_leadJet_phi",ptCuts[i],eta_cut.c_str())] = f["Jet_phi"][k];
                        *f[Form("AddJets%i%s_leadJet_btagCSV",ptCuts[i],eta_cut.c_str())] = f["Jet_btagCMVA"][k];
                    } 
                }
            }
            //std::cout<<Form("nAddJets%i%s_puid",ptCuts[i],eta_cut.c_str())<<std::endl;
            *in[Form("nAddJets%i%s",ptCuts[i],eta_cut.c_str())] = nAddJet_tmp;
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
    if (nAddLep>= *f["nAddLeptonsCut"]) return false;

    if(*in["nAddJets252p5"] >= *f["nAddJetsCut"]) return false;
    
    //if(fabs(*f["JJVdPhi"]) < *f["JJVDPhiCut"]) reurn false;

    return sel;
}

void EWKAnalysis::FinishEvent(){
   
    //if (bool(*f["doCutFlow"])) {
    //    ofile->cd();
    //    outputTree->Fill();
    //    return;
    //} 
    // General use variables

    if (cursample->sampleNum == 49 || cursample->sampleNum == 491) {
        *f["nProcEvents"] = cursample->CountWeighted->GetBinContent(1);
    }
    else { 
        *f["nProcEvents"] = cursample->processedEvents;
    }
    if (*in["sampleIndex"] != 0) {
        for (int i=0; i<*in["nLHE_weights_scale"];i++) {
            TH1F *CountWeightedLHEWeightScale = cursample->CountWeightedLHEWeightScale;
            f["LHE_weights_scale_normwgt"][i] = *f["nProcEvents"] / CountWeightedLHEWeightScale->GetBinContent(CountWeightedLHEWeightScale->FindBin(i));
        }
        for (int i=0; i<*in["nLHE_weights_pdf"];i++) {
            TH1F *CountWeightedLHEWeightPdf = cursample->CountWeightedLHEWeightPdf;
            if (CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i)) != 0) {
                f["LHE_weights_pdf_normwgt"][i] = *f["nProcEvents"] / CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i));
            }
            else {
                f["LHE_weights_pdf_normwgt"][i] = 1.0;
            }
            //cout<<f["LHE_weights_pdf_normwgt"][i]<<" = "<<*f["nProcEvents"]<<" / "<<CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i))<<endl;
        }
    }
    //cout<<f["LHE_weights_pdf_normwgt"][5]<<endl;

    if(*in["sampleIndex"]!=0){
        if (*f["genWeight"] > 0) {
            *f["weight"] = cursample->intWeight;
        } else {
            *f["weight"] = cursample->intWeight * -1;
        }
        if (cursample->sampleNum == 49 || cursample->sampleNum == 491) {
            // special prescription for WJets_BGenFilter sample
            *f["weight"] = *f["weight"]*fabs(*f["genWeight"]);
        } 
    }
    else {
        *f["weight"] = 1.0;
    }
     
    *f["weight_ptQCD"] = 1.0;
    *f["weight_ptEWK"] = 1.0;
    if(*in["sampleIndex"]!=0){
            *f["weight_PU"] = *f["puWeight"];
            *f["weight_PUUp"] = *f["puWeightUp"] / *f["puWeight"];
            *f["weight_PUDown"] = *f["puWeightDown"] / *f["puWeight"];
        if (*in["nGenTop"]==0 && *in["nGenVbosons"]>0) {
            // only apply to Z/W+jet samples
            *f["weight_ptQCD"]=ptWeightQCD(*in["nGenVbosons"], *f["lheHT"], in["GenVbosons_pdgId"][0]);
            *f["weight_ptEWK"]=ptWeightEWK(*in["nGenVbosons"], f["GenVbosons_pt"][0], *f["VtypeSim"], in["GenVbosons_pdgId"][0]);
        }
    } else {
        *f["weight_PU"]=1;
        *f["weight_PUUp"]=1;
        *f["weight_PUUp"]=1;
        *f["puWeight"]=1;
    }


    //*f["Vtype_f"] = (float) *f["Vtype"];
    //*f["absDeltaPullAngle"] = fabs(*f["deltaPullAngle"]);
    //*f["selLeptons_pt_0"] = f["selLeptons_pt"][*in["lepInd"]];
    //*f["selLeptons_eta_0"] = f["selLeptons_eta"][*in["lepInd"]];
    //*f["selLeptons_phi_0"] = f["selLeptons_phi"][*in["lepInd"]];
    //*in["selLeptons_pdgId_0"] = in["selLeptons_pdgId"][*in["lepInd"]];    
    //*in["selLeptons_eleCutIdCSA14_25ns_v1_0"] = in["selLeptons_eleCutIdCSA14_25ns_v1"][*in["lepInd"]];
    //*in["selLeptons_tightId_0"] = in["selLeptons_tightId"][*in["lepInd"]];
    if (*in["isWmunu"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso04"][*in["lepInd"]];  
    }
    else if (in["isWenu"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso03"][*in["lepInd"]];  
    }
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
 

    
    // Reconstruct Higgs and W and recalculate variables ourselves
    if(debug>1000) std::cout<<"Making composite candidates"<<std::endl;
    TLorentzVector MET,Lep,W;
    MET.SetPtEtaPhiM(*f["met_pt"], 0., *f["met_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd"]], f["selLeptons_eta"][*in["lepInd"]], f["selLeptons_phi"][*in["lepInd"]], f["selLeptons_mass"][*in["lepInd"]]);
    W = MET + Lep;

    
    // We already calculate these in Analyze()
    //*d["H_mass"] = Hbb.M();
    //*d["H_pt"] = Hbb.Pt();
    //*d["V_pt"] = W.Pt();
    //*d["HVdPhi"] = Hbb.DeltaPhi(W);
    
    // Set variables used by the BDT
    //*f["V_pt_f"] = (float) *f["V_pt"];
    //*f["HVdPhi_f"] = (float) *f["HVdPhi"];
  
    //*f["Jet1_mt"] = HJ1.Mt();
    //*f["Jet2_mt"] = HJ2.Mt();
    //*f["H_dR"] = (float) HJ1.DeltaR(HJ2);   
    //*f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    //*f["met_sumEt_f"] = (float) *f["met_sumEt"]; // is this the right variable??
    *f["nAddJet_f"] = (float) *in["nAddJets252p5"];
    *f["nAddLep_f"] = (float) *in["nAddLeptons"];
    *f["isWenu_f"] = (float) *in["isWenu"];
    *f["isWmunu_f"] = (float) *in["isWmunu"];
    *f["softActivityVH_njets5_f"] = (float) *in["softActivityVH_njets5"];
    *f["Jet1_pt"] = f["Jet_pt"][*in["JetInd1"]];
    *f["Jet2_pt"] = f["Jet_pt"][*in["JetInd2"]];
    *f["Jet1_eta"] = f["Jet_eta"][*in["JetInd1"]];
    *f["Jet2_eta"] = f["Jet_eta"][*in["JetInd2"]];
    *f["Jet1_qgl"] = f["Jet_qgl"][*in["JetInd1"]];
    *f["Jet2_qgl"] = f["Jet_qgl"][*in["JetInd2"]];
    //*f["nPVs_f"] = (float) *in["nPVs"];

    if (f.count("bdt_1lep")>0) {
        if(debug>5000) {
            std::cout<<"Evaluating BDT..."<<std::endl;
            PrintBDTInfoValues(bdtInfos["bdt_1lep"]);
            std::cout<<"BDT evaluates to: "<<EvaluateMVA(bdtInfos["bdt_1lep"])<<std::endl;
        }
        std::string bdtname(bdtInfos["bdt_1lep"]->bdtname);
        if (cursyst->name != "nominal") {
            bdtname.append("_");
            bdtname.append(cursyst->name);
        }
        if (m("doVPtReweighting") > 0) {
            float V_pt_nom = m("V_pt");
            *f["V_pt"] = m("V_pt_VPtCorrUp");
            *f[bdtname+"_VPtCorrUp"] = EvaluateMVA(bdtInfos["bdt_1lep"]);
            *f["V_pt"] = m("V_pt_VPtCorrDown");
            *f[bdtname+"_VPtCorrDown"] = EvaluateMVA(bdtInfos["bdt_1lep"]);
            *f["V_pt"] = V_pt_nom;
        }
        *f[bdtname] = EvaluateMVA(bdtInfos["bdt_1lep"]);
    }
    for(std::map<std::string,BDTInfo*>::iterator iterBDT=bdtInfos.begin();
          iterBDT!=bdtInfos.end(); iterBDT++){
     
        if(debug>5000) {
            PrintBDTInfoValues(iterBDT->second);
            std::cout<<"BDT evaluates to: "<<EvaluateMVA(iterBDT->second)<<std::endl;
        }
        std::string bdtname(iterBDT->second->bdtname);
        if (cursyst->name != "nominal") {
            bdtname.append("_");
            bdtname.append(cursyst->name);
        }
        *f[bdtname] = EvaluateMVA(iterBDT->second);
    }

    if (*in["sampleIndex"]!=0) {
        if (*in["isWmunu"] == 1) {
                //*f["Lep_SF"] = ( (20.1/36.4) * f["SF_MuIDTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIDTightGH"][*in["lepInd"]]) * ( (20.1/36.4) * f["SF_MuIsoTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIsoTightGH"][*in["lepInd"]] ) ; 
                *f["Lep_SF"] = ( (20.1/36.4) * f["SF_MuIDTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIDTightGH"][*in["lepInd"]]) * ( (20.1/36.4) * f["SF_MuIsoTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIsoTightGH"][*in["lepInd"]] ) *  ( (20.1/36.4) * f["SF_MuTriggerBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuTriggerGH"][*in["lepInd"]]) * ( (20.1/36.4) * f["SF_MuTrackerBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuTrackerGH"][*in["lepInd"]]); 
                *f["Lep_SFUp"] = ( (20.1/36.4) * (f["SF_MuIDTightBCDEF"][*in["lepInd"]] + f["SF_MuIDTightBCDEF_err"][*in["lepInd"]]) + (16.3/36.4) * (f["SF_MuIDTightGH"][*in["lepInd"]] + f["SF_MuIDTightGH_err"][*in["lepInd"]]) ) * ( (20.1/36.4) * (f["SF_MuIsoTightBCDEF"][*in["lepInd"]] + f["SF_MuIsoTightBCDEF_err"][*in["lepInd"]] ) + (16.3/36.4) * (f["SF_MuIsoTightGH"][*in["lepInd"]] + f["SF_MuIsoTightGH_err"][*in["lepInd"]]) ) *  ( (20.1/36.4) * (f["SF_MuTriggerBCDEF"][*in["lepInd"]] + f["SF_MuTriggerBCDEF_err"][*in["lepInd"]] )+ (16.3/36.4) * (f["SF_MuTriggerGH"][*in["lepInd"]]) + f["SF_MuTriggerGH_err"][*in["lepInd"]]) * ( (20.1/36.4) * (f["SF_MuTrackerBCDEF"][*in["lepInd"]] + f["SF_MuTrackerBCDEF_err"][*in["lepInd"]])+ (16.3/36.4) * (f["SF_MuTrackerGH"][*in["lepInd"]] + f["SF_MuTrackerGH_err"][*in["lepInd"]] ));
                *f["Lep_SFUp"] = *f["Lep_SFUp"] / *f["Lep_SF"];
                *f["Lep_SFDown"] = ( (20.1/36.4) * (f["SF_MuIDTightBCDEF"][*in["lepInd"]] - f["SF_MuIDTightBCDEF_err"][*in["lepInd"]]) + (16.3/36.4) * (f["SF_MuIDTightGH"][*in["lepInd"]] - f["SF_MuIDTightGH_err"][*in["lepInd"]]) ) * ( (20.1/36.4) * (f["SF_MuIsoTightBCDEF"][*in["lepInd"]] - f["SF_MuIsoTightBCDEF_err"][*in["lepInd"]] ) + (16.3/36.4) * (f["SF_MuIsoTightGH"][*in["lepInd"]] - f["SF_MuIsoTightGH_err"][*in["lepInd"]]) ) *  ( (20.1/36.4) * (f["SF_MuTriggerBCDEF"][*in["lepInd"]] - f["SF_MuTriggerBCDEF_err"][*in["lepInd"]] )+ (16.3/36.4) * (f["SF_MuTriggerGH"][*in["lepInd"]]) - f["SF_MuTriggerGH_err"][*in["lepInd"]]) * ( (20.1/36.4) * (f["SF_MuTrackerBCDEF"][*in["lepInd"]] - f["SF_MuTrackerBCDEF_err"][*in["lepInd"]])+ (16.3/36.4) * (f["SF_MuTrackerGH"][*in["lepInd"]] - f["SF_MuTrackerGH_err"][*in["lepInd"]] )); 
                *f["Lep_SFDown"] = *f["Lep_SFDown"] / *f["Lep_SF"];
         
        }
        else if (*in["isWenu"] == 1) {
      
               *f["Lep_SF"] = f["SF_SingleElTrigger"][*in["lepInd"]] * f["SF_ElIdIso"][*in["lepInd"]] *  f["SF_egammaEffi_tracker"][*in["lepInd"]];
               *f["Lep_SFUp"] = (f["SF_SingleElTrigger"][*in["lepInd"]] + f["SF_SingleElTrigger_err"][*in["lepInd"]] )* (f["SF_ElIdIso"][*in["lepInd"]] + f["SF_ElIdIso_err"][*in["lepInd"]] ) *  (f["SF_egammaEffi_tracker"][*in["lepInd"]] + f["SF_egammaEffi_tracker_err"][*in["lepInd"]] );
                *f["Lep_SFUp"] = *f["Lep_SFUp"] / *f["Lep_SF"];
               *f["Lep_SFDown"] = (f["SF_SingleElTrigger"][*in["lepInd"]] - f["SF_SingleElTrigger_err"][*in["lepInd"]] )* (f["SF_ElIdIso"][*in["lepInd"]] - f["SF_ElIdIso_err"][*in["lepInd"]] ) *  (f["SF_egammaEffi_tracker"][*in["lepInd"]] - f["SF_egammaEffi_tracker_err"][*in["lepInd"]] );
                *f["Lep_SFDown"] = *f["Lep_SFDown"] / *f["Lep_SF"];
                
               // failing SF to correct for data/MC differences in lepton veto yields 
               *f["LepFail_SF"] = 1.0;
               *f["LepFail_SFUp"] = 1.0;
               *f["LepFail_SFDown"] = 1.0;
               if (*in["nselLeptons"] > 1) {
                   for (int i=0; i < *in["nselLeptons"]; i++) {
                       if (i == *in["lepInd"]) continue;
                       *f["LepFail_SF"] = f["SF_EleVeto"][i];
                       *f["LepFail_SFUp"] = f["SF_EleVeto"][i] + f["SF_EleVeto_err"][i];
                       *f["LepFail_SFUp"] = *f["LepFail_SFUp"] / *f["LepFail_SF"];
                       *f["LepFail_SFDown"] = f["SF_EleVeto"][i] - f["SF_EleVeto_err"][i];
                       *f["LepFail_SFDown"] = *f["LepFail_SFDown"] / *f["LepFail_SF"];
                       break; // apply for highest pT lepton besides selected one
                   }
               }

               //*f["Lep_SF"] = f["SF_ElIdMVATrigWP80"][*in["lepInd"]] * f["SF_egammaEffi_tracker"][*in["lepInd"]];
               //*f["Lep_SF"] = 1.0;
           
        }
      
        if (*in["sampleIndex"]!=0) {
            *f["corrQGL"] = getQGLCorr(in["Jet_partonFlavour"][*in["JetInd1"]], f["Jet_eta"][*in["JetInd1"]], f["Jet_qgl"][*in["JetInd1"]]);
            *f["corrQGL"] = *f["corrQGL"] *  getQGLCorr(in["Jet_partonFlavour"][*in["JetInd2"]], f["Jet_eta"][*in["JetInd2"]], f["Jet_qgl"][*in["JetInd2"]]);
            *f["corrQGLUp"] = *f["corrQGL"];
            *f["corrQGLDown"] = 1. / *f["corrQGL"]; 
        }

        *f["weight"] = *f["weight"] * *f["weight_PU"] * *f["weight_ptEWK"] * *f["Lep_SF"] * *f["corrQGL"];
        //*f["weight"] = *f["weight"] * *f["weight_PU"] * *f["weight_ptQCD"] * *f["weight_ptEWK"] * *f["Lep_SF"];

        // Add NLO to LO W+jet re-weighting from Z(ll)H(bb)
        float deta_bb = fabs(f["Jet_eta"][*in["JetInd1"]] - f["Jet_eta"][*in["JetInd2"]]);
        int sampleIndex = *in["sampleIndex"];
        float WJetNLOWeight = 1.0;
        float weight_ptQCD = *f["weight_ptQCD"];
        if (sampleIndex<4100 || sampleIndex>4902) { WJetNLOWeight = 1.0; }
        else if (sampleIndex==4100 || sampleIndex==4200 || sampleIndex==4300 || sampleIndex==4400 || sampleIndex==4500 || sampleIndex==4600 || sampleIndex==4700 || sampleIndex==4800 || sampleIndex==4900 || sampleIndex==48100 || sampleIndex==49100) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 0);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        else if (sampleIndex==4101 || sampleIndex==4201 || sampleIndex==4301 || sampleIndex==4401 || sampleIndex==4501 || sampleIndex==4601 || sampleIndex==4701 || sampleIndex==4801 || sampleIndex==4901 || sampleIndex==48101 || sampleIndex==49101) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 1);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        else if (sampleIndex==4102 || sampleIndex==4202 || sampleIndex==4302 || sampleIndex==4402 || sampleIndex==4502 || sampleIndex==4602 || sampleIndex==4702 || sampleIndex==4802 || sampleIndex==4902 || sampleIndex==48102 || sampleIndex==49102) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 2);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        *f["weight"] = *f["weight"] * WJetNLOWeight;
        *f["WJetNLOWeight"] = WJetNLOWeight;
        
    }
    
    // preserve normalization after applying QGL corrections
    std::map<std::string, float> QGLMap;
    QGLMap["EWKWJets"] = 0.9461;
    QGLMap["EWKWJets_herwig"] = 0.9461;
    QGLMap["IntEWKWJets"] = 8.973620e-01;
    QGLMap["ZJets_0J"] = 0.9693;
    QGLMap["ZJets_1J"] = 0.9610;
    QGLMap["ZJets_2J"] = 0.9524;
    QGLMap["DYToLL_HT100to200"] = 0.9536;
    QGLMap["DYToLL_HT200to400"] = 0.9563;
    QGLMap["DYToLL_HT400to600"] = 0.9530;
    QGLMap["DYToLL_HT600to800"] = 0.9441;
    QGLMap["DYToLL_HT800to1200"] = 0.9362;
    QGLMap["DYToLL_HT1200to2500"] = 0.9227;
    QGLMap["DYToLL_HT2500toInf"] = 0.9063;
    QGLMap["DYToLL_madgraph"] = 0.9628;
    QGLMap["TToLeptons_t_powheg"] = 0.9844;
    QGLMap["TT_powheg"] = 0.9813;
    QGLMap["TToLeptons_s"] = 0.9898;
    QGLMap["TBarToLeptons_t_powheg"] = 0.9821;
    QGLMap["T_tW"] = 0.9746;
    QGLMap["Tbar_tW"] = 0.9750;
    QGLMap["WJets"] = 0.9683;
    QGLMap["WJets_0J"] = 0.9828;
    QGLMap["WJets_1J"] = 0.9674;
    QGLMap["WJets_2J"] = 0.9545;
    QGLMap["WW_fil"] = 0.9600;
    QGLMap["WZ_fil"] = 0.9600;
    QGLMap["ZZ_fil"] = 0.9591;
    QGLMap["WJets-HT100to200"] = 0.9585;
    QGLMap["WJets-HT200to400"] = 0.9586;
    QGLMap["WJets-HT400to600"] = 0.9558;
    QGLMap["WJets-HT600to800"] = 0.9489;
    QGLMap["WJets-HT800to1200"] = 0.9424;
    QGLMap["WJets-HT1200to2500"] = 0.9285;
    QGLMap["WJets-HT2500toInf"] = 0.9116;
    QGLMap["WJets_madgraph"] = 0.9614;
    QGLMap["WJetsMadInc"] = 0.9614;
    QGLMap["QCD_Pt_120to170_EMEnriched"] = 0.994663;
    QGLMap["QCD_Pt_120to170_MuEnrichedPt5"] = 0.994221;
    QGLMap["QCD_Pt_170to300_EMEnriched"] = 0.982772;
    QGLMap["QCD_Pt_170to300_MuEnrichedPt5"] = 1.002450;
    QGLMap["QCD_Pt_30to50_EMEnriched"] = 1.000000;
    QGLMap["QCD_Pt_30to50_MuEnrichedPt5"] = 1.031954;
    QGLMap["QCD_Pt_50to80_EMEnriched"] = 0.980436;
    QGLMap["QCD_Pt_50to80_MuEnrichedPt5"] = 0.949880;
    QGLMap["QCD_Pt_80to120_EMEnriched"] = 1.004608;
    QGLMap["QCD_Pt_80to120_MuEnrichedPt5"] = 1.032725;
    QGLMap["QCD_HT100to200"] = 1.0; 
 
    *f["corrQGL_norm"] = QGLMap[cursample->sampleName];

    // FIXME nominal must be last
    if(cursyst->name=="nominal"){
        ofile->cd();
        outputTree->Fill();
    }
    return;
}

void EWKAnalysis::TermAnalysis(){
    if(debug>10) std::cout<<"START TermAnalysis()"<<std::endl;
    ofile->cd();
    outputTree->Write();
    ofile->Close();
    if(debug>10) std::cout<<"DONE TermAnalysis()"<<std::endl;
    return;
}

bool EWKAnalysis::ElectronSelection(){
    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"Running Wenu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<f["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*f[\"met_pt\"] "<<*f["met_pt"]<<std::endl;
        std::cout<<"*f[selLeptons_eleSieie_0] = "<<*f["selLeptons_eleSieie_0"]<<std::endl;
        std::cout<<"*f[Jet2_pt] = "<<*f["Jet2_pt"]<<std::endl;
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
            //&& in["selLeptons_eleCutIdSpring15_25ns_v1"][i] >= *f["elidcut"]
            //&& in["selLeptons_tightId"][i] >= *f["elidcut"]
            //&& in["selLeptons_eleMVAIdSpring15Trig"][i] >= *f["elidcut"]
            && in["selLeptons_eleMVAIdSppring16GenPurp"][i] >= *f["elidcut"]
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

bool EWKAnalysis::MuonSelection(){
    
    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"Running Wmunu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso04\"] "<<f["selLeptons_relIso04"][0]<<std::endl;
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
            && f["selLeptons_relIso04"][i]< *f["murelisocut"]
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




// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float EWKAnalysis::ptWeightQCD(int nGenVbosons, float lheHT, int GenVbosons_pdgId){
  float SF = 1.;
  if (lheHT>100 && nGenVbosons==1){
    if (GenVbosons_pdgId == 23){ // Z
    SF =   ((lheHT>=100 && lheHT<200)*1.588 * ( 280.35 / (409.860000) ) + (lheHT>=200 && lheHT<400)*1.438 * ( 77.67 / ( 110.880000 )) + (lheHT>=400 && lheHT<600)*1.494 * (10.73 / (13.189 )) + (lheHT>=600)*1.139 * ( 4.116 / (4.524300) ));
    }
    if (abs(GenVbosons_pdgId) == 24){
      SF =   ((lheHT>=100 && lheHT<200)*1.588 * ( 1345 / (1.23 *  1.29e3) ) + (lheHT>=200 && lheHT<400)*1.438 * ( 359.7 / ( 1.23 *  3.86e2)) + (lheHT>=400 && lheHT<600)*1.494 * (48.91 / (1.23 * 47.9 )) + (lheHT>=600)*1.139 * ( 18.77 / (1.23 * 19.9) ));
    }
  }
  return SF>0?SF:0;
}

// weights correction for EWK NLO correction
// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float EWKAnalysis::ptWeightEWK(int nGenVbosons,float GenVbosons_pt,int VtypeSim,int GenVbosons_pdgId){
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


// W-jet NLO to LO re-weighting function from Z(ll)H(bb)
double EWKAnalysis::LOtoNLOWeightBjetSplitEtabb(double etabb, int njets){

    double SF = 1.;
    if(etabb < 5){
        if(njets < 1){
            SF =   0.935422 + 0.0403162*etabb -0.0089026*etabb*etabb +0.0064324*etabb*etabb*etabb -0.000212443*etabb*etabb*etabb*etabb;
        }else if(njets == 1){
            SF =   0.962415 +0.0329463*etabb -0.0414479*etabb*etabb +0.0240993*etabb*etabb*etabb -0.00278271*etabb*etabb*etabb*etabb;
        }else if(njets >=2){
            SF =   (0.721265 -0.105643*etabb -0.0206835*etabb*etabb +0.00558626*etabb*etabb*etabb)*TMath::Exp(0.450244*etabb);
        }
    }

    return SF;

}
void EWKAnalysis::smearJets(float JERScale) {
    for (int i=0; i<*in["nJet"]; i++) {
            float corr_nominal = f["Jet_corr_JER"][i];
            if (corr_nominal == -99) { continue; }
            //std::cout<<"original smearing correction: "<<corr_nominal<<std::endl;
            f["Jet_corr_JER"][i] = 1 + JERScale * (f["Jet_corr_JER"][i]-1);
            f["Jet_corr_JERUp"][i] = 1 + JERScale * (f["Jet_corr_JERUp"][i]-1);
            f["Jet_corr_JERDown"][i] = 1 + JERScale * (f["Jet_corr_JERDown"][i]-1);
            //std::cout<<"smearing jet "<<JERScale<<" times the nominal JER by factor "<<(f["Jet_corr_JER"][i] / corr_nominal)<<std::endl;
            //std::cout<<"Jet_pt = "<<f["Jet_pt"][i]<<std::endl;
            f["Jet_pt"][i] = f["Jet_pt"][i] * ( f["Jet_corr_JER"][i] / corr_nominal);
           // std::cout<<"Re-smeared jet_pt = "<<f["Jet_pt"][i]<<std::endl;
        
            // re-evaluate jet energy regression on re-smeared jets
            //std::cout<<"Jet_pt_reg was: "<<f["Jet_pt"][i]<<std::endl;
            f["Jet_pt_reg_Heppy"][i] = f["Jet_pt"][i];
            f["Jet_pt"][i] = evaluateRegression(i);
            //std::cout<<"now it has been re-evaluated to: "<<f["Jet_pt"][i]<<std::endl;
    }
}
float EWKAnalysis::evaluateRegression(int i) {
            if (f["Jet_pt"][i] == -99) { return -99; }
            *f["Jet1_pt"] = float(f["Jet_pt"][i]);
            *f["Jet1_eta"] = float(f["Jet_eta"][i]);
            TLorentzVector tmp;
            tmp.SetPtEtaPhiM(f["Jet_pt"][i],f["Jet_eta"][i],f["Jet_phi"][i],f["Jet_mass"][i]);
            *f["Jet1_mt"] = tmp.Mt();
            //std::cout<<"4-vector Mt() is "<<tmp.Mt()<<std::endl;
            //float mt = TMath::Sqrt( TMath::Power(tmp.Et(),2) - TMath::Power(tmp.Pt(),2) );
            //std::cout<<"by-hand Mt is "<<mt<<std::endl;
            //*f["Jet1_mt"] = mt;
            *f["Jet1_leadTrackPt"] = float(f["Jet_leadTrackPt"][i]);
            if (f["Jet_leptonPtRel"][i] > 0) {
                *f["Jet1_leptonPtRel"] = float(f["Jet_leptonPtRel"][i]);
            }
            else {
                *f["Jet1_leptonPtRel"] = 0.;
            }
            if (f["Jet_leptonPt"][i] > 0) {
                *f["Jet1_leptonPt"] = float(f["Jet_leptonPt"][i]);
            }
            else {
                *f["Jet1_leptonPt"] =  0.;
            }
            if (f["Jet_leptonDeltaR"][i] > 0) {
                *f["Jet1_leptonDeltaR"] = float(f["Jet_leptonDeltaR"][i]);
            }
            else {
                *f["Jet1_leptonDeltaR"] = 0.;
            }
            *f["Jet1_neHEF"] = float(f["Jet_neHEF"][i]);
            *f["Jet1_neEmEF"] = float(f["Jet_neEmEF"][i]);
            *f["Jet1_vtxMass"] = float(f["Jet_vtxMass"][i]);
            *f["Jet1_vtxPt"] = float(f["Jet_vtxPt"][i]);
            if (f["Jet_vtx3DVal"][i] > 0) {
                *f["Jet1_vtx3dL"] = f["Jet_vtx3DVal"][i];
            }
            else {
                *f["Jet1_vtx3dL"] = 0.;
            } 
            *f["Jet1_vtxNtracks"] = float(f["Jet_vtxNtracks"][i]);
            //*f["hJets_vtx3deL_0"] = f["Jet_vtx3DSig"][i];
            if (f["Jet_vtx3DSig"][i] > 0) {
                *f["Jet1_vtx3deL"] = f["Jet_vtx3DVal"][i] / f["Jet_vtx3DSig"][i]; 
            }
            else {
                *f["Jet1_vtx3deL"] = 0;
            }
            //*f["nPVs_f"] = (float) *in["nPVs"];
            //PrintBDTInfoValues(jet1EnergyRegression);
            return EvaluateRegression(bdtInfos["bdt_bjetreg"]);
}

float EWKAnalysis::getQGLCorr_q(float x) {
    return (-0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581);
}

float EWKAnalysis::getQGLCorr_g(float x) {
    return (-55.7067*x*x*x*x*x*x*x + 113.218*x*x*x*x*x*x -21.1421*x*x*x*x*x -99.927*x*x*x*x + 92.8668*x*x*x -34.3663*x*x + 6.27*x + 0.612992);
}

float EWKAnalysis::getQGLCorr(float partonFlav, float eta, float qgl) {
    float corr = 1.0;
    if (partonFlav != 0 && fabs(eta) < 2.0 && qgl > 0) {
        if (fabs(partonFlav) < 4) {
            corr = getQGLCorr_q(qgl);
        }
        else if (fabs(partonFlav) == 21) {
            corr = getQGLCorr_g(qgl);
        }
    }
    return corr;
}

TLorentzVector EWKAnalysis::getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET)
{
  TLorentzVector Lepton;
  Lepton.SetPxPyPzE(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
  TLorentzVector MET;
  MET.SetPxPyPzE(TMET.Px(), TMET.Py(), 0., TMET.E());

  double  mW = 80.38;

  //std::vector<math::XYZTLorentzVector> result;
  std::vector<TLorentzVector> result;

  //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );
    
  double MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py());
  double mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py();
  double a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
  double a2 = TMath::Power(a,2);
  double b  = (TMath::Power(Lepton.Energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.Energy(),2) - TMath::Power(Lepton.Pz(),2));
  double pz1(0),pz2(0),pznu(0);
  //int nNuSol(0);

  //math::XYZTLorentzVector p4nu_rec;
  TLorentzVector p4nu_rec;
  TLorentzVector p4W_rec;
  TLorentzVector p4b_rec;
  TLorentzVector p4Top_rec;
  TLorentzVector p4lep_rec;

  p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());
  
  //math::XYZTLorentzVector p40_rec(0,0,0,0);

  if(a2-b > 0 ){
    //if(!usePositiveDeltaSolutions_)
    //  {
    //    result.push_back(p40_rec);
    //    return result;
    //  }
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
    //nNuSol = 2;

    //if(usePzPlusSolutions_)pznu = pz1;    
    //if(usePzMinusSolutions_)pznu = pz2;
    //if(usePzAbsValMinimumSolutions_){
      pznu = pz1;
      if(fabs(pz1)>fabs(pz2)) pznu = pz2;
    //}
    double Enu = sqrt(MisET2 + pznu*pznu);

    p4nu_rec.SetPxPyPzE(MET.Px(), MET.Py(), pznu, Enu);

    result.push_back(p4nu_rec);

  }
  else{
    //if(!useNegativeDeltaSolutions_){
    //  result.push_back(p40_rec);
    //  return result;
    //}
    //    double xprime = sqrt(mW;

    double ptlep = Lepton.Pt(),pxlep=Lepton.Px(),pylep=Lepton.Py(),metpx=MET.Px(),metpy=MET.Py();

    double EquationA = 1;
    double EquationB = -3*pylep*mW/(ptlep);
    double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
    double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

    std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

    std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);

    double deltaMin = 14000*14000;
    double zeroValue = -mW*mW/(4*pxlep); 
    double minPx=0;
    double minPy=0;

        //std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl; 
      
    //if(usePxMinusSolutions_){
      for( int i =0; i< (int)solutions.size();++i){
      if(solutions[i]<0 ) continue;
      double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep); 
      double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
      double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 

            //std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl; 

      if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
      minPx=p_x;
      minPy=p_y;}
           //std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    //if(usePxPlusSolutions_){
      for( int i =0; i< (int)solutions2.size();++i){
        if(solutions2[i]<0 ) continue;
        double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep); 
        double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
        double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 
          //std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
        if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
          minPx=p_x;
          minPy=p_y;
        }
                //std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
    double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);

    if(deltaMin==14000*14000) return TLorentzVector(0,0,0,0);
    //if(deltaMin==14000*14000) return result.front();
    //    else std::cout << " test " << std::endl;

    if(delta2ZeroValue < deltaMin){
      deltaMin = delta2ZeroValue;
      minPx=zeroValue;
      minPy=pyZeroValue;}

    //    std::cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
    ///    ////Y part   

    double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
    double a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
    pznu = a_Minimum;

    //if(!useMetForNegativeSolutions_){
      double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
      p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
    //}
    //else{
    //  pznu = a;
    //  double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
    //  p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
    //}
    result.push_back(p4nu_rec);
  }
  return result.front();
}
