//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "VHbbTrigger.h"
#include <TLorentzVector.h>

// initialize parameters
VHbbTrigger::VHbbTrigger(){
    if(debug>10) std::cout<<"Constructing VHbbTrigger"<<std::endl;
}

// remove any 'new'
VHbbTrigger::~VHbbTrigger(){
    if(debug>10) std::cout<<"Deconstructing VHbbTrigger"<<std::endl;
}

//probably want to do bdtSetup here
void VHbbTrigger::InitAnalysis(){
    //SetupBDT();a//FIXME need to update BDT
    return;
}

//get a few entries and then preselect
//if sel, then analyzeevent
//default to false in the future
bool VHbbTrigger::Preselection(){
    bool sel=false;

    if(cursample->sampleNum<0){
        if(debug>100000){
            std::string str[] = {"b1Pt_gen","genJ1PtCut","b1Phi_gen", "b2Pt_gen","genJ2PtCut","b2Phi_gen" };
            std::vector<std::string> names(str, str + ( sizeof ( str ) /  sizeof ( std::string ) ) );  
            for(int ind=0; ind< (int)names.size(); ind++){
                std::cout<<names[ind]<<" "<<*f[names[ind]]<<std::endl;
            }
        }

        if( *f["b1Pt_gen"]          <*f["genJ1PtCut"]
            ||fabs(*f["b1Phi_gen"]) >*f["genJ1PhiCut"]
            ||*f["b2Pt_gen"]        <*f["genJ2PtCut"]
            ||fabs(*f["b2Phi_gen"]) >*f["genJ2PhiCut"]
            ||*f["genElePt"]        <*f["genElPtCut"]){
            return sel;
        }
    }

    int npresell1Jets=0;
    for(int il1Jet=0; il1Jet<*in["nl1jets"]; il1Jet++){
        if(f["l1jetPt"][il1Jet] > *f["jetPtPreCut"]) npresell1Jets++;
    }
    if(npresell1Jets>=(int)*f["preselNJets"]) sel=true;

    return sel;
}

bool VHbbTrigger::Analyze(){
    bool sel=false;

    *b["L1_EG20"]       =   PassEGL1(20.0,false,false);
    *b["L1_HTT175"] = (*f["l1htTot"]>175);  
    *b["L1_ETM70"] = (*f["l1metPt"]>70); 
 
    if( *b["L1_EG20"] || *b["L1_HTT175"] || *b["L1_ETM70"] ) sel=true;
    else if ((*f["genElePt"]>20&&abs(*f["genEleEta"])<2.5&&*f["b2Pt_gen"]>25&&*f["b1Pt_gen"]>25&&abs(*f["b2Eta_gen"])<2.5&&abs(*f["b1Eta_gen"])<2.5&& *f["genHPt"]>100)){
        sel=true;
    } else return sel;

    *b["L1_EG40"]       =   PassEGL1(40.0,false,false);
    *b["L1_EG40er"]     =   PassEGL1(40.0,false,true);
    *b["L1_EG40iso"]    =   PassEGL1(40.0,true, false);
    *b["L1_EG35"]       =   PassEGL1(35.0,false,false);
    *b["L1_EG35er"]     =   PassEGL1(35.0,false,true);
    *b["L1_EG35iso"]    =   PassEGL1(35.0,true, false);
    *b["L1_EG35isoer"]  =   PassEGL1(35.0,true, true);
    *b["L1_EG30"]       =   PassEGL1(30.0,false,false);
    *b["L1_EG30iso"]    =   PassEGL1(30.0,true, false);
    *b["L1_EG30isoer"]  =   PassEGL1(30.0,true, true);
    *b["L1_EG25"]       =   PassEGL1(25.0,false,false);
    *b["L1_EG25er"]     =   PassEGL1(25.0,false,true);
    *b["L1_EG25iso"]    =   PassEGL1(25.0,true, false);
    *b["L1_EG25isoer"]  =   PassEGL1(25.0,true, true);
    *b["L1_EG22isoer"]  =   PassEGL1(22.0,true, true);
    *b["L1_EG20iso"]    =   PassEGL1(20.0,true, false);
    *b["L1_EG20isoer"]  =   PassEGL1(20.0,true, true);
    

    *b["L1EG20iso_L1Jet25Eta35"]   =   PassEGPlusJetL1(20., true, false, 25, 3.5);
    *b["L1EG25_L1Jet40Eta25"]      =   PassEGPlusJetL1(25., false,false, 40, 2.5);
    *b["L1EG25er_L1Jet40Eta25"]    =   PassEGPlusJetL1(25., false,true,  40, 2.5);
    *b["L1EG25iso_L1Jet40Eta25"]   =   PassEGPlusJetL1(25., true, false, 40, 2.5);
    *b["L1EG25iso_L1Jet50Eta25"]   =   PassEGPlusJetL1(25., true, false, 50, 2.5);
    *b["L1EG20iso_L1Jet50Eta25"]   =   PassEGPlusJetL1(20., true, false, 50, 2.5);
    *b["L1EG20iso_L1Jet40Eta25"]   =   PassEGPlusJetL1(20., true, false, 40, 2.5);
    *b["L1EG25eriso_L1Jet40Eta25"] =   PassEGPlusJetL1(25., true, true , 40, 2.5);
    *b["L1EG30iso_L1Jet30Eta30"]   =   PassEGPlusJetL1(30., true, false, 30, 3.0);
    *b["L1EG30iso_L1Jet30Eta25"]   =   PassEGPlusJetL1(30., true, false, 30, 2.5);


    return sel;
}

void VHbbTrigger::FinishEvent(){
    *in["sampleIndex"] = cursample->sampleNum;
    *f["weight"] = cursample->intWeight;
    
    if(*f["doGenMatching"]==1) GenMatchElectrons();
    if(*f["recomputeWPs"]==1) RecomputeWPs();

  
    std::vector<std::string> L1seeds;
    L1seeds.clear();
    L1seeds.push_back("L1_EG30isoer");
    *b["HLT_Ele32_WP75"] = PassEleHLT(32., "WP75", L1seeds, true, in["HLT_Ele32_WP75_ElInd"]);
    *b["HLT_Ele32_WP75_NoER"] = PassEleHLT(32., "WP75", L1seeds, false);
    *b["HLT_Ele32_WP85"] = PassEleHLT(32., "WP85", L1seeds, true);
    L1seeds.push_back("L1_EG40");
    *b["HLT_Ele32_WP75_ADD40"] = PassEleHLT(32., "WP75", L1seeds, true);
    *b["HLT_Ele32_WP85_ADD40"] = PassEleHLT(32., "WP85", L1seeds, true);
    L1seeds.clear();
    L1seeds.push_back("L1_EG30isoer");
    L1seeds.push_back("L1_EG40");
    *b["HLT_Ele35_WP85_NoER"] = PassEleHLT(35., "WP85", L1seeds, false);
    *b["HLT_Ele35_WP85"] = PassEleHLT(35., "WP85", L1seeds, true);
    *b["HLT_Ele40_WP85"] = PassEleHLT(40., "WP85", L1seeds, true, in["HLT_Ele40_WP85_ElInd"]);
    *b["HLT_Ele45_WP85"] = PassEleHLT(45., "WP85", L1seeds, true);
    *b["HLT_Ele40_WP85_NoER"] = PassEleHLT(40., "WP85", L1seeds, false);
    *b["HLT_Ele45_WP85_NoER"] = PassEleHLT(45., "WP85", L1seeds, false);
    *b["HLT_Ele45_WP75_NoER"] = PassEleHLT(45., "WP75", L1seeds, false);
    *b["HLT_Ele40_WP75_NoER"] = PassEleHLT(40., "WP75", L1seeds, false);
    *b["HLT_Ele50_WP85"] = PassEleHLT(50., "WP85", L1seeds, true, in["HLT_Ele50_WP85_ElInd"]);
    *b["HLT_Ele40_WP75"] = PassEleHLT(40., "WP75", L1seeds, true);
    *b["HLT_Ele45_WP75"] = PassEleHLT(45., "WP75", L1seeds, true);
       
    L1seeds.clear();
    L1seeds.push_back("L1_EG30isoer");
    L1seeds.push_back("L1_EG35er");
    *b["HLT_Ele40_WP85_OLDL1"] = PassEleHLT(40., "WP85", L1seeds, true);
    *b["HLT_Ele40_WP85_OLDL1_NoER"] = PassEleHLT(40., "WP85", L1seeds, false);
    *b["HLT_Ele35_WP85_OLDL1"] = PassEleHLT(35., "WP85", L1seeds, true);
    *b["HLT_Ele35_WP85_OLDL1_NoER"] = PassEleHLT(35., "WP85", L1seeds, false);
    
    L1seeds.clear();
    L1seeds.push_back("L1_EG30isoer");
    *b["HLT_Ele50_WP85_No40"] = PassEleHLT(50., "WP85", L1seeds, true);
        
    L1seeds.clear();
    L1seeds.push_back("L1EG30iso_L1Jet30Eta30");
    //"L1EG20iso_L1Jet25Eta35"
    
    *b["HLT_Ele30_WP85_PFJet30"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 30, 3.0, in["HLT_Ele30_WP85_PFJet30_ElInd"], in["HLT_Ele30_WP85_PFJet30_JetInd"]);
    *b["HLT_Ele30_WP85_PFJet50"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 50, 3.0, in["HLT_Ele30_WP85_PFJet50_ElInd"], in["HLT_Ele30_WP85_PFJet50_JetInd"]);
   
    L1seeds.clear();
    L1seeds.push_back("L1EG25_L1Jet40Eta25");
    *b["HLT_Ele30_WP85_PFJet50_L1EG25_L1Jet40Eta25"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 50, 3.0, in["HLT_Ele30_WP85_PFJet50_L1EG25_L1Jet40Eta25_ElInd"], in["HLT_Ele30_WP85_PFJet50_L1EG25_L1Jet40Eta25_JetInd"]);
    
    L1seeds.clear();
    L1seeds.push_back("L1EG25er_L1Jet40Eta25");
    *b["HLT_Ele30_WP85_PFJet50_L1EG25er_L1Jet40Eta25"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 50, 3.0, in["HLT_Ele30_WP85_PFJet50_L1EG25er_L1Jet40Eta25_ElInd"], in["HLT_Ele30_WP85_PFJet50_L1EG25er_L1Jet40Eta25_JetInd"]);
    
    L1seeds.clear();
    L1seeds.push_back("L1EG25iso_L1Jet40Eta25");
    *b["HLT_Ele30_WP85_PFJet50_L1EG25iso_L1Jet40Eta25"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 50, 3.0, in["HLT_Ele30_WP85_PFJet50_L1EG25iso_L1Jet40Eta25_ElInd"], in["HLT_Ele30_WP85_PFJet50_L1EG25iso_L1Jet40Eta25_JetInd"]);
    
    L1seeds.clear();
    L1seeds.push_back("L1EG25eriso_L1Jet40Eta25");
    *b["HLT_Ele30_WP85_PFJet50_L1EG25eriso_L1Jet40Eta25"] = PassElePlusJetHLT(L1seeds, 30, "WP85", true, 50, 3.0, in["HLT_Ele30_WP85_PFJet50_L1EG25eriso_L1Jet40Eta25_ElInd"], in["HLT_Ele30_WP85_PFJet50_L1EG25eriso_L1Jet40Eta25_JetInd"]);
    
    *in["eventClass"]=0;
 
    L1seeds.clear();
    //L1seeds.push_back("L1_EG20");
    L1seeds.push_back("L1_EG30isoer");
    L1seeds.push_back("L1_HTT175");
    L1seeds.push_back("L1_ETM70");
    int eg20ind=-1;
    bool HLT_Ele20=PassEleHLT(20., "WP85", L1seeds, false, &eg20ind);
    *b["HLT_40CSV0p5_ETA2p4_BESTPT"]=FindBestCSVHLTJet(50, 2.4, 0.5);//, std::string discriminator, int* jetInd)
    *b["HLT_Ele20_WP85_40CSV0p5_ETA2p4_BESTPT"] = *b["HLT_40CSV0p5_ETA2p4_BESTPT"] && HLT_Ele20;
    std::pair<int,int> pair=(eg20ind!=-1)?HighestPtPFJets(30,2.5,f["hltElePhi"][eg20ind]):HighestPtPFJets(30,2.5);

    TLorentzVector dijet(0,0,0,0);
    if(pair.first!=-1 && pair.second!=-1){
        TLorentzVector jet1(0,0,0,0);
        TLorentzVector jet2(0,0,0,0);

        jet1.SetPtEtaPhiE(f["hltpfjetPt"][pair.first],f["hltpfjetEta"][pair.first],f["hltpfjetPhi"][pair.first],f["hltpfjetE"][pair.first]);
        jet2.SetPtEtaPhiE(f["hltpfjetPt"][pair.second],f["hltpfjetEta"][pair.second],f["hltpfjetPhi"][pair.second],f["hltpfjetE"][pair.second]);
        dijet=jet1+jet2;

        //std::cout<<"dijet mass "<<dijet.M()<<std::endl;
    }
    
    *b["HLT_Ele20_WP85_DijetPT50"]  = dijet.Pt()>50 && HLT_Ele20;
    *b["HLT_Ele20_WP85_DijetPT100"] = dijet.Pt()>100 && HLT_Ele20;
    *b["HLT_Ele20_WP85_DijetPT150"] = dijet.Pt()>150 && HLT_Ele20;
    
    *b["HLT_Ele30_WP85_DijetPT50"]  = false;
    *b["HLT_Ele30_WP85_DijetPT100"] = false;
    *b["HLT_Ele30_WP85_DijetPT150"] = false;
    if(HLT_Ele20){
        if(f["hltElePt"][eg20ind]>30){
            *b["HLT_Ele30_WP85_DijetPT50"]  = dijet.Pt()>50 && HLT_Ele20;
            *b["HLT_Ele30_WP85_DijetPT100"] = dijet.Pt()>100 && HLT_Ele20;
            *b["HLT_Ele30_WP85_DijetPT150"] = dijet.Pt()>150 && HLT_Ele20;
        }
    }


    ofile->cd();
    outputTree->Fill();
    return;
}

void VHbbTrigger::TermAnalysis(){
    if(debug>10) std::cout<<"START TermAnalysis()"<<std::endl;
    ofile->cd();
    outputTree->Write();
    ofile->Close();
    if(debug>10) std::cout<<"DONE TermAnalysis()"<<std::endl;
    return;
}



std::pair<int,int> VHbbTrigger::HighestPtJets(float jetPtCut, float jetEtaCut, float vetoPhi, float vetoEta){
    std::pair<int,int> pair(-1,-1);


    for(int ijet=0; ijet<*in["nl1jets"]; ijet++){
        if(debug>100000) std::cout<<"bjet search at L1 "<<ijet<<std::endl;
        if(f["l1jetPt"][ijet]>jetPtCut && fabs(f["l1jetEta"][ijet])<jetEtaCut) {
            if(debug>100000) std::cout<<"pass pt cut "<<f["l1jetPt"][ijet]<<std::endl;
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["l1jetEta"][ijet],f["l1jetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["l1jetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.first == -1 ) {
                pair.first = ijet;
            } else if(f["l1jetPt"][pair.first]<f["l1jetPt"][ijet]){
                pair.first = ijet;
            }
            if(debug>100000) std::cout<<"jet1 "<<pair.first<<std::endl;
        }
    }

    for(int ijet=0; ijet<*in["nl1jets"]; ijet++){
        if(debug>100000) std::cout<<"2nd bjet search at L1 "<<ijet<<std::endl;
        if(f["l1jetPt"][ijet]>jetPtCut && fabs(f["l1jetEta"][ijet])<jetEtaCut) {
            if(debug>100000) std::cout<<"pass pt cut "<<f["l1jetPt"][ijet]<<std::endl;
            if(ijet==pair.first) continue;
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["l1jetEta"][ijet],f["l1jetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["l1jetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.second == -1 ) {
                pair.second = ijet;
            } else if(f["l1jetPt"][pair.second]<f["l1jetPt"][ijet]){
                pair.second = ijet;
            }
            if(debug>100000) std::cout<<"jet1 jet2 "<<pair.first<<" "<<pair.second<<std::endl;
        }
    }

    return pair;
}


std::pair<int,int> VHbbTrigger::HighestPtPFJets(float jetPt, float jetEta, float vetoPhi, float vetoEta){
    std::pair<int,int> pair(-1,-1);


    for(int ijet=0; ijet<*in["nhltpfjets"]; ijet++){
        if(debug>100000) std::cout<<"jet search at HLT "<<ijet<<std::endl;
        if(f["hltpfjetPt"][ijet]>jetPt && fabs(f["hltpfjetEta"][ijet])<jetEta ){
            if(debug>100000) std::cout<<"pass pt cut "<<f["hltpfjetPt"][ijet]<<std::endl;
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["hltpfjetEta"][ijet],f["hltpfjetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["hltpfjetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.first == -1 ) {
                pair.first = ijet;
            } else if(f["hltpfjetPt"][pair.first]<f["hltpfjetPt"][ijet]){
                pair.first = ijet;
            }
            if(debug>100000) std::cout<<"jet1 "<<pair.first<<std::endl;
        }
    }

    for(int ijet=0; ijet<*in["nhltpfjets"]; ijet++){
        if(debug>100000) std::cout<<"2nd jet search at HLT "<<ijet<<std::endl;
        if(f["hltpfjetPt"][ijet]>jetPt && fabs(f["hltpfjetEta"][ijet])<jetEta ){
            if(debug>100000) std::cout<<"pass pt cut "<<f["hltpfjetPt"][ijet]<<std::endl;
            if(ijet==pair.first) continue;
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["hltpfjetEta"][ijet],f["hltpfjetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["hltpfjetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.second == -1 ) {
                pair.second = ijet;
            } else if(f["hltpfjetPt"][pair.first]<f["hltpfjetPt"][ijet]){
                pair.first = ijet;
            }

            if(debug>100000) std::cout<<"jet1 jet2 "<<pair.first<<" "<<pair.second<<std::endl;
        }
    }

    return pair;
}


bool VHbbTrigger::PassEGL1(float ptCut, bool isolated, bool etaRestricted, int* L1ElInd, float* L1ElPhi){

    if(debug>100000) std::cout<<"PassEGL1 ptCut isolated etaRestricted "<<ptCut<<" "<<isolated<<" "<<etaRestricted<<std::endl;
    bool pass=false;
    int selInd=-1;
    int selPt=0;
    int selPhi=-10;


    if(debug>100000) std::cout<<"finding highest pt isolated electron"<<std::endl;
    for(int iel=0; iel<*in["nl1eleIso"]; iel++){
        if(f["l1eleIsoPt"][iel] < ptCut ) continue;
        if(etaRestricted && fabs(f["l1eleIsoEta"][iel])>2.17) continue;

        if(selInd==-1 || (selPt<f["l1eleIsoPt"][iel])) {
            selInd=iel;
            selPt=f["l1eleIsoPt"][iel];
            selPhi=f["l1eleIsoPhi"][iel];
        }
    }

    //if(selInd!=-1){
    //    *in["elIsoInd"]=selInd;
    //    *f["elIsoPt"]=selPt;
    //}
     
    if(!isolated){
        //*in["elInd"]=-1;
        //*f["elPt"]=-20;

        if(debug>100000) std::cout<<"finding highest pt non-isolated electron"<<std::endl;
        for(int iel=0; iel<*in["nl1eleNonIso"]; iel++){
            if(f["l1eleNonIsoPt"][iel] < ptCut ) continue;
            if(etaRestricted && fabs(f["l1eleNonIsoEta"][iel])>2.17) continue;
            if(selInd==-1 || (selPt<f["l1eleNonIsoPt"][iel])) {
                selInd=iel;
                selPt=f["l1eleNonIsoPt"][iel];
                selPhi=f["l1eleNonIsoPhi"][iel];
            }
        }
        //if(selInd!=-1){
        //    *in["elInd"]=selInd;
        //    *f["elPt"]=selPt;
        //}
    }
    
   
    if(L1ElInd!=0) *L1ElInd=selInd;
    if(L1ElPhi!=0) *L1ElPhi=selPhi;

    pass = selInd!=-1;
    return pass;
}


bool VHbbTrigger::PassEGPlusJetL1(float EGPtCut, bool EGiso, bool EGer, float jetPtCut, float jetEtaCut){
    int L1ElInd=-99;
    float L1ElPhi=-10;
    bool passEGL1 = PassEGL1(EGPtCut, EGiso, EGer, &L1ElInd, &L1ElPhi);

    if(debug>100000) std::cout<<"passing? "<<passEGL1<<" EG pt iso er "<<EGPtCut<<" "<<EGiso<<" "<<EGer<<" "<<L1ElInd<<std::endl;
    
    if(!passEGL1) return false;  

    // FIXME only works for isolated L1 seeds right now
    
    std::pair<int,int> l1jets=HighestPtJets(jetPtCut,jetEtaCut,L1ElPhi);
    if(debug>100000) std::cout<<"jet1 jet2 "<<l1jets.first<<" "<<l1jets.second<<std::endl;
  
    bool passJet1=false;
    if(l1jets.first!=-1){
        if(debug>100000) std::cout<<"found a jet with pt "<<f["l1jetPt"][l1jets.first]<<"  cut "<<jetPtCut<<std::endl;
        if(f["l1jetPt"][l1jets.first]>jetPtCut){
            if(debug>100000) std::cout<<"found a jet with eta "<<f["l1jetEta"][l1jets.first]<<"  cut "<<jetEtaCut<<std::endl;
            if(fabs(f["l1jetEta"][l1jets.first])<jetEtaCut){
                if(debug>100000) std::cout<<"passing!"<<std::endl;
                passJet1=true;
            }
        }
    }
    return passJet1;
}


bool VHbbTrigger::PassEleHLT(float ptCut, std::string WP, std::vector<std::string> L1seeds, bool etaRestricted, int* elInd){
    bool passHLT=false;
    bool passL1=false;


    for(unsigned iL1=0; iL1<L1seeds.size(); iL1++){
        if(L1seeds[iL1]=="noseed"){
            passL1=true;
            if(debug>100000) std::cout<<"no seed required--passing l1"<<std::endl;
            break;
        }
        if(debug>100000){
            std::cout<<"checking "<<iL1;
            std::cout<<" "<<L1seeds[iL1]<<std::endl;
        }
        if(*b[L1seeds[iL1]]==true){
            if(debug>100000){
                std::cout<<"passing! "<<iL1;
                std::cout<<" "<<L1seeds[iL1]<<std::endl;
            }
            passL1=true;
            break;
        }
    }

    if(!passL1) return passHLT;

    int selInd=-1;
    float selPt=-1;

    for(int iEl=0; iEl<*in["nhltele"]; iEl++){
        if(f["hltElePt"][iEl]<ptCut) continue;
        if(debug>100000)  std::cout<<iEl<<" with "<<f["hltElePt"][iEl]<<std::endl;
        if(etaRestricted && fabs(f["hltEleEta"][iEl])>2.17) continue;
        if(debug>100000)  std::cout<<iEl<<" with "<<f["hltEleEta"][iEl]<<std::endl;
        if(*f["recomputeWPs"]==1){
            if(WP=="WP85" && b["hltEleWP85ReComp"][iEl]==false) continue;
            if(WP=="WP75" && b["hltEleWP75ReComp"][iEl]==false) continue;
        } else {
            if(WP=="WP85" && b["hltEleWP85"][iEl]==false) continue;
            if(WP=="WP75" && b["hltEleWP75"][iEl]==false) continue;
        }
        if(debug>100000)  std::cout<<"WP75 WP85  "<<b["hltEleWP75ReComp"][iEl]<<" "<<b["hltEleWP85ReComp"][iEl]<<std::endl;
         
        if(selInd==-1 || selPt< f["hltElePt"][iEl]){
            selInd=iEl;
            selPt=f["hltElePt"][iEl];
            break;
        }
    }

    passHLT=(selInd!=-1);
    if(elInd) *elInd=selInd;

    return passHLT;
}


bool VHbbTrigger::FindBestCSVHLTJet(float ptCut, float etaCut, float CSVCut, std::string discriminator, int* jetInd){
    int selInd=-1;
   
    if(debug>100000000) std::cout<<"nCSVhltPF  "<<*in["nCSVhltPF"]<<std::endl;
    for(int ibjet=0; ibjet<*in["nCSVhltPF"]; ibjet++){
        if(f["CSVhltPF_jetPt"][ibjet]>ptCut){
            if(fabs(f["CSVhltPF_jetEta"][ibjet])<etaCut){
                if(debug>100000000) std::cout<<"CSVhltPF_jetCSV "<<f["CSVhltPF_jetCSV"][ibjet]<<std::endl;
                if(f["CSVhltPF_jetCSV"][ibjet]>CSVCut){
                    selInd=ibjet;
                    if(discriminator=="PT") ptCut=f["CSVhltPF_jetPt"][ibjet];
                    if(discriminator=="CSV") CSVCut=f["CSVhltPF_jetCSV"][ibjet];
                }
            }
        }
    }

    if(jetInd) *jetInd=selInd;

    return selInd!=-1;
}



bool VHbbTrigger::PassElePlusJetHLT(std::vector<std::string> L1seeds, float elPtCut, std::string elWP, bool elEtaRestricted, float jetPtCut, float jetEtaCut, int* elInd, int* jetInd){ 
    int hltElInd=-1;
    if(debug>10000) std::cout<<"look for hlt electron"<<std::endl;
    bool passL1 = false;
    if(debug>10000) std::cout<<"size of l1 seeds "<<L1seeds.size()<<std::endl;
    for(unsigned iL1=0; iL1<L1seeds.size(); iL1++){
        if(debug>10000) std::cout<<"iL1 "<<iL1<<" "<<L1seeds[iL1]<<std::endl;
        if(*b[L1seeds[iL1]]) {
            passL1=true;
            break;
        }
    }

    if(!passL1) return false;
    std::vector<std::string> seeds;
    seeds.clear();
    seeds.push_back("noseed");
    bool passEleHLT = PassEleHLT(elPtCut, elWP, seeds, elEtaRestricted, &hltElInd);
    if(debug>10000) std::cout<<"found electron "<<hltElInd<<" look for hlt jets"<<std::endl;
    if(!passEleHLT) return false;
    std::pair<int,int> pair=HighestPtPFJets(jetPtCut,jetEtaCut,f["hltElePhi"][hltElInd]);

    if(pair.first==-1) {
        return false;
    } else {
        if(elInd) *elInd=hltElInd;
        if(jetInd) *jetInd=pair.first;
    }
     
    return true; 
}


void VHbbTrigger::RecomputeWPs(){

    bool printEvent=false;

    if(debug>1000000) printEvent=true;

    for(int ihlt=0; ihlt<*in["nhltele"]; ihlt++){
        b["hltEleWP75ReComp"][ihlt]=false;
        b["hltEleWP85ReComp"][ihlt]=false;
        if(fabs(f["hltEleEta"][ihlt])<1.5) { //EB cuts
            if(printEvent) std::cout<<"EB el PT "<<f["hltElePt"][ihlt]<<std::endl;
            if(f["hltEleSieie"][ihlt]<0.011){
            if(printEvent) std::cout<<"pass sieie"<<std::endl;
            if(f["hltEleECalIso"][ihlt]/f["hltElePt"][ihlt]<0.21){
            if(printEvent) std::cout<<"pass ecaliso/pt"<<std::endl;
            if(f["hltEleHCalIso"][ihlt]/f["hltElePt"][ihlt]<0.11){
            if(printEvent) std::cout<<"pass hcaliso/pt"<<std::endl;
            if(f["hltEleOneOESuperMinusOneOP"][ihlt]<0.032){
            if(printEvent) std::cout<<"pass 1/e-1/p"<<std::endl;
            if(f["hltEleChi2"][ihlt]<3){
            if(printEvent) std::cout<<"pass chi2"<<std::endl;
            if(f["hltEleDEta"][ihlt]<0.0035){
            if(printEvent) std::cout<<"pass DEta"<<std::endl;
            if(f["hltEleDPhi"][ihlt]<0.021){
            if(printEvent) std::cout<<"pass DPhi"<<std::endl;
            if(f["hltEleTrkIso"][ihlt]/f["hltElePt"][ihlt]<0.05){
            if(printEvent) std::cout<<"pass tkiso/pt"<<std::endl;
            if(printEvent) std::cout<<"pass all, but h-0.01*e "<<f["hltEleHoE"][ihlt]<<"*"<<f["hltEleE"][ihlt]<<" - 0.01*"<<f["hltEleE"][ihlt]<<"="<<f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<<std::endl;
            //if(f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<4)
            if(f["hltEleHoE"][ihlt]<4){
            if(printEvent) std::cout<<"pass h-0.01*e "<<f["hltEleHoE"][ihlt]<<"*"<<f["hltEleE"][ihlt]<<" - 0.01*"<<f["hltEleE"][ihlt]<<"="<<f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<<std::endl;
                b["hltEleWP75ReComp"][ihlt]=true;
            }}}}}}}}}

            if(f["hltEleSieie"][ihlt]<0.011){
            if(f["hltEleECalIso"][ihlt]/f["hltElePt"][ihlt]<0.16){
            if(f["hltEleHCalIso"][ihlt]/f["hltElePt"][ihlt]<0.20){
            if(f["hltEleOneOESuperMinusOneOP"][ihlt]<0.012){
            if(in["hltEleMissingHits"][ihlt]<2){
            if(f["hltEleDEta"][ihlt]<0.005){
            if(f["hltEleDPhi"][ihlt]<0.03){
            if(f["hltEleTrkIso"][ihlt]/f["hltElePt"][ihlt]<0.05){
            if(printEvent) std::cout<<"pass all but h/e?"<<std::endl;
            //if(f["hltEleHoE"][ihlt]<0.15)
            if(f["hltEleHoE"][ihlt]/f["hltEleE"][ihlt]+0.01<0.15){
            if(printEvent) std::cout<<"pass all"<<std::endl;
                b["hltEleWP85ReComp"][ihlt]=true;
            }}}}}}}}}
       } else {  //EE cuts
            if(printEvent) std::cout<<"EE el PT "<<f["hltElePt"][ihlt]<<std::endl;
            if(f["hltEleSieie"][ihlt]<0.031){
            if(printEvent) std::cout<<"pass sieie"<<std::endl;
            if(f["hltEleECalIso"][ihlt]/f["hltElePt"][ihlt]<0.14){
            if(printEvent) std::cout<<"pass ecaliso/pt"<<std::endl;
            if(f["hltEleHCalIso"][ihlt]/f["hltElePt"][ihlt]<0.21){
            if(printEvent) std::cout<<"pass hcaliso/pt"<<std::endl;
            if(f["hltEleOneOESuperMinusOneOP"][ihlt]<0.032){
            if(printEvent) std::cout<<"pass 1/e-1/p"<<std::endl;
            if(f["hltEleDEta"][ihlt]<0.0065){
            if(printEvent) std::cout<<"pass DEta"<<std::endl;
            if(f["hltEleDPhi"][ihlt]<0.035){
            if(printEvent) std::cout<<"pass DPhi"<<std::endl;
            if(f["hltEleTrkIso"][ihlt]/f["hltElePt"][ihlt]<0.05){
            if(printEvent) std::cout<<"pass tkiso/pt"<<std::endl;
            //if(in["hltEleMissingHits"][ihlt]<1)
            if(printEvent) std::cout<<"pass missing hits"<<std::endl;
            if(f["hltEleChi2"][ihlt]<2.8){
            if(printEvent) std::cout<<"pass chi2"<<std::endl;
            if(printEvent) std::cout<<"pass? h-0.01*e "<<f["hltEleHoE"][ihlt]<<"*"<<f["hltEleE"][ihlt]<<" - 0.01*"<<f["hltEleE"][ihlt]<<"="<<f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<<std::endl;
            //if(f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<13)
            if(f["hltEleHoE"][ihlt]<13){
            if(printEvent) std::cout<<"pass h-0.01*e "<<f["hltEleHoE"][ihlt]<<"*"<<f["hltEleE"][ihlt]<<" - 0.01*"<<f["hltEleE"][ihlt]<<"="<<f["hltEleHoE"][ihlt]*f["hltEleE"][ihlt] - 0.01*f["hltEleE"][ihlt]<<std::endl;
                b["hltEleWP75ReComp"][ihlt]=true;
                if(printEvent) std::cout<<"ihlt hltEleWP75ReComp "<<ihlt<<" "<<b["hltEleWP75ReComp"][ihlt]<<std::endl;
            }}}}}}}}}

            if(f["hltEleSieie"][ihlt]<0.033){
            if(f["hltEleECalIso"][ihlt]/f["hltElePt"][ihlt]<0.12){
            if(f["hltEleHCalIso"][ihlt]/f["hltElePt"][ihlt]<0.30){
            if(f["hltEleOneOESuperMinusOneOP"][ihlt]<0.009){
            if(in["hltEleMissingHits"][ihlt]<2){
            if(f["hltEleDEta"][ihlt]<0.010){
            if(f["hltEleDPhi"][ihlt]<0.03){
            if(f["hltEleTrkIso"][ihlt]/f["hltElePt"][ihlt]<0.05){
            if(printEvent) std::cout<<"pass all but h/e?"<<std::endl;
            if(f["hltEleHoE"][ihlt]/f["hltEleE"][ihlt]+0.01<0.20){
            if(printEvent) std::cout<<"pass all"<<std::endl;
                b["hltEleWP85ReComp"][ihlt]=true;
            }}}}}}}}}
        }
    }
}
    
void VHbbTrigger::GenMatchElectrons(){
    if(debug>10000000) std::cout<<"nhltele "<<*in["nhltele"]<<std::endl;
    *f["hltEleBestMatchDR"]=9999;
    *in["hltEleBestMatch"]=-1;
    for(int ihlt=0; ihlt<*in["nhltele"]; ihlt++){
        if(debug>10000000) {
            std::cout<<"ihlt "<<ihlt<<std::endl;
            std::cout<<"eta phi "<<f["hltEleEta"][ihlt]<<" "<<f["hltElePhi"][ihlt]<<std::endl;
        }
        if(*in["sampleIndex"]<0){
            f["hltEleMatchDR"][ihlt]=EvalDeltaR(f["hltEleEta"][ihlt],f["hltElePhi"][ihlt],*f["genEleEta"],*f["genElePhi"]);
            b["hltEleMatch"][ihlt]=(f["hltEleMatchDR"][ihlt]<0.15);

            if(f["hltEleMatchDR"][ihlt]<*f["hltEleBestMatchDR"]){
                *f["hltEleBestMatchDR"]=f["hltEleMatchDR"][ihlt];
                *in["hltEleBestMatch"]=ihlt;
            }
        } else {
            f["hltEleMatchDR"][ihlt]=9999;
            b["hltEleMatch"][ihlt]=0;
        }
    }
}
