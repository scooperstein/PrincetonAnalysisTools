//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "VHbbTrigger.h"

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

    int npreseljets=0;
    for(int ijet=0; ijet<*in["njets"]; ijet++){
        if(f["jetPt"][ijet] > *f["jetPtPreCut"]) npreseljets++;
    }
    if(npreseljets>=(int)*f["preselNJets"]) sel=true;

    return sel;
}

bool VHbbTrigger::Analyze(){
    bool sel=false;

    if(debug>100000) std::cout<<"finding highest pt isolated electron"<<std::endl;
    *in["elIsoInd"]=-1;
    *f["elIsoPt"]=-20;
    for(int iel=0; iel<*in["neleIso"]; iel++){
        if(f["eleIsoPt"][iel] < *f["elIsoPtCut"] ) continue;
        if(*in["elIsoInd"]==-1){
            *in["elIsoInd"]=iel;
            *f["elIsoPt"]=f["eleIsoPt"][iel];
        } else if(*f["elIsoPt"]<f["eleIsoPt"][iel]) {
            *in["elIsoInd"]=iel;
            *f["elIsoPt"]=f["eleIsoPt"][iel];
        }
    }
    
     
    if(debug>100000) std::cout<<"finding highest pt non-isolated electron"<<std::endl;
    *in["elInd"]=-1;
    *f["elPt"]=-20;
    for(int iel=0; iel<*in["neleNonIso"]; iel++){
        if(f["eleNonIsoPt"][iel] < *f["elPtCut"] ) continue;
        if(*in["elInd"]==-1){
            *in["elInd"]=iel;
            *f["elPt"]=f["eleNonIsoPt"][iel];
        } else if(*f["elPt"]<f["eleNonIsoPt"][iel]) {
            *in["elInd"]=iel;
            *f["elPt"]=f["eleNonIsoPt"][iel];
        }
    }
   
    if(debug>100000) std::cout<<"setting some iso electron branches"<<std::endl;
    if(*in["elIsoInd"]!=-1){
        *f["elIsoEta"]=f["eleIsoEta"][*in["elIsoInd"]];
        *f["elIsoPhi"]=f["eleIsoPhi"][*in["elIsoInd"]];
        *f["dr_eiso"]=EvalDeltaR(f["eleIsoEta"][*in["elIsoInd"]],f["eleIsoPhi"][*in["elIsoInd"]],*f["genEleEta"],*f["genElePhi"]);
    } else {
        *f["elIsoEta"]=-20;
        *f["elIsoPhi"]=-20;
        *f["dr_eiso"]=-20;
    }    

    if(debug>100000) std::cout<<"setting some non-iso electron branches"<<std::endl;
    if(*in["elInd"]!=-1){
        *f["elEta"]=f["eleNonIsoEta"][*in["elInd"]];
        *f["elPhi"]=f["eleNonIsoPhi"][*in["elInd"]];
        *f["dr_e"]=EvalDeltaR(f["eleNonIsoEta"][*in["elInd"]],f["eleNonIsoPhi"][*in["elInd"]],*f["genEleEta"],*f["genElePhi"]);
    } else {
        *f["dphijet1e"]=-20;
        *f["elEta"]=-20;
        *f["elPhi"]=-20;
        *f["dr_e"]=-20;
    }    

    if(debug>100000) std::cout<<"re-setting some non-iso electron branches"<<std::endl;
    if(*in["elIsoInd"]!=-1 && (*f["elPt"]<*f["elIsoPt"]) && *f["elIsoPt"]>*f["elPtCut"]){
        *f["elEta"]=    *f["elIsoEta"];
        *f["elPhi"]=    *f["elIsoPhi"];
        *f["elPt"]=     *f["elIsoPt"];
        *f["dr_e"]=     *f["dr_eiso"];
    }

     
    if(debug>1000) {
        std::cout<<"selecting bjets"<<std::endl;
    }
    std::pair<int,int> jets= (*f["elPt"]>0)?HighestPtJets(*f["elPhi"]):HighestPtJets();
    //std::pair<int,int> jets= (*f["elPt"]>0)?HighestPtJets(*f["elPhi"],*f["elEta"]):HighestPtJets();
    //std::pair<int,int> jets= (*in["elIsoInd"]!=-1)?HighestPtJets(f["eleIsoPhi"][*in["elIsoInd"]],f["eleIsoEta"][*in["elIsoInd"]]):HighestPtJets();
   
    // there aren't two acceptable jets
    if(*f["preselNJets"]==2 && (jets.first==-1 || jets.second==-1)) return sel;
    if(*f["preselNJets"]==1 && jets.first==-1) return sel;
    if(*f["preselNJets"]==0 && jets.first==-1) {
        //dummy to protect assumption of jet below
        jets.first=0;
    }
    *in["jetInd1"]=jets.first;
    *in["jetInd2"]=jets.second;

    if(*f["elPhi"]!=-20){
        *f["dphijet1e"]=EvalDeltaPhi(*f["elPhi"],f["jetPhi"][jets.first]);
    } else {
        *f["dphijet1e"]=-20;
    }

    if(*in["elIsoInd"]!=-1){
        *f["dphijet1eiso"]=EvalDeltaPhi(*f["elIsoPhi"],f["jetPhi"][jets.first]);
    } else {
        *f["dphijet1eiso"]=-20;
    }    
    
    //FIXME
    if(*in["elIsoInd"]!=-1 && (*f["elPt"]<*f["elIsoPt"]) && *f["elIsoPt"]>*f["elPtCut"]){
        *f["dphijet1e"]=*f["dphijet1eiso"];
    }

    if(debug>1000) {
        std::cout<<"found two bjets with pt ";
        
        if(jets.first!=-1) std::cout<<f["jetPt"][jets.first]<<" ";
        else std::cout<<" no first jet ";
        if(jets.second!=-1) std::cout<<f["jetPt"][jets.second]<<" ";
        else std::cout<<" no second jet ";
           
        std::cout<<std::endl;
    }

    if(jets.first!=-1) {
        *f["jetPt1"]=f["jetPt"][jets.first];    
        *f["jetPhi1"]=f["jetPhi"][jets.first];    
    } else {
        *f["jetPt1"]=-1;
        *f["jetPhi1"]=-1;
    }    

    if(jets.second!=-1) {
        *f["jetPt2"]=f["jetPt"][jets.second];    
        *f["jetPhi2"]=f["jetPhi"][jets.second];    
    } else {
        *f["jetPt2"]=-1;
        *f["jetPhi2"]=-1;
    }    
 
    // make composite of jets for classification
    // need jet energy or mass

    *in["muInd"]=-1;
    *f["muPt"]=-1;
    
    for(int imu=0; imu<*in["nmuons"]; imu++){
        if(f["muonPt"][imu] < 10 || fabs(f["muonEta"][imu]) > 2.4) continue;
        if(jets.first!=-1) {
            if(EvalDeltaR(f["jetEta"][jets.first],f["jetPhi"][jets.first],f["muonEta"][imu],f["muonPhi"][imu])<0.5) continue;
        }
        if(jets.second!=-1) {
            if(EvalDeltaR(f["jetEta"][jets.second],f["jetPhi"][jets.second],f["muonEta"][imu],f["muonPhi"][imu])<0.5) continue;
        }
        if(*in["muInd"]==-1){
            *in["muInd"]=imu;
            *f["muPt"]=f["muonPt"][imu];
        } else if(*f["muPt"]<f["muonPt"][imu]) {
            *in["muInd"]=imu;
            *f["muPt"]=f["muonPt"][imu];
        }
    }
        
    // could be wrong if elIsoPtCut>elPtCut      
    if(     (*f["elIsoPt"] <  *f["elIsoPtCut"]  
        || *in["elIsoInd"]==-1)
        &&  (*f["elPt"] <  *f["elPtCut"]  
        || *in["elInd"]==-1)){
        return sel;
    }
    sel=true;
   

    if(debug>1000) {
        std::cout<<"selecting event"<<std::endl;
    }
    
    if(jets.first!=-1) {
        *f["dr_jet1"]=std::min(EvalDeltaR(f["jetEta"][jets.first],f["jetPhi"][jets.first],*f["b1Eta_gen"],*f["b1Phi_gen"]),EvalDeltaR(f["jetEta"][jets.first],f["jetPhi"][jets.first],*f["b2Eta_gen"],*f["b2Phi_gen"]));
        *f["dphijet1met"]=EvalDeltaPhi(f["jetPhi"][jets.first],*f["metPhi"]);
    } else {
        *f["dr_jet1"]=-1;
        *f["dphijet1met"]=-20;
    }

    if(jets.second!=-1) {
        *f["dr_jet2"]=std::min(EvalDeltaR(f["jetEta"][jets.second],f["jetPhi"][jets.second],*f["b1Eta_gen"],*f["b1Phi_gen"]),EvalDeltaR(f["jetEta"][jets.second],f["jetPhi"][jets.second],*f["b2Eta_gen"],*f["b2Phi_gen"]));
    } else {
        *f["dr_jet2"]=-1;
    } 

    if(jets.first!=-1 && jets.second!=-1) {
        *f["dphijets"]=EvalDeltaPhi(f["jetPhi"][jets.first],f["jetPhi"][jets.second]);
    
        float dijetdPhi= f["jetPhi"][jets.first]-f["jetPhi"][jets.second];
        float dijetPhi= (f["jetPhi"][jets.first] + f["jetPhi"][jets.second])/2;
        if(abs(dijetdPhi)>PI){
            //std::cout<<"need to fix phijj "<<dijetPhi<<std::endl;
            if(dijetPhi>0){ 
                dijetPhi = dijetPhi-PI;
            } else{
                dijetPhi = dijetPhi+PI;
            }
            //std::cout<<"PI "<<PI<<std::endl;
            //std::cout<<"subtracted pi to "<<dijetPhi<<std::endl;
        } 

        if(dijetPhi>PI) {
            dijetPhi=dijetPhi-2*PI;
        }
        //    std::cout<<"phi1 phi2 phijj "<<f["jetPhi"][jets.first]<<" "<<f["jetPhi"][jets.second]<<" "<<dijetPhi<<std::endl;
        //}

        *f["phijj"]=dijetPhi;
        
    } //end dijet if
    

    if(*in["muInd"]!=-1){
        //*f["dphijetselmu"]=EvalDeltaPhi(dijetPhi,f["muonPhi"][*in["muInd"]]);
        *f["muEta"]=f["muonEta"][*in["muInd"]];
    } else {
        *f["muEta"]=-20;
    }    


    *in["eventClass"]=0;

    return sel;
}

void VHbbTrigger::FinishEvent(){
    *in["sampleIndex"] = cursample->sampleNum;
    *f["weight"] = cursample->intWeight;
  
    *b["pass_eleHLT"] = 
        ((*f["elPt"]>35 || *f["elIsoPt"]>30) 
        && *f["hltEleWP85Pt"]>35 && fabs(*f["hltEleWP85Eta"])<2.1); 
        //&& *f["hltElePt"]>35 && fabs(*f["hltEleEta"])<2.1); 
    *b["pass_elePlusMetHLT"] = 
        ((*f["elPt"]>35 || *f["elIsoPt"]>30 || *f["mhtPt"]>70) 
        && *f["hltEleWP85Pt"]>25 && fabs(*f["hltEleWP85Eta"])<2.1
        //&& *f["hltElePt"]>25 && fabs(*f["hltEleEta"])<2.1
        && *f["pfmet"]>80 ); 
   
    *b["passCurrentL1"] = 
        (*f["elPt"]>35 || *f["elIsoPt"]>30 || *f["mhtPt"]>70); 
 
    *b["notPassBaseline"]=(!*b["pass_eleHLT"] && !*b["pass_elePlusMetHLT"]);
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



std::pair<int,int> VHbbTrigger::HighestPtJets(float vetoPhi, float vetoEta){
    std::pair<int,int> pair(-1,-1);

    for(int ijet=0; ijet<*in["njets"]; ijet++){
        if(f["jetPt"][ijet]>*f["j1ptCut"]) {
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["jetEta"][ijet],f["jetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["jetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.first == -1 ) {
                pair.first = ijet;
            } else if(f["jetPt"][pair.first]<f["jetPt"][ijet]){
                pair.first = ijet;
            }
        }
    }

    for(int ijet=0; ijet<*in["njets"]; ijet++){
        if(f["jetPt"][ijet]>*f["j2ptCut"]) {
            if(ijet==pair.first) continue;
            if(vetoEta!=-10 && vetoPhi!=-10){
                if(EvalDeltaR(f["jetEta"][ijet],f["jetPhi"][ijet],vetoEta,vetoPhi)<0.5) continue;
            } else if(vetoEta==-10 && vetoPhi!=-10){
                if(EvalDeltaPhi(f["jetPhi"][ijet],vetoPhi)<0.3) continue;
            }
            if( pair.second == -1 ) {
                pair.second = ijet;
            } else if(f["jetPt"][pair.second]<f["jetPt"][ijet]){
                pair.first = ijet;
            }
        }
    }

    return pair;
}

