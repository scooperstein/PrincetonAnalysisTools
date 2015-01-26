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
        if( *f["b1Pt_gen"]<*f["genJ1PtCut"]
            ||*f["b2Pt_gen"]<*f["genJ2PtCut"]
            ||*f["genElePt"]<*f["genElPtCut"]){
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

    *in["elIsoInd"]=-1;
    *f["elIsoPt"]=-1;
    for(int iel=0; iel<*in["neleIso"]; iel++){
        if(f["eleIsoPt"][iel] < 10 || fabs(f["eleIsoEta"][iel]) > 2.4) continue;
        //if(EvalDeltaR(f["jetEta"][jets.first],f["jetPhi"][jets.first],f["eleIsoEta"][iel],f["eleIsoPhi"][iel])<0.5) continue;
        //if(EvalDeltaR(f["jetEta"][jets.second],f["jetPhi"][jets.second],f["eleIsoEta"][iel],f["eleIsoPhi"][iel])<0.5) continue;
        if(*in["elIsoInd"]==-1){
            *in["elIsoInd"]=iel;
            *f["elIsoPt"]=f["eleIsoPt"][iel];
        } else if(*f["elIsoPt"]<f["eleIsoPt"][iel]) {
            *in["elIsoInd"]=iel;
            *f["elIsoPt"]=f["eleIsoPt"][iel];
        }
    }
    
     
    if(debug>1000) {
        std::cout<<"selecting bjets"<<std::endl;
    }
    std::pair<int,int> jets= (*in["elIsoInd"]!=-1)?HighestPtJets(f["eleIsoPhi"][*in["elIsoInd"]],f["eleIsoEta"][*in["elIsoInd"]]):HighestPtJets();
   
    // there aren't two acceptable jets
    if(*f["preselNJets"]==2 && (jets.first==-1 || jets.second==-1)) return sel;
    if(*f["preselNJets"]==1 && jets.first==-1) return sel;
    *in["jetInd1"]=jets.first;
    *in["jetInd2"]=jets.second;

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
        
        
    if(     *f["elIsoPt"] <  *f["elIsoPtCut"]  
        &&  *f["muPt"] <  *f["muPtCut"]  
        &&  *f["mhtPt"] <  *f["metPtCut"]  ) {
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
    

    if(*in["elIsoInd"]!=-1){
        *f["elIsoEta"]=f["eleIsoEta"][*in["elIsoInd"]];
        *f["elIsoPhi"]=f["eleIsoPhi"][*in["elIsoInd"]];
        *f["elIsoIso"]=f["eleIsoIso"][*in["elIsoInd"]];
        *f["dphijet1e"]=EvalDeltaPhi(*f["elIsoPhi"],f["jetPhi"][jets.first]);
        *f["dr_eiso"]=EvalDeltaR(f["eleIsoEta"][*in["elIsoInd"]],f["eleIsoPhi"][*in["elIsoInd"]],*f["genEleEta"],*f["genElePhi"]);
    } else {
        *f["dphijet1e"]=-20;
        *f["elIsoEta"]=-20;
        *f["elIsoPhi"]=-20;
        *f["elIsoIso"]=-20;
    }    

    if(*f["dphijet1e"]<0 && *f["dphijet1e"]>-19) {
        std::cout<<"e j dphi "<<*f["elIsoPhi"]<<" "<<f["jetPhi"][jets.first]<<" "<<*f["dphijet1e"]<<std::endl;
    }
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
            if( pair.second == -1 ) {
                pair.second = ijet;
            } else if(f["jetPt"][pair.second]<f["jetPt"][ijet]){
                pair.first = ijet;
            }
        }
    }

    return pair;
}

