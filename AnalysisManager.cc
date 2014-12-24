#define AnalysisManager_cxx
#include "AnalysisManager.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TMath.h>

AnalysisManager::AnalysisManager(){
}

void AnalysisManager::Initialize(const char* fileName) {
// used to generate this class and read the Tree.
    TChain *chain = new TChain("tree");
    chain->Add(fileName);
    InitChain(chain);

    ui.clear();
    in.clear();
    f.clear();
    d.clear();
    b.clear();
    
    branches.clear();
    branchInfos.clear();

    debug=0;
    outputTreeName="condensed_tree";
    safemode=1;
    
    return;
}


AnalysisManager::~AnalysisManager()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();

    for(std::map<std::string,unsigned int*>::iterator uiit=ui.begin();
            uiit!=ui.end();  ++uiit){
        if(debug>100) std::cout<<"I'm deleting "<<uiit->first<<std::endl;
        delete uiit->second;
    }

    for(std::map<std::string,int*>::iterator iit=in.begin();
            iit!=in.end();  ++iit){
        if(debug>100) std::cout<<"I'm deleting "<<iit->first<<std::endl;
        delete iit->second;
    }

    for(std::map<std::string,float*>::iterator fit=f.begin();
            fit!=f.end();  ++fit){
        if(debug>100) std::cout<<"I'm deleting "<<fit->first<<std::endl;
        delete fit->second;
    }

    for(std::map<std::string,double*>::iterator dit=d.begin();
            dit!=d.end();  ++dit){
        if(debug>100) std::cout<<"I'm deleting "<<dit->first<<std::endl;
        delete dit->second;
    }

    for(std::map<std::string,bool*>::iterator bit=b.begin();
            bit!=b.end();  ++bit){
        if(debug>100) std::cout<<"I'm deleting "<<bit->first<<std::endl;
        delete bit->second;
    }
}


void AnalysisManager::AddSample(SampleContainer sample){
    samples.push_back(sample);
}


Int_t AnalysisManager::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t AnalysisManager::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
   }
   return centry;
}


void AnalysisManager::InitChain(TChain *tree)
{
    // The InitChain() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    if(debug>100) std::cout<<"Resetting branches"<<std::endl;
    ResetBranches();

    // Set special branches
    fChain->SetBranchAddress("H", &H, &b_H);
    fChain->SetBranchAddress("V", &V, &b_V);
    fChain->SetBranchAddress("METtype1corr", &METtype1corr, &b_METtype1corr);
}


void AnalysisManager::SetupBranch(std::string name, int type, int length){
    branches[name] = new TBranch;
    branchInfos[name] = new BranchInfo(name,type,length,"existing");

    // Only 0-9 are setup with types for the moment.
    if(type>9 || type<0) {
        std::cout<<"Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }

    if(type==0) {
        ui[name] = new unsigned int;
        fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
    } else if(type==1) {
        in[name] = new int;
        fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
    } else if(type==2) {
        f[name] = new float;
        fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
    } else if(type==3) {
        d[name] = new double;
        fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
    } else if(type==4) {
        b[name] = new bool;
        fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
    }

    
    if(type>4 && type<10 && length<0) {
        std::cout<<"Types 5-9 are arrays and need a length greater than 0... not "
            <<length<<std::endl;
        return;
    }

    if(type==5) {
        ui[name] = new unsigned int[length];
        fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
    } else if(type==6) {
        in[name] = new int[length];
        fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
    } else if(type==7) {
        f[name] = new float[length];
        fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
    } else if(type==8) {
        d[name] = new double[length];
        fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
    } else if(type==9) {
        b[name] = new bool[length];
        fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
    }

    return;
}

// call after adding all pre-existing branches, but before adding any new branches
void AnalysisManager::ConfigureOutputTree() {
    outputTree = fChain->CloneTree(0);
}

void AnalysisManager::SetupNewBranch(std::string name, int type, int length, bool newmem){
    if(debug>1000) {
        std::cout<<"SetupNewBranch "<<name<<std::endl;
        std::cout<<"\ttype length newmem "<<type<<" "<<length<<" "<<newmem<<std::endl;
    }
    if(newmem) branchInfos[name] = new BranchInfo(name,type,length,"new");
    if(debug>1000) std::cout<<"BranchInfo instaniated"<<std::endl; 
       

    if(type>9 || type<0) {
        std::cout<<"New Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }

    if(type==0) {
        if(newmem) ui[name] = new unsigned int;
        branches[name] = outputTree->Branch(name.c_str(), ui[name], Form("%s/i",name.c_str()));
    } else if(type==1) {
        if(newmem) in[name] = new int;
        branches[name] = outputTree->Branch(name.c_str(), in[name], Form("%s/I",name.c_str()));
    } else if(type==2) {
        if(newmem) f[name] = new float;
        branches[name] = outputTree->Branch(name.c_str(), f[name], Form("%s/F",name.c_str()));
    } else if(type==3) {
        if(newmem) d[name] = new double;
        branches[name] = outputTree->Branch(name.c_str(), d[name], Form("%s/D",name.c_str()));
    } else if(type==4) {
        if(newmem) b[name] = new bool;
        branches[name] = outputTree->Branch(name.c_str(), b[name], Form("%s/O",name.c_str()));
    }


    if(type>4 && type<10 && length<0) {
        std::cout<<"Types 5-9 are arrays and need a length greater than 0... not "
            <<length<<std::endl;
        return;
    }

    if(type==5) {
        if(newmem) ui[name] = new unsigned int[length];
        branches[name] = outputTree->Branch(Form("%s[%i]",name.c_str(),length), ui[name], Form("%s/i",name.c_str()));
    } else if(type==6) {
        if(newmem) in[name] = new int[length];
        branches[name] = outputTree->Branch(Form("%s[%i]",name.c_str(),length), in[name], Form("%s/I",name.c_str()));
    } else if(type==7) {
        if(newmem) f[name] = new float[length];
        branches[name] = outputTree->Branch(Form("%s[%i]",name.c_str(),length), f[name], Form("%s/F",name.c_str()));
    } else if(type==8) {
        if(newmem) d[name] = new double[length];
        branches[name] = outputTree->Branch(Form("%s[%i]",name.c_str(),length), d[name], Form("%s/D",name.c_str()));
    } else if(type==9) {
        if(newmem) b[name] = new bool[length];
        branches[name] = outputTree->Branch(Form("%s[%i]",name.c_str(),length), b[name], Form("%s/O",name.c_str()));
    }

    return;
}

void AnalysisManager::ResetBranches(){
    
    for(std::map<std::string,BranchInfo*>::iterator biit=branchInfos.begin();
            biit!=branchInfos.end(); ++biit) {
        if(debug>100) std::cout<<"Branch "<<biit->second->name<<" of type, prov "<<biit->second->type<<" "<<biit->second->prov<<std::endl;
        if(biit->second->prov == "existing"){
            std::string name(biit->first);
            if(biit->second->type>9 || biit->second->type<0){
                std::cout<<"Branch "<<name<<" of unknown type "<<biit->second->type<<std::endl;
                continue;
            }
            if(biit->second->type%5==0) {
                fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
            } else if(biit->second->type%5==1) {
                fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
            } else if(biit->second->type%5==2) {
                fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
            } else if(biit->second->type%5==3) {
                fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
            } else if(biit->second->type%5==4) {
                fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
            } 
        }
    }
   
}


void AnalysisManager::PrintBranches(){
    std::cout<<"Branches in branch map WHY"<<std::endl;
    for(std::map<std::string,TBranch*>::iterator ibranch=branches.begin(); 
            ibranch!=branches.end(); ++ibranch){
        std::cout<<ibranch->first<<std::endl;
    }
}


void AnalysisManager::SetBranches(){
    // only copy branches specified
    if(debug>10) std::cout<<"SetBranchStatus of existing branches"<<std::endl;
    fChain->SetBranchStatus("*", 0);
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin(); 
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "existing") {
            if(debug>100) std::cout<<"fChain->SetBranchStatus("<<ibranch->first.c_str()<<", 1);"<<std::endl;
            fChain->SetBranchStatus(ibranch->first.c_str(), 1);
        }
        if(ibranch->second->prov == "new") {
            
        }
    }
}

void AnalysisManager::SetNewBranches(){
    //  
    if(debug>10) std::cout<<"SetupNewBranches"<<std::endl;
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin(); 
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "new") {
            if(debug>100) std::cout<<"SetupNewBranch "<<ibranch->first.c_str()<<std::endl;
            SetupNewBranch(ibranch->first, ibranch->second->type, ibranch->second->length, false); //newmem is false
        }
    }
}


void AnalysisManager::Loop(){

    if(debug>10) std::cout<<"Starting Loop"<<std::endl;
    SetBranches();
    
    // add new branches
    ofile = new TFile(Form("%s.root",outputTreeName.c_str()),"recreate");
    ofile->cd();
    //TTree *newtree = fChain->CloneTree(0);
    // for now we will use a default file to set the structure for the output tree
    outputTree = fChain->CloneTree(0); // need this one to only keep the branches you want 
    
    // it would probably be smarter to check also if SetupBranch() has been called since the last time outputTree was initialized
    //if (!outputTree) {
    //    std::cout<<"Did not call ConfigureOutputTree() before Loop(). Assuming there are no new branches..."<<std::endl;
    //    outputTree = fChain->CloneTree(0);
    //}  
    // let's add some of our own branches
    SetNewBranches();
    
    
    if(debug>10) std::cout<<"Done setting up branches; about to Init"<<std::endl;
    // eventually we'd like to set up the identifiers unique to each BDT as members of the class. These identifiers are:
    // 1. BDT name
    // 2. List of names of variables used to train BDT
    // 3. List of references to variables used to train BDT
    // 4. Reference to BDT?
    // ...etc
    InitAnalysis();
    
    
    if(debug>10) std::cout<<"About to loop over samples"<<std::endl;
    // loop through one sample at a time
    for (int i = 0; i < (int)samples.size(); i++) { 
        cursample = &samples[i];
        // set fChain to the TChain for the current sample
        InitChain(cursample->sampleChain);
    
        // FIXME should have a sample name, but doesn't right now
        //if(debug>100) std::cout<<"About to loop over events in "<<cursample->sampleName<<std::endl;
        // loop through the events
        Long64_t nentries = fChain->GetEntries();
        if(debug>1) std::cout<<"looping over "<<nentries<<std::endl;
        Long64_t nbytes = 0, nb = 0;
        // FIXME need a loop over systematics
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            bool presel = Preselection();
            if(!presel) continue;
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;
            if(ientry%10000==0 && debug>1)  std::cout<<"entry "<<ientry<<std::endl;
  
            bool select = Analyze();//SampleWenuAnalysis();
            if(select){
                FinishEvent();
            }
        } // end event loop
    } // end sample loop
    
    
    TermAnalysis();
    
    /*// Now let's mess around with re-evaluating the BDTs
    std::cout<<"Let's try to re-evaluate the BDTs..."<<std::endl;
    WriteBDTs("./", Form("%s.root",outputTreeName), "./", Form("%s_reevaluated.root",outputTreeName),"");
    */
}


void AnalysisManager::InitAnalysis(){
    if(debug>100) std::cout<<"InitAnalysis"<<std::endl;
}


bool AnalysisManager::Preselection(){
    bool sel=false;
    if(1) sel=true;
    return sel;
} 


bool AnalysisManager::Analyze(){
    bool sel=false;
    return sel;
}

            
void AnalysisManager::FinishEvent(){
    //need to fill the tree, hist, or whatever here.
    ofile->cd();
    outputTree->Fill();
    return;
}


void AnalysisManager::TermAnalysis() {
    // save tree here
    ofile->cd();
    outputTree->Write();
    ofile->Close();
}


// Set up all the BDT branches and configure the BDT's with the same input variables as used in training. Run before looping over events.
void AnalysisManager::SetupBDT( ) {

    // let's do 8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg as an example
    if(debug>100){
        thereader = new TMVA::Reader( "!Color:!Silent" );
    } else {
        thereader = new TMVA::Reader( "!Color:Silent" );
    }
    int bdtType = 0;
    if (bdtType == 0) {
        thereader->AddVariable("H.massCorr",                &testMass); 
        thereader->AddVariable("H.ptCorr",                  &H.ptCorr);
        thereader->AddVariable("V.pt",                      &V.pt);
        thereader->AddVariable("hJet_csvReshapedNew[0]",    &(f["hJet_csvReshapedNew"][0]) );
        thereader->AddVariable("hJet_csvReshapedNew[1]",    &(f["hJet_csvReshapedNew"][1]) );
        thereader->AddVariable("HVdPhi",                    f["HVdPhi"]);
        thereader->AddVariable("H.dEta",                    &H.dEta);
        thereader->AddVariable("Sum$(aJet_pt>20 && abs(aJet_eta)<4.5)", f["naJetsPassingCuts"]);
        thereader->AddVariable("H.dR",                      &H.dR);
        thereader->AddVariable("abs(deltaPullAngle)",       f["absDeltaPullAngle"]);
        thereader->AddVariable("hJet_ptCorr[0]",            &(f["hJet_ptCorr"][0]) );
        thereader->AddVariable("hJet_ptCorr[1]",            &(f["hJet_ptCorr"][1]) );
        thereader->AddVariable("MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)", f["highestCSVaJet"]);
        thereader->AddVariable("MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)", f["minDeltaRaJet"]);   
 
        thereader->AddSpectator("Vtype", f["Vtype_f"]);
        thereader->AddSpectator("Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )", f["naLeptonsPassingCuts"]);
        thereader->AddSpectator("METtype1corr.et", &METtype1corr.et);

        thereader->BookMVA("BDT method","aux/TMVA_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4_BDT.weights.xml");
        //thereader->BookMVA("BDT method",bdtInfo.xmlFile);
    }

    /*for (std::map<char*,Float_t&>::iterator it = bdtInfo.inputs.begin(); it != bdtInfo.inputs.end(); ++it) {
        thereader->AddVariable(it->first, it->second);
    }

    for (std::map<char*, Float_t&>::iterator it = bdtInfo.secInputs.begin(); it!= bdtInfo.secInputs.end(); ++it) {
        thereader->AddSpectator(it->first, it->second);
    }

    thereader->BookMVA(bdtInfo.method, bdtInfo.xmlFile);*/

}


void AnalysisManager::m(std::string key){
    if(branchInfos.count(key)==0){
        std::cout<<"There is no branch with name "<<key<<std::endl;
        if(safemode){
            std::cout<<"The program must be terminated..."<<std::endl;
            std::exit(0);
        } else {
            if(debug>1) std::cout<<"Returning 0 and hoping for the best."<<std::endl;
        }
    } else {
        if(debug>1000) std::cout<<"here is where I should return the right branch value"<<std::endl;
        switch(branchInfos[key]->type)
        {
        case 0:
            if(debug>1000) std::cout<<"unsigned int "<<*ui[key]<<std::endl; 
            //return *ui[key];
            break;
        case 1:
            if(debug>1000) std::cout<<"int "<<*in[key]<<std::endl; 
            //return *in[key];
            break;
        case 2:
            if(debug>1000) std::cout<<"float "<<*f[key]<<std::endl; 
            //return *f[key];
            break;
        default:
            if(debug>10) std::cout<<"I don't know type "<<branchInfos[key]->type<<" yet..."<<std::endl;
            //return;
        }
    }
}


//Value AnalysisManager::RetrieveValue(std::string key)
//{
//    //get value
//    //std::string value = *f[key];
//    std::string value("123");
//    return { value };
//}


double AnalysisManager::EvalDeltaR(double eta0, double phi0, double eta1, double phi1)
{
  double dEta = fabs(eta0-eta1);
  double dPhi = fabs(phi0-phi1);

  if(dPhi > 3.14159)
    dPhi = 2.0*3.14159 - dPhi;

  return TMath::Sqrt(TMath::Power(dEta,2)+TMath::Power(dPhi,2));
}
