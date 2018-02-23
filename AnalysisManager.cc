#define AnalysisManager_cxx
#include "AnalysisManager.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TMath.h>

AnalysisManager::AnalysisManager(){
    intL=20000; // pb^-1
    settingsTree = new TTree("settings","settings");
    debug = 0;
    BDTisSet = false;
    jet1EnergyRegressionIsSet = false;
    jet2EnergyRegressionIsSet = false;
}

void AnalysisManager::Initialize(std::string filename) {
// used to generate this class and read the Tree.
    InitChain(filename);

    ui.clear();
    in.clear();
    f.clear();
    d.clear();
    b.clear();
    
    branches.clear();
    branchInfos.clear();

    debug=0;
    if(outputTreeName==""){
        outputTreeName="condensed_tree";
    }
    safemode=1;
    
    return;
}


AnalysisManager::~AnalysisManager()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();

    if(debug>10000) std::cout<<"uints"<<std::endl;
    for(std::map<std::string,unsigned int*>::iterator uiit=ui.begin();
            uiit!=ui.end();  ++uiit){
        if(debug>1000) std::cout<<"I'm deleting "<<uiit->first<<std::endl;
        delete uiit->second;
    }

    if(debug>10000) std::cout<<"ints"<<std::endl;
    for(std::map<std::string,int*>::iterator iit=in.begin();
            iit!=in.end();  ++iit){
        if(debug>1000) std::cout<<"I'm deleting "<<iit->first<<std::endl;
        delete iit->second;
    }

    if(debug>10000) std::cout<<"floats"<<std::endl;
    for(std::map<std::string,float*>::iterator fit=f.begin();
            fit!=f.end();  ++fit){
        if(debug>1000) std::cout<<"I'm deleting "<<fit->first<<std::endl;
        delete fit->second;
    }
    
    if(debug>10000) std::cout<<"doubles"<<std::endl;
    for(std::map<std::string,double*>::iterator dit=d.begin();
            dit!=d.end();  ++dit){
        if(debug>1000) std::cout<<"I'm deleting "<<dit->first<<std::endl;
        delete dit->second;
    }

    if(debug>10000) std::cout<<"bools"<<std::endl;
    for(std::map<std::string,bool*>::iterator bit=b.begin();
            bit!=b.end();  ++bit){
        if(debug>1000) std::cout<<"I'm deleting "<<bit->first<<std::endl;
        delete bit->second;
    }
    
    //if(debug>1000) std::cout<<"deleting settingsTree"<<std::endl;
    //delete settingsTree;
}


void AnalysisManager::AddSample(SampleContainer sample){
    samples.push_back(sample);
}

void AnalysisManager::AddSystematic(SystematicContainer syst){
    systematics.push_back(syst);
}

void AnalysisManager::AddScaleFactor(SFContainer sf) {
    scaleFactors.push_back(sf);
    std::cout<<"added sf: "<<sf.name<<", now going to test it"<<std::endl;
    float sf_err = 1.0;
    float sfac = sf.getScaleFactor(100.,1.2,sf_err);
    std::cout<<sfac<<std::endl;
}

void AnalysisManager::AddBDT(BDTInfo bdt) {
    BDTisSet = true;
    bdtInfos.push_back(bdt);
    SetupBDT(bdt);
}

void AnalysisManager::SetJet1EnergyRegression(BDTInfo reg1) {
    jet1EnergyRegression = reg1;
    jet1EnergyRegressionIsSet = true;
    if(debug>10000) {
        PrintBDTInfoValues(reg1);
    }
    SetupBDT(reg1);
}
void AnalysisManager::SetJet2EnergyRegression(BDTInfo reg2) {
    jet2EnergyRegression = reg2;
    jet2EnergyRegressionIsSet = true;
    if(debug>10000) {
        PrintBDTInfoValues(reg2);
    }
    SetupBDT(reg2);
}

void AnalysisManager::PrintBDTInfoValues(BDTInfo bdt) {
    std::cout<<"Printing information for BDT "<<bdt.bdtname<<"..."<<std::endl;
    /*for (unsigned int i=0; i<bdt.inputNames.size(); i++) {
         std::cout<<"Input variable: "<<bdt.inputNames[i].c_str()<<", reference in tree: "<<bdt.localVarNames[i].c_str()<<", current value: "<<*f[bdt.localVarNames[i]]<<std::endl;
     }
     for (unsigned int i=0; i<bdt.inputSpectatorNames.size(); i++) {
         std::cout<<"Spectator variable: "<<bdt.inputSpectatorNames[i].c_str()<<", reference in tree: "<<bdt.localSpectatorVarNames[i].c_str()<<", current value: "<<*f[bdt.localSpectatorVarNames[i]]<<std::endl;
     }*/
     for (unsigned int i=0; i < bdt.bdtVars.size(); i++) {
         BDTVariable bdtvar = bdt.bdtVars[i];
         std::cout<<"Input variable: "<<bdtvar.varName.c_str()<<", reference in tree: "<<bdtvar.localVarName.c_str()<<", current value: "<<*f[bdtvar.localVarName]<<", isSpec: "<<bdtvar.isSpec<<std::endl;
     }
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


void AnalysisManager::InitChain(std::string filename)
{
    // The InitChain() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    //fChain = new TChain("tree");
    fChain = new TChain("Events");
    //std::cout<<"opening "<<filename.c_str()<<std::endl;
    //TFile* tf = TFile::Open(filename.c_str());
    //std::cout<<"adding to chain"<<std::endl;
    //fChain->Add(tf);
    fChain->Add(filename.c_str());
    fCurrent = -1;
    fChain->SetMakeClass(1);

    if(debug>100) std::cout<<"Resetting branches"<<std::endl;
    ResetBranches();

    // Set special branches
    //fChain->SetBranchAddress("H", &H, &b_H);
    //fChain->SetBranchAddress("V", &V, &b_V);
    //fChain->SetBranchAddress("METtype1corr", &METtype1corr, &b_METtype1corr);
}


void AnalysisManager::SetupBranch(std::string name, int type, int length, int onlyMC, std::string prov, std::string lengthBranch){
    branches[name] = new TBranch;
    if(debug>10000) std::cout<<"new TBranch"<<std::endl;
    branchInfos[name] = new BranchInfo(name,type,length,onlyMC,prov,-999.,lengthBranch); 
    if(debug>10000) std::cout<<"new BranchInfo"<<std::endl;

    // Only 0-9 are setup with types for the moment.
    if(type>9 || type<0) {
        std::cout<<"Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }
    if(debug>10000) std::cout<<"checking type and setting branch address.  type "<<type<<std::endl;

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
            <<length<<" name "<<name.c_str()<<std::endl;
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
    if(debug>100000) std::cout<<"getting entries of outputTree"<<std::endl;
    if(debug>100000) std::cout<<"entries "<<outputTree->GetEntries()<<std::endl;
}

void AnalysisManager::SetupNewBranch(std::string name, int type, int length, bool newmem, std::string treetype, float val){
    if(debug>1000) {
        std::cout<<"SetupNewBranch "<<name<<std::endl;
        std::cout<<"treetype name type length val \t"<<treetype<<" "<<name<<" "<<type<<" "<<length<<" "<<val<<std::endl;
    }
    
    TTree* treeptr;
    if(treetype=="output") { // outputtree
        treeptr=outputTree;
    } else if(treetype=="settings") { // settingstree
        treeptr=settingsTree;
    } else {
        std::cout<<"treetype "<<treetype<<" is unknown.  Not setting up "<<name<<std::endl;
        return;
    }
    
    if(newmem) {
        if(treetype=="output") { // outputtree
            branchInfos[name] = new BranchInfo(name,type,length,false,"new");
        } else if(treetype=="settings") { // settingstree
            branchInfos[name] = new BranchInfo(name,type,length,false,"settings",val);
        }
    }
    if(debug>1000) std::cout<<"BranchInfo instaniated"<<std::endl; 
       

    if(type>9 || type<0) {
        std::cout<<"New Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }

    if(type==0) {
        if(newmem) ui[name] = new unsigned int;
        branches[name] = treeptr->Branch(name.c_str(), ui[name], Form("%s/i",name.c_str()));
    } else if(type==1) {
        if(newmem) in[name] = new int;
        branches[name] = treeptr->Branch(name.c_str(), in[name], Form("%s/I",name.c_str()));
    } else if(type==2) {
        if(newmem) f[name] = new float;
        branches[name] = treeptr->Branch(name.c_str(), f[name], Form("%s/F",name.c_str()));
    } else if(type==3) {
        if(newmem) d[name] = new double;
        branches[name] = treeptr->Branch(name.c_str(), d[name], Form("%s/D",name.c_str()));
    } else if(type==4) {
        if(newmem) b[name] = new bool;
        branches[name] = treeptr->Branch(name.c_str(), b[name], Form("%s/O",name.c_str()));
    }


    if(type>4 && type<10 && length<0) {
        std::cout<<"Types 5-9 are arrays and need a length greater than 0... not "
            <<length<<" name "<<name.c_str()<<std::endl;
        return;
    }

    if(type==5) {
        if(newmem) ui[name] = new unsigned int[length];
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), ui[name], Form("%s[%i]/i",name.c_str(),length));
    } else if(type==6) {
        if(newmem) in[name] = new int[length];
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), in[name], Form("%s[%i]/I",name.c_str(),length));
    } else if(type==7) {
        if(newmem) f[name] = new float[length];
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), f[name], Form("%s[%i]/F",name.c_str(),length));
    } else if(type==8) {
        if(newmem) d[name] = new double[length];
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), d[name], Form("%s[%i]/D",name.c_str(),length));
    } else if(type==9) {
        if(newmem) b[name] = new bool[length];
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), b[name], Form("%s[%i]/O",name.c_str(),length));
    }

    return;
}

void AnalysisManager::ResetBranches(){
    
    for(std::map<std::string,BranchInfo*>::iterator biit=branchInfos.begin();
            biit!=branchInfos.end(); ++biit) {
        if(debug>100) std::cout<<"Branch "<<biit->second->name<<" of type, prov "<<biit->second->type<<" "<<biit->second->prov<<std::endl;
        if(biit->second->prov == "existing" || biit->second->prov == "early"){
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
    std::cout<<"Branches in branch map"<<std::endl;
    for(std::map<std::string,TBranch*>::iterator ibranch=branches.begin(); 
            ibranch!=branches.end(); ++ibranch){
        std::cout<<ibranch->first<<" "<<branchInfos[ibranch->first]->prov<<std::endl;
    }
}


void AnalysisManager::SetBranches(){
    // only copy branches specified
    if(debug>10) std::cout<<"SetBranchStatus of existing branches"<<std::endl;
    fChain->SetBranchStatus("*", 0);
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin(); 
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "existing" || ibranch->second->prov == "early") {
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
        if(ibranch->second->prov == "settings") {
            if(debug>100) std::cout<<"Setting value "<<ibranch->first.c_str()<<std::endl;
            // Branch is already setup properly
            //SetupNewBranch(ibranch->first, ibranch->second->type, ibranch->second->length, false, "settings", ibranch->second->val); //newmem is false
            // Set the value to val which is read from settings.txt
            *f[ibranch->second->name]=ibranch->second->val;
            if(debug>100) std::cout<<"Setting value done "<<ibranch->first.c_str()<<std::endl;
        }
    }
}

void AnalysisManager::GetEarlyEntries(Long64_t entry, bool isData){
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin(); 
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "early" && !(isData && ibranch->second->onlyMC)) {
            if(debug>100000) std::cout<<"Getting entry for early branch "<<ibranch->first<<std::endl;
            branches[ibranch->first]->GetEntry(entry);
            if(debug>100000) std::cout<<"Got entry for early branch "<<ibranch->first<<std::endl;
        }
    }
}

// Return a std::vector of the names of the samples
std::vector<std::string> AnalysisManager::ListSampleNames() {
    std::vector<std::string> snlist = std::vector<std::string>();
    for(int i=0; i < (int)samples.size(); i++) {
        snlist.push_back(samples[i].sampleName);
    }
    return snlist;
}


// Process all input samples and all events
void AnalysisManager::Loop(std::string sampleName, std::string filename, std::string ofilename, bool doSkim){
    // Specify sample name if we want to run on only a particular sample, specify
    // filenames if we want to run only on specific files from that sample.
    std::vector<std::string> filenames;
    if (!filename.empty()) {
        std::stringstream ss(filename);
        std::string file;
        while (std::getline(ss, file, ','))
        {
            filenames.push_back(file);
        }
    }
    if(!sampleName.empty()){
        bool sampleFound=false;
        for(int i=0; i<(int)samples.size(); i++) {
            if(sampleName == samples[i].sampleName) {
                ofile = new TFile(ofilename.c_str(), "recreate"); // moved here by Jan
                sampleFound=true;
                SampleContainer* onlySample = new SampleContainer(samples[i]);
                if (filenames.size() > 0) {
                    // keep track of total number of processed events for the sample
                    // so we still get the weight right
                    //int processedEvents = onlySample->processedEvents;
                    std::vector<std::string> sampleFiles = onlySample->files;
                    onlySample->files.clear();
                    //onlySample->sampleChain = new TChain("tree");
                    onlySample->sampleChain = new TChain("Events");
                    for (int j=0; j<(int)filenames.size(); j++) {
                        if (std::find(sampleFiles.begin(), sampleFiles.end(), filenames[j]) != sampleFiles.end() ) {
                            onlySample->AddFile(filenames[j].c_str(),1);
                        }
                        else {
                            std::cout<<"Analysis Manager tried to run on file "<<filenames[j]<<" in sample "<<sampleName<<", but this file is not in the sample's list of files. Skipping..."<<std::endl;
                            std::cout<<"Let's print the full list of files for this sample..."<<std::endl;
                            for (int k=0; k<(int)sampleFiles.size(); k++) {
                                std::cout<<sampleFiles[k]<<std::endl;
                            }
                        }
                    }
                    //onlySample->processedEvents = processedEvents;
                    //ofile = new TFile(Form("%s_%s_%i.root",outputTreeName.c_str(),samples[0].sampleName.c_str(),fNum),"recreate");
                }
                else {
                    //ofile = new TFile(Form("%s_%s.root",outputTreeName.c_str(),samples[0].sampleName.c_str()),"recreate");
                }
                //ofile = new TFile(ofilename.c_str(), "recreate");
                samples.clear();
                samples.push_back(*onlySample);
                break;
            } 
        }
        if(!sampleFound) {
            std::cout<<"Analysis Manager tried to loop over the individual sample "<<sampleName<<", but this sample is not in the samples list. Skipping..."<<std::endl;
        }
        if(debug>10) std::cout<<"Starting Loop over individual sample %s"<<sampleName<<std::endl;
    }
    else {
        ofile = new TFile(Form("%s.root",outputTreeName.c_str()),"recreate");
    }

    if(debug>10) std::cout<<"Starting Loop"<<std::endl;
    SetBranches();
    
    // add new branches
    ofile->cd();
    //TTree *newtree = fChain->CloneTree(0);
    // for now we will use a default file to set the structure for the output tree
    outputTree = new TTree(); 
    outputTree = fChain->CloneTree(0); // need this one to only keep the branches you want 
    
    // it would probably be smarter to check also if SetupBranch() has been called since the last time outputTree was initialized
    //if (!outputTree) {
    //    std::cout<<"Did not call ConfigureOutputTree() before Loop(). Assuming there are no new branches..."<<std::endl;
    //    outputTree = fChain->CloneTree(0);
    //}  
    // let's add some of our own branches
    SetNewBranches();
    SetupSystematicsBranches();
    // FIXME add branches to settings regarding splitting
    //SetupNewBranch("jobNum", 3, -1, true, "settings", jobNum);
    settingsTree->Fill();
    settingsTree->Write();
    delete settingsTree;
 
    
    if(debug>10) std::cout<<"Done setting up branches; about to Init"<<std::endl;
    InitAnalysis();
    
    
    if(debug>10) std::cout<<"About to loop over samples"<<std::endl;
    // loop through one sample at a time
    for (int i = 0; i < (int)samples.size(); i++) { 
        cursample = &samples[i];
        cursample->ComputeWeight(*f["intL"]);
         
        for(int ifile=0; ifile<(int)(cursample->files.size()); ifile++){
            if(std::find(readFiles.begin(), readFiles.end(), cursample->files[ifile]) != readFiles.end() ) {
                // file has already been processed, skip it
                continue;
            }
            readFiles.push_back(cursample->files[ifile]);
            // set fChain to the TChain for the current sample
            InitChain(cursample->files[ifile]);
    
            // FIXME should have a sample name, but doesn't right now
            if(debug>0) std::cout<<"About to loop over events in "<<cursample->files[ifile]<<std::endl;
            // loop through the events
            Long64_t nentries = round(fChain->GetEntries()*cursample->procEff);
            if(debug>1) std::cout<<"looping over "<<nentries<<std::endl;
            Long64_t nbytes = 0, nb = 0;
            int saved=0;
            // FIXME need a loop over systematics
            for (Long64_t jentry=0; jentry<nentries;jentry++) {
                if((jentry%1000==0 && debug>0) || debug>100000)  std::cout<<"entry saved weighted "<<jentry<<" "<<saved<<" "<<saved*cursample->intWeight<<std::endl;
                //if((jentry%10000==0 && debug>0) || debug>100000)  std::cout<<"entry saved weighted "<<jentry<<" "<<saved<<" "<<saved*cursample->intWeight<<std::endl;
                
                
                GetEarlyEntries(jentry, cursample->sampleNum==0);
                bool anyPassing=false;
                for(int iSyst=0; iSyst<systematics.size(); iSyst++){
                    cursyst=&(systematics[iSyst]);
                    if (cursample->sampleNum == 0 && cursyst->name != "nominal") continue;
                    //ApplySystematics(true);

                    if(debug>100000) std::cout<<"checking preselection"<<std::endl;
                    bool presel = Preselection();
                    if(presel) anyPassing=true;
                }
                
                if(!anyPassing) continue;
  

                if(debug>1000) std::cout<<"passed presel; Loading tree"<<std::endl;
                Long64_t ientry = LoadTree(jentry);
                if (ientry < 0) break;
                //nb = fChain->GetEntry(jentry);   
                //nbytes += nb;
                
                anyPassing=false;
                for(int iSyst=0; iSyst<systematics.size(); iSyst++){
                    nb = fChain->GetEntry(jentry);   
                    nbytes += nb;
                    
                    if (!doSkim) {

                        cursyst=&(systematics[iSyst]);
                        if (cursample->sampleNum == 0 && cursyst->name != "nominal") continue;
                        ApplySystematics();
                        if(debug>1000) std::cout<<"running analysis"<<std::endl;
                        bool select = Analyze();
                        *b[Form("Pass_%s",cursyst->name.c_str())] = false;
                        if(select) { 
                            anyPassing=true;
                            *b[Form("Pass_%s",cursyst->name.c_str())] = true;
                        }
                        if(select || (cursyst->name=="nominal" && anyPassing)){
                            if(debug>1000) std::cout<<"selected event; Finishing"<<std::endl;
                            //for (int i=0; i < scaleFactors.size(); i++) {
                            //    SFContainer sf = scaleFactors[i];
                            //    float sf_err = 0.0;
                            //    if (cursample->sampleNum != 0) {
                            //        for (int j=0; j < *in[sf.length]; j++) {
                            //            if (sf.binning.find("abs") == -1) {
                            //                f[sf.branchname][j] = sf.getScaleFactor(f[sf.branches[0]][j], f[sf.branches[1]][j], sf_err);
                            //            }
                            //            else {
                            //                 f[sf.branchname][j] = sf.getScaleFactor(fabs(f[sf.branches[0]][j]), fabs(f[sf.branches[1]][j]), sf_err);
                            //            }
                            //            f[Form("%s_err",sf.branchname.c_str())][j] = sf_err; 
                            //        }
                            //    }
                            //    else {
                            //        // data event, scale factor should just be 1.0
                            //        for (int j=0; j < *in[sf.length]; j++) {
                            //            f[sf.branchname][j] = 1.0;
                            //            f[Form("%s_err",sf.branchname.c_str())][j] = 0.0; 
                            //        }
                            //    }
                            //    
                            //}
                            FinishEvent();
                            if(cursyst->name=="nominal") saved++;
                        }
                    }
                    else {
                        // running skim
                        ofile->cd();
                        outputTree->Fill();
                        saved++;
                    }
                }
            } // end event loop
        } // end file loop
        ofile->cd();
        if (!doSkim) {
            //cursample->CountWeightedLHEWeightScale->Write(Form("CountWeightedLHEWeightScale_%s",cursample->sampleName.c_str()));
            //cursample->CountWeightedLHEWeightPdf->Write(Form("CountWeightedLHEWeightPdf_%s",cursample->sampleName.c_str()));
            cursample->CountWeighted->Write(Form("CountWeighted_%s",cursample->sampleName.c_str()));
            //cursample->CountFullWeighted->Write(Form("CountFullWeighted_%s",cursample->sampleName.c_str()));
        }
        else {
            //cursample->CountWeightedLHEWeightScale->Write();
            //cursample->CountWeightedLHEWeightPdf->Write();
            cursample->CountWeighted->Write();
            //cursample->CountFullWeighted->Write();
        }
    } // end sample loop
    if(debug>1000) std::cout<<"Finished looping"<<std::endl;
    
    
    TermAnalysis();
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
  
    // FIXME nominal must be last!
    if(cursyst->name=="nominal"){
        ofile->cd();
        outputTree->Fill();
    }
    return;
}


void AnalysisManager::TermAnalysis() {
    // save tree here
    ofile->cd();
    outputTree->Write();
    ofile->Close();
}


// Set up all the BDT branches and configure the BDT's with the same input variables as used in training. Run before looping over events.
void AnalysisManager::SetupBDT(BDTInfo bdtInfo) {
    //TMVA::Reader *thereader = bdtInfo.reader;
    
    for (unsigned int i=0; i < bdtInfo.bdtVars.size(); i++) {
        BDTVariable bdtvar = bdtInfo.bdtVars[i];
        if (!bdtvar.isSpec) {
            bdtInfo.reader->AddVariable(bdtvar.varName, f[bdtvar.localVarName]);
        } else {
            bdtInfo.reader->AddSpectator(bdtvar.varName, f[bdtvar.localVarName]);
        }
    }

    std::cout<<"booking MVA for bdt with name...  "<<bdtInfo.bdtname<<std::endl;
    bdtInfo.reader->BookMVA(bdtInfo.bdtmethod, bdtInfo.xmlFile);
}


void AnalysisManager::SetupSystematicsBranches(){
    std::cout<<"loop through systematics "<<systematics.size()<<std::endl;
    for(int iSyst=0; iSyst<systematics.size(); iSyst++){
        std::cout<<"syst name "<<systematics[iSyst].name<<std::endl;
        for(int iBrnch=0; iBrnch<systematics[iSyst].branchesToEdit.size(); iBrnch++){
            std::string newName(systematics[iSyst].branchesToEdit[iBrnch]);
            newName.append("_");
            newName.append(systematics[iSyst].name);
            std::cout<<"\tbranch name "<<newName<<std::endl;
            std::cout<<"\t\ttype, length,  "
                <<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->type<<", "
                <<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length<<std::endl;
                //<<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length<<" "
            SetupNewBranch(newName, 
                branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->type,
                branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length);
        }
        for(int iBDT=0; iBDT<bdtInfos.size(); iBDT++) {
            std::string bdtname(bdtInfos[iBDT].bdtname);
            if (systematics[iSyst].name != "nominal") {
                bdtname.append("_");
                bdtname.append(systematics[iSyst].name);
                SetupNewBranch(bdtname, 2);
            }
        }
        if (systematics[iSyst].name != "nominal") {
            SetupNewBranch(Form("H_mass_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("H_pt_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("Jet_btagCSV_%s", systematics[iSyst].name.c_str()), 7, 100);
            SetupNewBranch(Form("nAddJets252p9_puid_%s", systematics[iSyst].name.c_str()), 1);
        }
    }
}

// applies scales to floats and doubles (and arrays)
// smearing can be added; other types can be added.
void AnalysisManager::ApplySystematics(bool early){
    for(int iBrnch=0; iBrnch<cursyst->branchesToEdit.size(); iBrnch++){
        //std::cout<<"iBrnch "<<iBrnch<<std::endl;
        std::string oldBranchName(cursyst->branchesToEdit[iBrnch]);
        BranchInfo* existingBranchInfo = branchInfos[oldBranchName.c_str()];
        //std::cout<<"early? "<< branchInfos[cursyst->branchesToEdit[iBrnch]]->prov<<std::endl;
        if(!early || (early && branchInfos[cursyst->branchesToEdit[iBrnch]]->prov == "early")){
            std::string systBranchName(oldBranchName);
            systBranchName.append("_");
            systBranchName.append(cursyst->name);
            //std::cout<<"new branch "<<systBranchName.c_str()<<std::endl;
            int thisType=existingBranchInfo->type;
            // just doing floats and doubles
            // smearing can be added at any time
            float jetPtSplit = 100.; // classify high/low jets by this pT threshold
            float jetEtaSplit = 1.4; // classify central/forward jets by this eta threshold
            int ptSplit = 0; // if 0 don't split, if 1 split low, if 2 splight high
            int etaSplit = 0; // if 0 don't split, if 1 split low, if 2 splight high
            if (cursyst->name.find("High")!=std::string::npos) { ptSplit = 2; }
            if (cursyst->name.find("Low")!=std::string::npos) { ptSplit = 1; }
            if (cursyst->name.find("Forward")!=std::string::npos) { etaSplit = 2; }
            if (cursyst->name.find("Central")!=std::string::npos) { etaSplit = 1; }
            //std::cout<<cursyst->name<<std::endl;
            //std::cout<<ptSplit<<", "<<etaSplit<<std::endl;
            if(thisType==2){
                // scale the current branch
                if (cursyst->scaleVar[iBrnch] == "") {
                    // flat scaling
                    *f[oldBranchName.c_str()]=*f[oldBranchName.c_str()] * cursyst->scales[iBrnch];
                }
                else {
                    // dynamic scaling
                    *f[oldBranchName.c_str()]=*f[oldBranchName.c_str()] * *f[cursyst->scaleVar[iBrnch]];
                }
                // copy the value to the new branch
                *f[systBranchName.c_str()]=*f[oldBranchName.c_str()];
            } else if(thisType==3){
                // scale the current branch
                if (cursyst->scaleVar[iBrnch] == "") {
                    // flat scaling
                    *d[oldBranchName.c_str()]=*d[oldBranchName.c_str()] * cursyst->scales[iBrnch];
                }
                else {
                    // dynamic scaling
                    *d[oldBranchName.c_str()]=*d[oldBranchName.c_str()] * *d[cursyst->scaleVar[iBrnch]];
                }
                // copy the value to the new branch
                *d[systBranchName.c_str()]=*d[oldBranchName.c_str()];
            } else if(thisType==7){
                //std::cout<<"length branch "<<existingBranchInfo->lengthBranch<<std::endl;
                //std::cout<<"length "<<*in[existingBranchInfo->lengthBranch]<<std::endl;
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){// scale the current branch
                    //std::cout<<"Jet pt, eta: "<<f["Jet_pt_reg"][ind]<<", "<<f["Jet_eta"][ind]<<std::endl;
                    //std::cout<<"old val "<<f[oldBranchName.c_str()][ind]<<std::endl;
                    // FIXME: we are assuming that the only array variables we scale are jet variables, we should in principle make this more generic
                    if (ptSplit==0 || (ptSplit==1 && f["Jet_pt_reg"][ind]<jetPtSplit) || (ptSplit==2 && f["Jet_pt_reg"][ind]>jetPtSplit)) { 
                        if (etaSplit==0 || (etaSplit==1 && fabs(f["Jet_eta"][ind])<jetEtaSplit) || (etaSplit==2 && fabs(f["Jet_eta"][ind])>jetEtaSplit)) { 
                            //std::cout<<"got through"<<std::endl;
                            if (cursyst->scaleVar[iBrnch] == "") {
                                // flat scaling
                                f[oldBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind] * cursyst->scales[iBrnch];
                            }
                            else {
                                //std::cout<<f[oldBranchName.c_str()][ind]<<" * "<<f[cursyst->scaleVar[iBrnch]][ind]<<std::endl;
                                // dynamic scaling
                                f[oldBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind] * f[cursyst->scaleVar[iBrnch]][ind];
                            }
                            //std::cout<<"new val "<<f[oldBranchName.c_str()][ind]<<std::endl;
                            // copy the value to the new branch
                            f[systBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind];
                            //std::cout<<"new branch "<<f[systBranchName.c_str()][ind]<<std::endl;
                        }
                    }
                }
            } else if(thisType==8){
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){// scale the current branch
                    if (cursyst->scaleVar[iBrnch] == "") {
                         // flat scaling
                        d[oldBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind] * cursyst->scales[iBrnch];
                    }
                    else {
                        // dynamic scaling
                        d[oldBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind] * d[cursyst->scaleVar[iBrnch]][ind];
                    }
                    // copy the value to the new branch
                    d[systBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind];
                }
            }
        }
    }
}


double AnalysisManager::m(std::string key, int index){
    if(debug>1000) std::cout<<"looking for key "<<key<<" with index "<<index<<std::endl;
    if(branchInfos.count(key)==0){
        std::cout<<"There is no branch with name "<<key<<std::endl;
        if(safemode){
            std::cout<<"The program must be terminated..."<<std::endl;
            std::exit(0);
        } else {
            if(debug>1) std::cout<<"Returning -999 and hoping for the best."<<std::endl;
            return -999;
        }
    } else {
        if(debug>1000) std::cout<<"here is where I should return the right branch value"<<std::endl;
        switch(branchInfos[key]->type)
        {
        case 0:
            if(debug>1000) std::cout<<"unsigned int "<<*ui[key]<<std::endl; 
            return (double)*ui[key];
        case 1:
            if(debug>1000) std::cout<<"int "<<*in[key]<<std::endl; 
            return (double)*in[key];
        case 2:
            if(debug>1000) std::cout<<"float "<<*f[key]<<std::endl; 
            return (double)*f[key];
        case 3:
            if(debug>1000) std::cout<<"double "<<*d[key]<<std::endl; 
            return (double)*d[key];
        case 4:
            if(debug>1000) std::cout<<"boolean "<<*b[key]<<std::endl; 
            return (double)*b[key];
        case 5:
            if(debug>1000) std::cout<<"unsigned int "<<ui[key][index]<<std::endl; 
            return (double)ui[key][index];
        case 6:
            if(debug>1000) std::cout<<"int "<<in[key][index]<<std::endl; 
            return (double)in[key][index];
        case 7:
            if(debug>1000) std::cout<<"float "<<f[key][index]<<std::endl; 
            return (double)f[key][index];
        case 8:
            if(debug>1000) std::cout<<"double "<<d[key][index]<<std::endl; 
            return (double)d[key][index];
        case 9:
            if(debug>1000) std::cout<<"boolean "<<b[key][index]<<std::endl; 
            return (double)b[key][index];
        default:
            if(debug>10) std::cout<<"I don't know type "<<branchInfos[key]->type<<" yet..."<<std::endl;
            return -999;
        }
    }
}


double AnalysisManager::EvalDeltaPhi(double phi0, double phi1){
    double dPhi = fabs(phi0-phi1);
    //std::cout<<"dPhi PI "<<dPhi<<" "<<PI<<std::endl;

    if(dPhi > PI)
        dPhi = 2.0*PI - dPhi;
    
    return dPhi;
}


double AnalysisManager::EvalDeltaR(double eta0, double phi0, double eta1, double phi1)
{
  double dEta = fabs(eta0-eta1);
  double dPhi = EvalDeltaPhi(phi0, phi1);


  return TMath::Sqrt(TMath::Power(dEta,2)+TMath::Power(dPhi,2));
}
