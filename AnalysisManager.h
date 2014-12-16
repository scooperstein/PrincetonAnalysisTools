//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 27 11:04:51 2014 by ROOT version 5.27/06b
// from TTree tree/myTree
// found on file: /eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/WH_125_lumiWeighted.root
//////////////////////////////////////////////////////////

#ifndef AnalysisManager_h
#define AnalysisManager_h

#define LEPTONMAX 50
#define JETMAX 200


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "HelperClasses/BDTInfo.h"
#include "HelperClasses/BranchInfo.h"
#include "HelperClasses/infoStructs.h"
#include "HelperClasses/SampleContainer.h"
#include "HelperClasses/SampleContainer.cc"

#include <string>

struct Value
{
    std::string _value;

    template<typename T>
    operator T() const   //implicitly convert into T
    {
       std::stringstream ss(_value);
       T convertedValue;
       if ( ss >> convertedValue ) return convertedValue;
       else throw std::runtime_error("conversion failed");
    }
};



class AnalysisManager {
public :
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    TTree           *outputTree; // what will be the condensed output tree
    std::string     outputTreeName;
    TMVA::Reader    *thereader; // for evaluating the BDT
    BDTInfo         bdtInfo;
    std::vector<SampleContainer> samples; 

    int debug;
    bool safemode;

    std::map<std::string,BranchInfo* > branchInfos;
    
    std::map<std::string,TBranch*> branches;

    std::map<std::string,unsigned int*> ui;
    std::map<std::string,int*> in;  // if is it "i" then the map is const somehow
    std::map<std::string,float*> f;
    std::map<std::string,double*> d;
    std::map<std::string,bool*> b;

    // derived variables used in BDT training
    float testMass; //what is this?!?!
    HiggsInfo H;
    TrackInfo V;
    METInfo METtype1corr;
   
    TBranch* b_H;
    TBranch* b_V;
    TBranch* b_METtype1corr;
 

    AnalysisManager(const char* fname);
    virtual          ~AnalysisManager();
    double           EvalDeltaR(double eta0, double phi0, double eta1, double phi1);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TChain *tree);
    virtual void     SetupBranch(std::string name, int type, int length=-1);
    virtual void     SetupNewBranch(std::string name, int type, int length=-1, bool newmem=true);
    virtual void     SetNewBranches();
    virtual void     ResetBranches();
    virtual void     SetBranches();
    virtual void     PrintBranches();
    virtual void     Loop();
    //virtual void     WriteBDTs(std::string indirname, std::string infilename, std::string outdirname, std::string outfilename, std::string cutstring);
    virtual void     setupBDT();
    Value            RetrieveValue(std::string key);
    void             m(std::string key);
    bool             SampleWenuAnalysis();
    void             AddSample(SampleContainer sample);
    void             ConfigureOutputTree();
};

#endif
