//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 27 11:04:51 2014 by ROOT version 5.27/06b
// from TTree tree/myTree
// found on file: /eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/WH_125_lumiWeighted.root
//////////////////////////////////////////////////////////

#ifndef AnalysisManager_h
#define AnalysisManager_h

#define PI 3.14159265358979323846
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "HelperClasses/BDTInfo.h"
#include "HelperClasses/BranchInfo.h"
#include "HelperClasses/InfoStructs.h"
#include "HelperClasses/SystematicContainer.cc"
#include "HelperClasses/SampleContainer.cc"
#include "HelperClasses/SFContainer.cc"
#include "HelperClasses/BTagCalibrationStandalone.h"

//#include "HelperClasses/BaseAnalysis.h"

#include <string>
#include <memory>

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
    AnalysisManager();
    virtual ~AnalysisManager();

    void Initialize(std::string filename);

    //File management
    TChain                          *fChain;        //!pointer to the analyzed TTree or TChain
    Int_t                           fCurrent;       //!current Tree number in a TChain
    TFile                           *ofile;         // what will be the trees saved here
    TTree                           *outputTree;    // what will be the condensed output tree
    TTree                           *settingsTree;  // contains analysis settings
    std::string                     outputTreeName;
    std::vector<std::string>        readFiles;       // keep track of which files we've already read
    std::map<std::string,BDTInfo*>  bdtInfos;
    std::map<std::string,BDTInfo*>::iterator  thisBDTInfo;
    std::vector<SampleContainer>    samples;
    SampleContainer*                cursample;
    void                            AddSample(SampleContainer sample);
    std::vector<SystematicContainer>  systematics;
    SystematicContainer*            cursyst;
    void                            AddSystematic(SystematicContainer syst);
    std::vector<SFContainer>        scaleFactors;
    void                            AddScaleFactor(SFContainer sf);
    std::unique_ptr<BTagCalibrationReader> bTagCalibReader;
    void                            InitializeBTagSF(const std::string & bTagCalibFile);
    void                            AddBDT(std::string name, BDTInfo bdtInfo);
    void                            SetJet1EnergyRegression(BDTInfo reg1);
    void                            ConfigureOutputTree();

    //General Physics Information
    float intL;

    int debug;
    bool safemode;


    //Branch management
    std::map<std::string,BranchInfo* > branchInfos;

    std::map<std::string,TBranch*> branches;

    std::map<std::string,unsigned int*> ui;
    std::map<std::string,int*> in;  // if is it "i" then the map is const somehow
    std::map<std::string,float*> f;
    std::map<std::string,double*> d;
    std::map<std::string,bool*> b;
    double m(std::string key,int index=-1); // scaffold for function to get values from maps
    int mInt(std::string key,int index=-1); 

    //Special branches
    // derived variables used in BDT training
    float testMass; //what is this?!?!
    HiggsInfo H;
    TrackInfo V;
    METInfo METtype1corr;

    TBranch* b_H;
    TBranch* b_V;
    TBranch* b_METtype1corr;

    Int_t           GetEntry(Long64_t entry);
    Long64_t        LoadTree(Long64_t entry);
    void            InitChain(std::string filename);

    void            SetupBranch(std::string name, int type, int length=-1, int onlyMC=0, std::string prov="existing", std::string lengthBranch="");
    void            SetupNewBranch(std::string name, int type, int length=-1, bool newmem=true, std::string treetype="output", float val=-999);
    void            SetNewBranches();
    void            ResetBranches();
    void            SetBranches();
    void            PrintBranches();
    void            GetEarlyEntries(Long64_t entry,bool isData=false);
    std::vector<std::string>      ListSampleNames();
    void            PrintBDTInfoValues(BDTInfo* bdt);

    //void            Loop(std::string sampleName="", std::string filename="", int fNum=1 );
    void            Loop(std::string sampleName="", std::string filename="", std::string ofilename="test.root", bool doSkim=false );
    //virtual void     WriteBDTs(std::string indirname, std::string infilename, std::string outdirname, std::string outfilename, std::string cutstring);
    //Value            RetrieveValue(std::string key);


    //AnalysisFuctions
    //BaseAnalysis    *analysis;
    //void            SetAnalysis(BaseAnalysis* scaffold);
    virtual void    InitAnalysis();
    virtual bool    Preselection();
    virtual bool    Analyze();
    virtual void    FinishEvent();
    virtual void    TermAnalysis();

    // general use functions
    double          EvalDeltaR(double eta0, double phi0, double eta1, double phi1);
    double          EvalDeltaPhi(double phi0, double phi1);
    void            SetupBDT(BDTInfo* bdtInfo);
    void            SetBDTVariables(BDTInfo* bdtInfo);
    float           EvaluateMVA(BDTInfo* bdtInfo);
    float           EvaluateRegression(BDTInfo* bdtInfo);
    void            SetupSystematicsBranches();
    void            ApplySystematics(bool early=false);

};

#endif
