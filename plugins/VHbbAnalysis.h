//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include <iostream>
#include "../AnalysisManager.h"
#include "TLorentzVector.h"

using namespace std;


class VHbbAnalysis : public AnalysisManager {
    
    public:
    VHbbAnalysis();
    virtual ~VHbbAnalysis();

    virtual void InitAnalysis();
    virtual bool Preselection();
    virtual bool Analyze();
    virtual void FinishEvent();
    virtual void TermAnalysis();

    bool ElectronSelection();
    bool MuonSelection();
    float ReWeightMC(int nPU=0);
    std::pair<int,int> HighestPtBJets();
    std::pair<int,int> HighestCSVBJets();
    std::pair<int,int> HighestPtJJBJets();
    double GetRecoTopMass(TLorentzVector Jet, bool isJet=true, int useMET=0, bool regPT=true);
    float ptWeightQCD(int nGenVbosons=0, float lheHT=0., int GenVbosons_pdgId=0);
    float ptWeightEWK(int nGenVbosons=0,float GenVbosons_pt=0.,int VtypeSim=0,int GenVbosons_pdgId=0);
};

