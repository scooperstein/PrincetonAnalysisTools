//
//  Written by Chris Palmer 
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include <iostream>
#include "../AnalysisManager.h"

using namespace std;


class VHbbTrigger : public AnalysisManager {
    
    public:
    VHbbTrigger();
    virtual ~VHbbTrigger();

    virtual void InitAnalysis();
    virtual bool Preselection();
    virtual bool Analyze();
    virtual void FinishEvent();
    virtual void TermAnalysis();

    std::pair<int,int> HighestPtJets(float jetPtCut=30, float jetEtaCut=2.5, float vetoPhi=-10, float vetoEta=-10);
    std::pair<int,int> HighestPtPFJets(float jetPtCut=20, float jetEtaCut=3.5, float vetoPhi=-10, float vetoEt=-10);

    bool PassEGL1(float ptCut=40, bool isolated=false, bool etaRestricted=false, int*  L1ElInd=0, float* L1ElPhi=0);
    bool PassEleHLT(float ptCut, std::string WP, std::vector<std::string> L1seeds, bool etaRestricted=false, int* elInd=0);
//
    bool PassEGPlusJetL1(float EGPtCut=20, bool EGiso=true, bool EGer=false, float jetPtCut=20, float jetEtaCut=3.5);
    bool PassElePlusJetHLT(std::vector<std::string> L1seeds, float elPtCut, std::string elWP, bool elEtaRestricted, float jetPtCut, float jetEtaCut, int* elInd=0, int* jetInd=0);
    void RecomputeWPs();
};

