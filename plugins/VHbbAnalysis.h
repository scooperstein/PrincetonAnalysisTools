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

    bool WenuHbbSelection();
    bool WmunuHbbSelection();
    std::pair<int,int> HighestPtBJets();
    std::pair<int,int> HighestCSVBJets();
    double GetRecoTopMass(TLorentzVector Jet, bool isJet=true);
};

