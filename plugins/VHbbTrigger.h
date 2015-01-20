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

    std::pair<int,int> HighestPtJets(float vetoPhi=-10, float vetoEta=-10);
};

