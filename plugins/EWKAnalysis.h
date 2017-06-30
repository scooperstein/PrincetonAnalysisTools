//  EWK VBF W+2jet analysis
//  Inheriting from Analysis Manager
//

#include <iostream>
#include "../AnalysisManager.h"
#include "TLorentzVector.h"

using namespace std;


class EWKAnalysis : public AnalysisManager {
    
    public:
    EWKAnalysis();
    virtual ~EWKAnalysis();

    virtual void InitAnalysis();
    virtual bool Preselection();
    virtual bool Analyze();
    virtual void FinishEvent();
    virtual void TermAnalysis();

    bool ElectronSelection();
    bool MuonSelection();
    float ptWeightQCD(int nGenVbosons=0, float lheHT=0., int GenVbosons_pdgId=0);
    float ptWeightEWK(int nGenVbosons=0,float GenVbosons_pt=0.,int VtypeSim=0,int GenVbosons_pdgId=0);
    double LOtoNLOWeightBjetSplitEtabb(double etabb=10., int njets=0);
    void  smearJets(float JERScale=1.0);
    float evaluateRegression(int i=0);
    TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET);
    float getQGLCorr(float partonFlav=0, float eta=0., float qgl=-10.); 
    float getQGLCorr_q(float x=-99.);
    float getQGLCorr_g(float x=-99.);
};

