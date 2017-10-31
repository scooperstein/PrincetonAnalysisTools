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

    bool ElectronSelection(int);
    bool MuonSelection(int);
    float ReWeightMC(int nPU=0);
    float puWeight_ichep(int i=0);
    float puWeight_ichep_up(int i=0);
    float puWeight_ichep_down(int i=0);
    std::pair<int,int> HighestPtBJets();
    std::pair<int,int> HighestCSVBJets();
    std::pair<int,int> HighestPtJJBJets();
    double GetRecoTopMass(TLorentzVector Jet, bool isJet=true, int useMET=0, bool regPT=true);
    float ptWeightQCD(int nGenVbosons=0, float lheHT=0., int GenVbosons_pdgId=0);
    float ptWeightEWK(int nGenVbosons=0,float GenVbosons_pt=0.,int VtypeSim=0,int GenVbosons_pdgId=0);
    TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET);
    double LOtoNLOWeightBjetSplitEtabb(double etabb=10., int njets=0);
    float getVPtCorrFactor(float V_pt=0.);
    float getVPtCorrFactorUp(float V_pt=0.);
    float getVPtCorrFactorDown(float V_pt=0.);
    void  smearJets(float JERScale=1.0);
    float evaluateRegression(int i=0);
    void SetupFactorizedJECs(std::string variation = "nominal");
};

