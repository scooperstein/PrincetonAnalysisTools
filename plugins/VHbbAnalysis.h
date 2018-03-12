//
//  Written by Chris Palmer
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#ifndef ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_
#define ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_

#include <iostream>
#include "TLorentzVector.h"
#include "../AnalysisManager.h"


class VHbbAnalysis : public AnalysisManager {

    public:
        VHbbAnalysis();
        virtual ~VHbbAnalysis();

        virtual void InitAnalysis();
        virtual bool Preselection();
        virtual bool Analyze();
        virtual void FinishEvent();
        virtual void TermAnalysis();

	std::string taggerName;
	void SetTaggerName(float taggerType);
		
        std::pair<int,int> HighestPtGoodElectronsOppCharge(float min_pt, float max_rel_iso, float idcut);
        std::pair<int,int> HighestPtGoodMuonsOppCharge(float min_pt, float max_rel_iso);
        bool ElectronSelection(int);
        bool MuonSelection(int);
        int UpdatedVType();
        bool PassVTypeAndTrigger(int vtype);
        float ReWeightMC(int nPU=0);
        float puWeight_ichep(int i=0);
        float puWeight_ichep_up(int i=0);
        float puWeight_ichep_down(int i=0);
        std::pair<int,int> HighestPtBJets();
        std::pair<int,int> HighestTaggerValueBJets(float j1ptCut, float j2ptCut, std::string taggerName);
        std::pair<int,int> HighestDeepCSVBJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestCMVABJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestCSVBJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestPtJJBJets();
        double GetRecoTopMass(TLorentzVector Jet, bool isJet=true, int useMET=0, bool regPT=true);
        float ptWeightQCD(int nGenVbosons=0, float lheHT=0., int GenVbosons_pdgId=0);
        float ptWeightEWK(int nGenVbosons=0, float GenVbosons_pt=0., int VtypeSim=0, int GenVbosons_pdgId=0);
        TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET);
        double LOtoNLOWeightBjetSplitEtabb(double etabb=10., int njets=0);
        float getVPtCorrFactor(float V_pt=0.);
        float getVPtCorrFactorUp(float V_pt=0.);
        float getVPtCorrFactorDown(float V_pt=0.);
        void  smearJets(float JERScale=1.0);
        float evaluateRegression(int i=0);
        void SetupFactorizedJECs(std::string variation="nominal");
};

#endif // ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_
