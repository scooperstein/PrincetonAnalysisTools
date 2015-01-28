#ifndef BDTInfo_h
#define BDTInfo_h

#include <map>
#include <string>

#include "TMVA/Reader.h"

class BDTInfo{

public:
  std::string bdtname;
  std::vector<std::string> inputNames;
  std::vector<std::string> localVarNames;
  std::vector<std::string> inputSpectatorNames;
  std::vector<std::string> localSpectatorVarNames;
  //std::string method;
  std::string xmlFile;
  TMVA::Reader *reader;

  BDTInfo(std::string, std::string);
  BDTInfo();
  void AddVariable(std::string, std::string);
  void AddSpectatorVariable(std::string, std::string);

};

inline BDTInfo::BDTInfo(std::string _bdtname, std::string _xmlFile) {
    bdtname = _bdtname;
    inputNames = std::vector<std::string>();
    localVarNames = std::vector<std::string>();
    //method = _method;
    xmlFile = _xmlFile;
    reader = new TMVA::Reader( "!Color:Silent" );    

   /* map<char*,float*> floatMap;
  floatMap["H.massCorr"] = &H.massCorr;
  floatMap["H.ptCorr"] = &H.ptCorr;
  floatMap["V.pt"] = &V.pt;
  floatMap["hJet_csvReshapedNew[0]"] = &(hJet_csvReshapedNew[0]);
  floatMap["hJet_csvReshapedNew[1]"] = &(hJet_csvReshapedNew[1]);
  floatMap["HVdPhi"] = &HVdPhi;
  floatMap["H.dEta"] = &H.dEta;
  floatMap["H.dR"] = &H.dR;
  floatMap["abs(deltaPullAngle)"] = &absDeltaPullAngle;
  floatMap["hJet_ptCorr[0]"] = &hJet_ptCorr[0];
  floatMap["hJet_ptCorr[1]"] = &hJet_ptCorr[1];
  floatMap["Vtype"] = &Vtype_f;
  floatMap["Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )"] = &naLeptonsPassingCuts;
  floatMap["Sum$(aJet_pt>20 && abs(aJet_eta)<4.5)"] = &naJetsPassingCuts;
  floatMap["METtype1corr.et"] = &METtype1corr.et;
  floatMap["MET.et"] = &METtype1corr.et;
  floatMap["MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)"] = &highestCSVaJet;
  floatMap["MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)"] = &minDeltaRaJet;
  */
}

inline BDTInfo::BDTInfo() {
    BDTInfo("", "");
}

inline void BDTInfo::AddVariable(std::string varName, std::string localVarName) {
    inputNames.push_back(varName);
    localVarNames.push_back(localVarName);
}

inline void BDTInfo::AddSpectatorVariable(std::string varName, std::string localVarName) {
    inputSpectatorNames.push_back(varName);
    localSpectatorVarNames.push_back(localVarName);
}

#endif
