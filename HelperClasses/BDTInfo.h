#ifndef BDTInfo_h
#define BDTInfo_h

#include <map>
#include <string>

class BDTInfo{

public:
  std::map<char*,float*> inputs;
  std::map<char*,float*> secInputs;
  std::string method;
  std::string xmlFile;
  std::string outputBranch;

  BDTInfo(std::map<char*,float*>, std::map<char*,float*>, char*, char*, char*);
  BDTInfo();

};

BDTInfo::BDTInfo(std::map<char*,float*> _inputs, std::map<char*,float*> _secInputs, char* _method, char* _xmlFile, char* _outputBranch) {
    inputs = _inputs;
    secInputs = _secInputs;
    method = _method;
    xmlFile = _xmlFile;
    outputBranch = _outputBranch;

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

BDTInfo::BDTInfo() {
    inputs = std::map<char*,float*>();
    secInputs = std::map<char*,float*>();
    method = "";
    xmlFile = "";
    outputBranch = "";
}


#endif
