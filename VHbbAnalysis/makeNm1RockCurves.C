void makeNm1RockCurves(const char* fname) {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L ../../../scripts/makeBDTRockCurve.C");

const int nbdts = 5;
char* bdtvars[nbdts] = {"wHVdEta","wHVdEta_4MET","wJJEtaBal_v2","wJJEtaBal_and_HVdEta_4MET","TMBestCSV"};
//const int nbdts = 23;
//char* bdtvars[nbdts] = {"top10Var","noMbb","noHPt","noWPt","noTopMass","noHVdPhi","noNAddJet","noLepMetDPhi","noJet1Pt","noJet2Pt","noLepPt","noLepIso","noJet1CSV","noLepEta","noJet1Eta","noJet2Eta","noPtBal","noAddJetCSV","noWMt","noMET"};
//char* bdtvars[nbdts] = {"top10Var","noMbb","noDRjj","noHPt","noWPt","noJet2CSV","noTopMass","noHVdPhi","noNAddJet","noLepMetDPhi","noNSA5Jet","noJet1Pt","noJet2Pt","noLepPt","noLepIso","noJet1CSV","noLepEta","noJet1Eta","noJet2Eta","noPtBal","noAddJetCSV","noWMt","noMET"};
//const int nbdts = 12;
//char* bdtvars[nbdts] = {"noMbb","noHPt","noWPt","noJet2CSV","noTopMass","noHVdPhi","noNAddJet","noLepMetDPhi","noNSA5Jet","noLepPt","noWMt","noMET"};

for (int i=0; i<nbdts; i++) {
    gROOT->ProcessLine("std::vector<string> bdtnames;");
    gROOT->ProcessLine("bdtnames.push_back(\"BDT_May16_450_3\")");
    gROOT->ProcessLine(Form("bdtnames.push_back(\"BDT_May16_%s\")",bdtvars[i]));
    gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,0,\"%s\")",fname,bdtvars[i]));
    gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,1,\"%s\")",fname,bdtvars[i]));
    gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,2,\"%s\")",fname,bdtvars[i]));
    gROOT->ProcessLine("bdtnames.clear()");
}

}
