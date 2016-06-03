void addNm1BDTs(const char* fname) {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L ../../scripts/addBDTToTree.C ");

//const int nbdts = 1;
//char* bdtvars[nbdts] = {"top10Var"};
//int varindex[nbdts] = {0};
//const int nbdts = 23;
//char* bdtvars[nbdts] = {"allVar","noMbb","noDRjj","noHPt","noWPt","noJet2CSV","noTopMass","noHVdPhi","noNAddJet","noLepMetDPhi","noNSA5Jet","noJet1Pt","noJet2Pt","noLepPt","noLepIso","noJet1CSV","noLepEta","noJet1Eta","noJet2Eta","noPtBal","noAddJetCSV","noWMt","noMET"};
//int varindex[nbdts] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
const int nbdts = 11;
char* bdtvars[nbdts] = {"allVar","noMbb","noHPt","noVPt","noJet2CSV","noTopMass","noHVdPhi","noNAddJet","noLepMetDPhi","noNSA5Jet","noMET"};
int varindex[nbdts] = {0,1,3,4,5,6,7,8,9,10,22};

gROOT->ProcessLine("std::vector<const char*> weights;");
gROOT->ProcessLine("std::vector<const char*> branches;");
for (int i=0; i<nbdts; i++) {
    gROOT->ProcessLine("weights = std::vector<const char*>();");
    gROOT->ProcessLine("branches = std::vector<const char*>();");
    gROOT->ProcessLine(Form("weights.push_back(\"../../scripts/weights/TMVA_13TeV_10Nm1_%s_150_3_H125Sig_0b1b2bWjetsTTbarBkg_Mjj_BDT.weights.xml\")",bdtvars[i]));
    gROOT->ProcessLine(Form("branches.push_back(\"BDT_10Nm1_%s\")",bdtvars[i]));
    gROOT->ProcessLine(Form("addBDTToTree(\"%s\",\"%s_2\",weights,branches,%i)",fname,fname,varindex[i]));
    //gROOT->ProcessLine(Form("mv %s_2 %s;",fname,fname));
    gROOT->ProcessLine(Form("gSystem->CopyFile(\"%s_2\",\"%s\",true)",fname,fname));
    gROOT->ProcessLine("weights.clear()");
    gROOT->ProcessLine("branches.clear()");
}

}
