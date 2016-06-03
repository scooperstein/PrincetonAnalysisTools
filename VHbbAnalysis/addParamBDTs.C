void addParamBDTs(const char* fname) {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L ../../scripts/addBDTToTree.C ");

const int nTreeIter = 4;
const int nMaxDepthIter = 3;
int nTreeArr[nTreeIter] = {650,700,750,800};
//int nTreeArr[nTreeIter] = {150,200,250,300,350,400,450,500,550,600,650,700,750,800};
int maxDepthArr[nMaxDepthIter] = {3,4,5};
//const int nTreeIter = 1;
//const int nMaxDepthIter = 1;
//int nTreeArr[nTreeIter] = {300};
//int maxDepthArr[nMaxDepthIter] = {3};

gROOT->ProcessLine("std::vector<const char*> weights;");
gROOT->ProcessLine("std::vector<const char*> branches;");
for (int i=0; i<nTreeIter; i++) {
    for (int j=0; j<nMaxDepthIter; j++) {
        gROOT->ProcessLine(Form("weights.push_back(\"../../scripts/weights/TMVA_13TeV_May23_%i_%i_H125Sig_0b1b2bWjetsTTbarBkg_Mjj_BDT.weights.xml\")",nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine(Form("branches.push_back(\"BDT_May23_%i_%i\")",nTreeArr[i],maxDepthArr[j]));
    }
}

gROOT->ProcessLine(Form("addBDTToTree(\"%s\",\"%s_2\",weights,branches,0)",fname,fname));
gROOT->ProcessLine(Form("gSystem->CopyFile(\"%s_2\",\"%s\",true)",fname,fname));

}
