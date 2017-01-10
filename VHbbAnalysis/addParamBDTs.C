void addParamBDTs(const char* fname) {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L ../../scripts/addBDTToTree.C ");

const int nTreeIter = 14;
//const int nTreeIter = 10;
const int nMaxDepthIter = 3;
//int nTreeArr[nTreeIter] = {650,700,750,800};
int nTreeArr[nTreeIter] = {150,200,250,300,350,400,450,500,550,600,650,700,750,800};
//int nTreeArr[nTreeIter] = {350,400,450,500,550,600,650,700,750,800};
//int nTreeArr[nTreeIter] = {150,200,250,300,350};
int maxDepthArr[nMaxDepthIter] = {3,4,5};
//const int nTreeIter = 1;
//const int nMaxDepthIter = 1;
//int nTreeArr[nTreeIter] = {300};
//int maxDepthArr[nMaxDepthIter] = {3};

gROOT->ProcessLine("std::vector<const char*> weights;");
gROOT->ProcessLine("std::vector<const char*> branches;");
for (int i=0; i<nTreeIter; i++) {
    for (int j=0; j<nMaxDepthIter; j++) {
        //if ((i==0 && j == 0) || (i==0 && j==1)) continue;
        gROOT->ProcessLine(Form("weights.push_back(\"../../scripts/weights/TMVA_13TeV_V24_Nov16_%i_%i_VVSig_0b1b2bWjetsTTbarVVBkg_Mjj_BDT.weights.xml\")",nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine(Form("branches.push_back(\"BDT_V24_Nov16_%i_%i_VVSig\")",nTreeArr[i],maxDepthArr[j]));
    }
}

gROOT->ProcessLine(Form("addBDTToTree(\"%s\",\"%s_2\",weights,branches,0)",fname,fname));
//gROOT->ProcessLine(Form("gSystem->CopyFile(\"%s_2\",\"%s\",true)",fname,fname));

}
