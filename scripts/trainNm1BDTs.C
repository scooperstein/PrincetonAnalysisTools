void trainNm1BDTs() {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L OLD_TMVAClassification_8TeV.C");

const int nTreeIter = 14;
const int nMaxDepthIter = 3;
int nTreeArr[nTreeIter] = {150,200,250,300,350,400,450,500,550,600,650,700,750,800};
int maxDepthArr[nMaxDepthIter] = {3,4,5};

for (int i=0; i<nTreeIter; i++) {
    for (int j=0; j<nMaxDepthIter; j++) {
        gROOT->ProcessLine(Form("OLD_TMVAClassification_8TeV(0, %i, 0.05, %i, 0.1, 0, \"\")",nTreeArr[i],maxDepthArr[j]));
    }
}

}
