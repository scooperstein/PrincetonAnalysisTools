void makeBDTParamRockCurves(const char* fname) {

gROOT->SetBatch(true);
gROOT->ProcessLine(".L ../../../scripts/makeBDTRockCurve.C");

const int nTreeIter = 14;
const int nMaxDepthIter = 1;
int nTreeArr[nTreeIter] = {150,200,250,300,350,400,450,500,550,600,650,700,750,800};
int maxDepthArr[nMaxDepthIter] = {3};

//const int nTreeIter = 3;
//const int nMaxDepthIter = 2;
//int nTreeArr[nTreeIter] = {150,300,500};
//int maxDepthArr[nMaxDepthIter] = {3,4};

gROOT->ProcessLine("std::vector<string> bdtnames;");
for (int i=0; i<nTreeIter; i++) {
    for (int j=0; j<nMaxDepthIter; j++) {
        gROOT->ProcessLine("bdtnames.push_back(\"BDT_May6_150_3\")");
        gROOT->ProcessLine(Form("bdtnames.push_back(\"BDT_May6_%i_%i\")",nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,0,\"%i_%i\")",fname,nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,1,\"%i_%i\")",fname,nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,2,\"%i_%i\")",fname,nTreeArr[i],maxDepthArr[j]));
        gROOT->ProcessLine("bdtnames.clear()");
    }
}
/*gROOT->ProcessLine("bdtnames.push_back(\"BDT_April28_150_3");
gROOT->ProcessLine("bdtnames.push_back(\"BDT_April28_300_3");
gROOT->ProcessLine("bdtnames.push_back(\"BDT_April28_500_4");
gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,0,\"all\")",fname));
gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,1,\"all\")",fname));
gROOT->ProcessLine(Form("makeBDTRockCurve(\"%s\",\"tree\",bdtnames,2,\"all\")",fname));
*/
}
