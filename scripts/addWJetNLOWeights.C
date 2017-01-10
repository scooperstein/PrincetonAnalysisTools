double LOtoNLOWeightBjetSplitEtabb(double etabb, int njets){

    double SF = 1.;
    if(etabb < 5){
        if(njets < 1){
            SF =   0.935422 + 0.0403162*etabb -0.0089026*etabb*etabb +0.0064324*etabb*etabb*etabb -0.000212443*etabb*etabb*etabb*etabb;
        }else if(njets == 1){
            SF =   0.962415 +0.0329463*etabb -0.0414479*etabb*etabb +0.0240993*etabb*etabb*etabb -0.00278271*etabb*etabb*etabb*etabb;
        }else if(njets >=2){
            SF =   (0.721265 -0.105643*etabb -0.0206835*etabb*etabb +0.00558626*etabb*etabb*etabb)*TMath::Exp(0.450244*etabb);
        }
    }

    return SF;

}


void addWJetNLOWeights(char* ifname, char* ofname) {

TFile *ifile = TFile::Open(ifname);
TTree *tree = (TTree*) ifile->Get("tree");
Int_t nentries = tree->GetEntries();
//Int_t nentries = 10;
Float_t WJetNLOWeight;
Int_t hJetInd1;
Int_t hJetInd2;
Int_t nJet;
Float_t Jet_eta[100];
Int_t sampleIndex;
Float_t weight_ptQCD;

tree->SetBranchAddress("hJetInd1",&hJetInd1);
tree->SetBranchAddress("hJetInd2",&hJetInd2);
tree->SetBranchAddress("nJet",&nJet);
tree->SetBranchAddress("sampleIndex",&sampleIndex);
tree->SetBranchAddress("Jet_eta",Jet_eta);
tree->SetBranchAddress("weight_ptQCD",&weight_ptQCD);

TFile *ofile = TFile::Open(ofname,"RECREATE");
TTree *otree = tree->CloneTree(0);
otree->Branch("WJetNLOWeight",&WJetNLOWeight,"WJetNLOWeight/F");

cout<<"Processing Total Number of Events: "<<nentries<<endl;

for (i=0; i<nentries; i++) {
    if (i%10000 == 0) {
        cout<<"Processing entry: "<<i<<endl;
    }
    tree->GetEntry(i);
    float deta_bb = fabs(Jet_eta[hJetInd1] - Jet_eta[hJetInd2]);
    if (sampleIndex<4100 || sampleIndex>4902) { WJetNLOWeight = 1.0; }
    else if (sampleIndex==4100 || sampleIndex==4200 || sampleIndex==4300 || sampleIndex==4400 || sampleIndex==4500 || sampleIndex==4600 || sampleIndex==4700 || sampleIndex==4800 || sampleIndex==4900) {
        WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 0);
        WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
    }
    else if (sampleIndex==4101 || sampleIndex==4201 || sampleIndex==4301 || sampleIndex==4401 || sampleIndex==4501 || sampleIndex==4601 || sampleIndex==4701 || sampleIndex==4801 || sampleIndex==4901) {
        WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 1);
        WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
    }
    else if (sampleIndex==4102 || sampleIndex==4202 || sampleIndex==4302 || sampleIndex==4402 || sampleIndex==4502 || sampleIndex==4602 || sampleIndex==4702 || sampleIndex==4802 || sampleIndex==4902) {
        WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 0);
        WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
    }
    else {
        cout<<"didn't fall into any category, setting to 1..."<<endl;
        WJetNLOWeight = 1.0;
    }
    ofile->cd();
    otree->Fill(); 
    
}

ofile->cd();
otree->Write();

}
