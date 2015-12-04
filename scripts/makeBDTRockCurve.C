void makeBDTRockCurve(char* filename, char *treename, char* bdtname, bool tree1FromTMVA=false, char* filename2="", char* treename2="", char* bdtname2="", bool tree2FromTMVA=false) {

TFile *ifile = new TFile(filename, "r");
TTree *tree = (TTree*) ifile->Get(treename);

char* presel = "hJets_btagCSV_1>0.6";

TTree *tree2 = new TTree();

if (filename2!="") {
    TFile *ifile2 = new TFile(filename2, "r");
    tree2 = (TTree*) ifile2->Get(treename2);    
}

int nBins = 1000;

TH1F *h_bdtSig = new TH1F("h_bdtSig","h_bdtSig",nBins,-1,1); 
TH1F *h_bdtBkg = new TH1F("h_bdtBkg","h_bdtBkg",nBins,-1,1);

TH1F *h2_bdtSig = new TH1F("h2_bdtSig","h2_bdtSig",nBins,-1,1); 
TH1F *h2_bdtBkg = new TH1F("h2_bdtBkg","h2_bdtBkg",nBins,-1,1);

if (tree1FromTMVA) {
    tree->Draw(Form("%s>>h_bdtSig",bdtname),Form("(classID==0&&(%s))*weight",presel));
    tree->Draw(Form("%s>>h_bdtBkg",bdtname),Form("(classID==1&&(%s))*weight",presel));
}
else {
    tree->Draw(Form("%s>>h_bdtSig",bdtname),Form("(sampleIndex<0&&(%s))*weight",presel));
    tree->Draw(Form("%s>>h_bdtBkg",bdtname),Form("(sampleIndex>0&&(%s))*weight",presel));
}
n2entries = tree2->GetEntries();
if (n2entries > 0) {
  if (tree2FromTMVA) {
      tree2->Draw(Form("%s>>h2_bdtSig",bdtname2),Form("(classID==0&&(%s))*weight",presel));
      tree2->Draw(Form("%s>>h2_bdtBkg",bdtname2),Form("(classID==1&&(%s))*weight",presel));
  }
  else {
      tree2->Draw(Form("%s>>h2_bdtSig",bdtname2),Form("(sampleIndex<0&&H_mass>90&&H_mass<150&&(%s))*weight",presel));
      tree2->Draw(Form("%s>>h2_bdtBkg",bdtname2),Form("(sampleIndex>0&&H_mass>90&&H_mass<150&&(%s))*weight",presel));
  }
}

Double_t *sig = h_bdtSig->GetIntegral();
Double_t *bkg = h_bdtBkg->GetIntegral();

Double_t *sig2 = h2_bdtSig->GetIntegral();
Double_t *bkg2 = h2_bdtBkg->GetIntegral();

// we want to cut > BDT value, GetIntegral() assumes < 
for (int i=0; i<nBins; i++) {
    sig[i] = 1 - sig[i];
    bkg[i] = 1 - bkg[i];
    if (sig[i] < 0.05 && sig[i] > 0) {
        std::cout<<"sig eff: "<<sig[i]<<", bkg eff: "<<bkg[i]<<std::endl;
    }
}

if (n2entries > 0) {

    // we want to cut > BDT value, GetIntegral() assumes < 
    for (int i=0; i<nBins; i++) {
        sig2[i] = 1 - sig2[i];
        bkg2[i] = 1 - bkg2[i];
        if (sig2[i] < 0.05 && sig2[i] > 0) {
            std::cout<<"sig2 eff: "<<sig2[i]<<", bkg2 eff: "<<bkg2[i]<<std::endl;
        }
    }
}

TMultiGraph *mg = new TMultiGraph();
TLegend *leg = new TLegend(0.1,0.5,0.5,0.9);
TGraph *gb = new TGraph(nBins, sig, bkg);
gb->SetLineColor(kBlue);
gb->SetMarkerColor(kMagenta);
gb->SetMarkerStyle(0);
gb->SetLineWidth(2);
mg->Add(gb,"LP");

// add some cut-based WP's
TGraph *gwp = new TGraph();
//gwp->SetPoint(0,0.00034,0.0072);
gwp->SetPoint(0,0.260, 0.0193);
gwp->SetPoint(1,0.1647,0.006803);
gwp->SetPoint(2,0.1488,0.005470);
gwp->SetPoint(3,0.1209,0.003662);
gwp->SetMarkerStyle(29);
gwp->SetMarkerSize(2);
mg->Add(gwp,"LP");

//leg->AddEntry(gb,filename,"elp");
leg->AddEntry(gb,"BDT","elp");
leg->AddEntry(gwp,"cut-based WP's", "p");

if (n2entries > 0) {
    TGraph *gb2 = new TGraph(nBins, sig2, bkg2);
    gb2->SetLineColor(kRed);
    gb2->SetMarkerColor(kGreen);
    gb2->SetMarkerStyle(0);
    gb2->SetLineWidth(2);
    mg->Add(gb2,"LP");  
    //leg->AddEntry(gb2,filename2,"elp");
    leg->AddEntry(gb2,"BDT no M_{jj}","elp");
}

//TCanvas *canv = new TCanvas();

mg->SetTitle("");
//mg->SetMinimum();
//mg->SetMaximum(0.1);

mg->Draw("ALP");

mg->GetYaxis()->SetTitle("Bkg Efficiency");
mg->GetXaxis()->SetTitle("Signal Efficiency");
//mg->GetXaxis()->SetLimits(0.0,0.001);
TGaxis::SetMaxDigits(3);
//mg->GetXaxis()->SetRange(990,1000);
leg->Draw("same");
//canv->Update();

}
