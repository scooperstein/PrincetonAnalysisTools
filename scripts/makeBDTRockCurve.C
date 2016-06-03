void makeBDTRockCurve(char* filename, char *treename, std::vector<std::string> bdtnames=std::vector<std::string>(), int zoom=0) {

TFile *ifile = new TFile(filename, "r");
TTree *tree = (TTree*) ifile->Get(treename);

//char* presel = "hJets_btagCSV_1>0.6";
//char* presel = "hJets_btagCSV_0>0.85 && hJets_btagCSV_1>0.6 && abs(HVdPhi)>2.5 && abs(lepMetDPhi)<2 && nAddJet_f<2 && nAddLep_f<2";
char* presel = "H_mass>90&&H_mass<150&&V_pt>150";

int nBins = 1000;
int nBDTs = (int) bdtnames.size();

std::vector<TH1F*> h_bdtSig_test;
std::vector<TH1F*> h_bdtBkg_test;
std::vector<TH1F*> h_bdtSig_train;
std::vector<TH1F*> h_bdtBkg_train;
std::vector<Double_t*> sig_test;
std::vector<Double_t*> bkg_test;
std::vector<Double_t*> sig_train;
std::vector<Double_t*> bkg_train;
for (int i=0; i<nBDTs; i++) {
    h_bdtSig_test.push_back(new TH1F(Form("h_bdtSig_test_%s",bdtnames[i].c_str()),Form("h_bdtSig_test_%s",bdtnames[i].c_str()),nBins,-1,1)); 
    h_bdtBkg_test.push_back(new TH1F(Form("h_bdtBkg_test_%s",bdtnames[i].c_str()),Form("h_bdtBkg_test_%s",bdtnames[i].c_str()),nBins,-1,1));
    h_bdtSig_train.push_back(new TH1F(Form("h_bdtSig_train_%s",bdtnames[i].c_str()),Form("h_bdtSig_train_%s",bdtnames[i].c_str()),nBins,-1,1)); 
    h_bdtBkg_train.push_back(new TH1F(Form("h_bdtBkg_train_%s",bdtnames[i].c_str()),Form("h_bdtBkg_train_%s",bdtnames[i].c_str()),nBins,-1,1));
    tree->Draw(Form("%s>>h_bdtSig_test_%s",bdtnames[i].c_str(),bdtnames[i].c_str()),Form("(sampleIndex<0&&evt%2==0&&(%s))*weight",presel));
    tree->Draw(Form("%s>>h_bdtBkg_test_%s",bdtnames[i].c_str(),bdtnames[i].c_str()),Form("(sampleIndex>0&&evt%2==0&&(%s))*weight",presel));
    tree->Draw(Form("%s>>h_bdtSig_train_%s",bdtnames[i].c_str(),bdtnames[i].c_str()),Form("(sampleIndex<0&&evt%2==1&&(%s))*weight",presel));
    tree->Draw(Form("%s>>h_bdtBkg_train_%s",bdtnames[i].c_str(),bdtnames[i].c_str()),Form("(sampleIndex>0&&evt%2==1&&(%s))*weight",presel));
    sig_test.push_back(h_bdtSig_test[i]->GetIntegral());
    bkg_test.push_back(h_bdtBkg_test[i]->GetIntegral());
    sig_train.push_back(h_bdtSig_train[i]->GetIntegral());
    bkg_train.push_back(h_bdtBkg_train[i]->GetIntegral());
}

Double_t sig_den = h_bdtSig_test[0]->Integral() + h_bdtSig_train[0]->Integral();
Double_t bkg_den = h_bdtBkg_test[0]->Integral() + h_bdtBkg_test[0]->Integral();

std::cout<<"sig_den = "<<sig_den<<std::endl;;
std::cout<<"bkg_den = "<<bkg_den<<std::endl;;

// we want to cut > BDT value, GetIntegral() assumes <
for (int i=0; i<nBDTs; i++) { 
    for (int j=0; j<nBins; j++) {
        sig_test[i][j] = 1 - sig_test[i][j];
        bkg_test[i][j] = 1 - bkg_test[i][j];
        sig_train[i][j] = 1 - sig_train[i][j];
        bkg_train[i][j] = 1 - bkg_train[i][j];
        //if (sig[i][j] > 0.01 && sig[i][j] < 1.) {
        //   //std::cout<<"sig eff: "<<sig[i]<<", bkg eff: "<<bkg[i]<<std::endl;
        //}
    }
}

TMultiGraph *mg = new TMultiGraph();
TLegend *leg = new TLegend(0.1,0.5,0.5,0.9);

//const TColor *colors = [kBlue, kRed, kGreen, kCyan, kViolet, kOrange, kMagenta, kYellow];
const Int_t colors[8] = {880, 432, 600, 394, 418, 616, 808, 400};

for (int i=0; i < nBDTs; i++) {
    TGraph *gb_test = new TGraph(nBins, sig_test[i], bkg_test[i]);
    TGraph *gb_train = new TGraph(nBins, sig_train[i], bkg_train[i]);
    gb_test->SetLineColor(colors[i]);
    gb_train->SetLineColor(colors[i]);
    //gb->SetMarkerStyle(0);
    gb_test->SetLineWidth(2);
    gb_train->SetLineWidth(2);
    gb_train->SetLineStyle(2);
    mg->Add(gb_test,"LP");
    mg->Add(gb_train,"LP");
    leg->AddEntry(gb_test, Form("Test %s",bdtnames[i].c_str()));
    leg->AddEntry(gb_train, Form("Train %s",bdtnames[i].c_str()));
}
// add some cut-based WP's
TGraph *gwp = new TGraph();
//gwp->SetPoint(0,0.00034,0.0072);
//gwp->SetPoint(0,0.260, 0.0193);
//gwp->SetPoint(1,0.1647,0.006803);
//gwp->SetPoint(2,0.1488,0.005470);
//gwp->SetPoint(3,0.1209,0.003662);

//gwp->SetPoint(0,0.6159,0.3502);
//gwp->SetPoint(1,0.536,0.2767);
//gwp->SetPoint(2,0.39,0.1236);
//gwp->SetPoint(3,0.3523,0.09935);
//gwp->SetPoint(4,0.2694,0.06652);
//gwp->SetPoint(5,0.2088,0.03458);
//gwp->SetPoint(6,0.1358,0.01482);
//gwp->SetPoint(7,0.0777,0.00476);
//gwp->SetPoint(8,0.0537,0.00218);
//gwp->SetPoint(9,0.0399,0.00127);
//gwp->SetPoint(10,0.0279,0.0007);
//gwp->SetPoint(11,0.0211,0.000527);
//gwp->SetPoint(12,0.0129,0.000260);

gwp->SetPoint(0,10.23/sig_den,2494.9/bkg_den);
gwp->SetPoint(1,8.91/sig_den,1971.3/bkg_den);
gwp->SetPoint(2,6.478/sig_den,880.44/bkg_den);
gwp->SetPoint(3,5.851/sig_den,707.84/bkg_den);
gwp->SetPoint(4,4.754/sig_den,473.97/bkg_den);
gwp->SetPoint(5,3.468/sig_den,246.4/bkg_den);
gwp->SetPoint(6,2.256/sig_den,105.6/bkg_den);
gwp->SetPoint(7,1.29/sig_den,33.88/bkg_den);
gwp->SetPoint(8,0.893/sig_den,15.51/bkg_den);
gwp->SetPoint(9,0.662/sig_den,9.022/bkg_den);
gwp->SetPoint(10,0.464/sig_den,4.987/bkg_den);
gwp->SetPoint(11,0.35/sig_den,3.757/bkg_den);
gwp->SetPoint(12,0.214/sig_den,1.853/bkg_den);

gwp->SetMarkerStyle(29);
gwp->SetMarkerSize(2);
//mg->Add(gwp,"LP");

//leg->AddEntry(gb,filename,"elp");
//leg->AddEntry(gb,"BDT","elp");
leg->AddEntry(gwp,"cut-based WP's", "p");

TCanvas *canv = new TCanvas();

mg->SetTitle("");
mg->SetMinimum(0.0);
if (zoom==1) {
  mg->SetMaximum(0.01);
}
else if (zoom==2) {
  mg->SetMaximum(0.1);
}

mg->Draw("ALP");

mg->GetYaxis()->SetTitle("Bkg Efficiency");
mg->GetXaxis()->SetTitle("Signal Efficiency");
if (zoom==1) {
  mg->GetXaxis()->SetLimits(0.0,0.1);
}
else if (zoom==2) {
  mg->GetXaxis()->SetLimits(0.0,0.3);
}
TGaxis::SetMaxDigits(3);
//mg->GetXaxis()->SetRange(990,1000);
leg->Draw("same");
//canv->Update();
canv->SaveAs(Form("%i.png",zoom));

}
