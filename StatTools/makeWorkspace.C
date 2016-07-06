//
// Perform mits on the Mjj distribution for a combination of signal
// and background templates. Used for study of alternative Mjj-fitting
// analysis strategy.
//
// Author: Stephane Cooperstein
//


void makeWorkspace(char* channel) {

TFile *ifile = TFile::Open(Form("hists_%s.root",channel), "r");
TH1F *hist_bkg = (TH1F*) ifile->Get(Form("BDT_%s_Bkg",channel));
TH1F *hist_wh = (TH1F*) ifile->Get(Form("BDT_%s_WH",channel));
TH1F *hist_zh = (TH1F*) ifile->Get(Form("BDT_%s_ZH",channel));

RooWorkspace w(Form("ws_%s",channel));

RooRealVar mass("mass","mass",50,200);
RooDataHist data("data","data",mass,hist_bkg);
RooDataHist wh("wh","wh",mass,hist_wh);
RooDataHist zh("zh","zh",mass,hist_zh);

RooRealVar a0("a0","a0",-0.009,-10000,10000);
RooRealVar a1("a1","a1",0.0002,-10000,10000);
RooRealVar a2("a2","a2",0.0001,-10000,10000);
RooRealVar a3("a3","a3",0.0001,-10000,10000);
RooRealVar a4("a4","a4",0.000001,-100,100);
RooRealVar a5("a5","a5",0.01,-100,100);
RooRealVar a6("a6","a6",0.01,-100,100);
RooRealVar a7("a7","a7",0.01,-100,100);
RooRealVar a8("a8","a8",0.01,-100,100);
RooRealVar a9("a9","a9",0.01,-100,100);
RooRealVar a10("a10","a10",0.5,-100,100);
RooRealVar a11("a11","a11",0.5,-100,100);
RooRealVar a12("a12","a12",0.5,-100,100);
RooRealVar a13("a13","a13",0.5,-100,100);
RooRealVar a14("a14","a14",0.5,-100,100);
RooRealVar a15("a15","a15",-0.5,-100,100);

RooRealVar peak("peak","peak",100,20,500);
RooRealVar width("width","width",40,1,100);
RooRealVar tail("tail","tail",0.00001,-100,100);

RooArgList arglist(a0,a1,a2,a3);
//RooArgList arglist(a0,a1,a2,a3,a4,a5,a6,a7,a8);
//arglist.add(a9);
//arglist.add(a10);
//arglist.add(a11);
//arglist.add(a12);
//arglist.add(a13);
RooPolynomial poly("poly","poly",mass,arglist);

RooRealVar lambda("lambda", "slope", -0.05, -0.2., 0.);
RooExponential expo("expo", "exponential PDF", mass, lambda);

RooProdPdf bkgmodel("bkgmodel","bkgmodel",poly,expo);

RooNovosibirsk nb("nb","nb",mass,peak,width,tail);

RooRealVar mean("mean", "mean" , 100, 20, 180.0) ;
RooRealVar sigma("sigma", "sigma" , 30, 1.0, 100.0) ;
RooRealVar alpha("alpha", "alpha", -1, -100, 100);
RooRealVar n("n","n",0,-100,100);
RooCBShape CBall("CBall", "Crystal Ball shape", mass, mean, sigma, alpha, n);

//RooAbsReal CBall_norm("CBall_norm","Cball_norm",1300,0,100000);

//poly.fitTo(data,SumW2Error(kFALSE));
//expo.fitTo(data,SumW2Error(kFALSE));
//bkgmodel.fitTo(data,SumW2Error(kFALSE));
nb.fitTo(data,SumW2Error(kFALSE));

RooPlot* frame = mass.frame();
data.plotOn(frame);
nb.plotOn(frame);
//expo.plotOn(frame);
//bkgmodel.plotOn(frame);
//CBall.plotOn(frame);
cout<<"chi2/ndof = "<<frame->chiSquare()<<endl;
frame->Draw();


w.import(data);
w.import(wh);
w.import(zh);
w.import(nb);

w.writeToFile(Form("workspace_%s.root",channel));

}
