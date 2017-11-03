#include <string>
#include <vector>
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
using namespace std;
  
#include <string>

void addBDTToTree_ewk(const char* oldfilename="testingtree_mh125_evenbkg_newmvas.root",const char* newfilename="testingtree_mh125_evenbkg_allmvas.root", std::vector<const char*> weights=std::vector<const char*>(),std::vector<const char*> newbranch=std::vector<const char*>(),int Nm1Var=0, bool refill=false) {

  using std::max;
  using std::min;
  
  int nBDTs = (int) weights.size();
  std::vector<TMVA::Reader*>  tmvaReader_BDTs; 

  ULong64_t evt;
  Float_t   Mjj,HJ1_HJ2_dR,PtJJ,J1_J2_dEta,Jet1_pt,Jet2_pt,Jet1_eta,Jet2_eta,Jet1_qgl,hJets_pt_0,hJets_pt_1,selLeptons_pt_0,selLeptons_relIso_0,hJets_btagCSV_0,hJets_eta_0,hJets_eta_1,jjWPtBalance,AddJets252p9_puid_leadJet_btagCSV,AddJets252p9_puid_leadJet_pt,selLeptons_eta_0,zW,Jet2_qgl_f,HVdEta,HVdEta_4MET,JJEtaBal,PtJJ_over_Mjj,J1_J2_dEta_over_Mjj,TopMass_over_Mjj,Jet1_eta_over_Mjj,zW_over_Mjj,Jet2_qgl;
  //Float_t  *Jet_btagCSV,*Jet_pt_reg,*selLeptons_pt,*selLeptons_eta,*Jet_eta;

  //Int_t isWenu;
  Float_t BDT[nBDTs];

  for (int i=0; i<nBDTs; i++) {
      /*tmvaReader_BDTs.push_back(new TMVA::Reader("V:!Color:Silent"));
      tmvaReader_BDTs[i]->AddVariable("Mjj",&Mjj);
      tmvaReader_BDTs[i]->AddVariable("PtJJ", &PtJJ);
      tmvaReader_BDTs[i]->AddVariable("J1_J2_dEta", &J1_J2_dEta);
      tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd1]", &hJets_btagCSV_0);
      tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd2]", &Jet1_pt);
      tmvaReader_BDTs[i]->AddVariable("Jet1_eta", &Jet1_eta);
      tmvaReader_BDTs[i]->AddVariable("Jet2_eta", &Jet2_eta);
      tmvaReader_BDTs[i]->AddVariable("nAddLep_f", &nAddLep_f);
      tmvaReader_BDTs[i]->AddVariable("Jet_pt[hJetInd1]", &hJets_pt_0);
      tmvaReader_BDTs[i]->AddVariable("Jet_pt[hJetInd2]", &hJets_pt_1);
      tmvaReader_BDTs[i]->AddVariable("zW", &zW);
      tmvaReader_BDTs[i]->AddVariable("Top1_mass_fromLepton", &Top1_mass_fromLepton);
      tmvaReader_BDTs[i]->AddVariable("Jet1_qgl", &Jet1_qgl);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_pt[lepInd]", &selLeptons_pt_0);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_eta[lepInd]", &selLeptons_eta_0);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_relIso03[lepInd]", &selLeptons_relIso03_0);
      tmvaReader_BDTs[i]->AddVariable("isWenu", &isWenu_f);
      tmvaReader_BDTs[i]->AddVariable("isWmunu_f", &isWmunu_f);*/

      tmvaReader_BDTs.push_back(new TMVA::Reader("V:!Color:Silent"));
      if (Nm1Var!=1) {
        tmvaReader_BDTs[i]->AddVariable("Mjj", &Mjj);
      }
      //if (Nm1Var!=2) {
      //  tmvaReader_BDTs[i]->AddVariable("HJ1_HJ2_dR", &HJ1_HJ2_dR);
      //}
      if (Nm1Var!=3) {
        tmvaReader_BDTs[i]->AddVariable("PtJJ",   &PtJJ);
      }
      if (Nm1Var!=4) {
        tmvaReader_BDTs[i]->AddVariable("J1_J2_dEta",   &J1_J2_dEta);
      }
      //if (Nm1Var!=5) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet1_pt", &Jet1_pt);
      //  //tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd2]", &Jet1_pt);
      //}
      //if (Nm1Var!=6) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet2_pt", &Jet2_pt);
      //  //tmvaReader_BDTs[i]->AddVariable("Jet2_pt", &Jet2_pt);
      //  //tmvaReader_BDTs[i]->AddVariable("Jet2_pt"", &Jet2_pt");
      //}
      //if (Nm1Var!=7) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet1_eta", &Jet1_eta);
     // }
      //if (Nm1Var!=8) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet2_eta", &Jet2_eta);
      //  //tmvaReader_BDTs[i]->AddVariable("nAddJets252p9_puid", &Jet2_eta);
     // }
      if (Nm1Var!=9) {
        tmvaReader_BDTs[i]->AddVariable("max(0.,Jet1_qgl)", &Jet1_qgl);
      }
      if (Nm1Var!=10) {
        //tmvaReader_BDTs[i]->AddVariable("Jet2_qgl", &Jet2_qgl);
        tmvaReader_BDTs[i]->AddVariable("max(0.,Jet2_qgl)", &Jet2_qgl_f);
      }
      //if (Nm1Var!=11) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet_pt_reg[hJetInd1]", &hJets_pt_0);
      //}
      //if (Nm1Var!=12) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet_pt_reg[hJetInd2]", &hJets_pt_1);
      //}
      //if (Nm1Var!=13) {
      //  tmvaReader_BDTs[i]->AddVariable("selLeptons_pt[lepInd]", &selLeptons_pt_0);
      //}
      //if (Nm1Var!=14) {
      //  tmvaReader_BDTs[i]->AddVariable("selLeptons_relIso_0", &selLeptons_relIso_0);
     // }
     // if (Nm1Var!=15) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd1]", &hJets_btagCSV_0);
      //}
      //if (Nm1Var!=16) {
      //  tmvaReader_BDTs[i]->AddVariable("selLeptons_eta[lepInd]", &selLeptons_eta_0);
      //}
      //if (Nm1Var!=17) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet_eta[hJetInd1]", &hJets_eta_0);
      //}
      //if (Nm1Var!=18) {
      //  tmvaReader_BDTs[i]->AddVariable("Jet_eta[hJetInd2]", &hJets_eta_1);
      //}
      //if (Nm1Var!=19) {
      //  tmvaReader_BDTs[i]->AddVariable("jjWPtBalance", &jjWPtBalance);
      //}
      //if (Nm1Var!=20) {
      //  tmvaReader_BDTs[i]->AddVariable("AddJets252p9_puid_leadJet_btagCSV", &AddJets252p9_puid_leadJet_btagCSV);
      //}
      //tmvaReader_BDTs[i]->AddVariable("AddJets252p9_puid_leadJet_pt", &AddJets252p9_puid_leadJet_pt);
      //if (Nm1Var!=21) {
      //  //tmvaReader_BDTs[i]->AddVariable("hJets_btagCSV_0", &hJets_btagCSV_0);
      //  tmvaReader_BDTs[i]->AddVariable("selLeptons_eta_0", &selLeptons_eta_0);
     // }
      if (Nm1Var!=22) {
        tmvaReader_BDTs[i]->AddVariable("min(abs(zW),5)", &zW);
      }
      //tmvaReader_BDTs[i]->AddVariable("HVdEta_4MET", &HVdEta_4MET);
      //tmvaReader_BDTs[i]->AddVariable("JJEtaBal", &JJEtaBal);
      tmvaReader_BDTs[i]->BookMVA("EWKWJets", weights[i]);
      cout<<"finished tmvaReader_BDT "<<i<<endl;
  }

  TFile *oldfile = TFile::Open(oldfilename);
  
  TTree *oldtree = (TTree*)oldfile->Get("tree");

  oldtree->SetBranchAddress("evt", &evt);
  oldtree->SetBranchAddress("Mjj", &Mjj);
  oldtree->SetBranchAddress("PtJJ", &PtJJ);
  oldtree->SetBranchAddress("J1_J2_dEta", &J1_J2_dEta);
  oldtree->SetBranchAddress("Jet1_pt", &Jet1_pt);
  oldtree->SetBranchAddress("Jet2_pt", &Jet2_pt);
  oldtree->SetBranchAddress("Jet1_eta", &Jet1_eta);
  oldtree->SetBranchAddress("Jet2_eta", &Jet2_eta);
  oldtree->SetBranchAddress("Jet1_qgl", &Jet1_qgl);
  oldtree->SetBranchAddress("Jet2_qgl", &Jet2_qgl);
  oldtree->SetBranchAddress("selLeptons_eta_0", &selLeptons_eta_0);
  oldtree->SetBranchAddress("zW", &zW);
  Long64_t nentries = oldtree->GetEntries();

  TFile *newfile = TFile::Open(newfilename, "recreate");
  TTree *newtree = oldtree->CloneTree(0);

  std::vector<TBranch*> bBDT;
  for (int i=0; i<nBDTs; i++) {
    if(refill){
      oldtree->SetBranchAddress(newbranch[i],  &BDT[i]);
    } else {  
      char nametype[100];
      strcpy(nametype,newbranch[i]);
      strcat(nametype,"/F");
    
      bBDT.push_back(newtree->Branch(newbranch[i],&BDT[i]));
    }
  }
  
  //nentries = 1000;
  cout<<"starting loop over "<<nentries<<" events"<<endl;
  for (Long64_t z=0; z<nentries; z++) {
  //for (Long64_t z=0; z<50000; z++) {
    oldtree->GetEntry(z);
    //if (evt!=1537680) continue;
    //isWenu_f = (float) isWenu;
    zW = fabs(zW);
    if (zW > 5.) zW = 5;
    if (Jet1_qgl < 0) Jet1_qgl = 0.;
    if (Jet2_qgl < 0) Jet2_qgl = 0.;

    //BDT  = tmvaReader_BDT->EvaluateMVA("Gradient");
    for (int i=0; i<nBDTs; i++) {       
        BDT[i]  = tmvaReader_BDTs[i]->EvaluateMVA("EWKWJets");       

        if(z%50000==1) cout<<z<<"  "<<BDT[i]<<std::endl;
        //cout<<"BDT = "<<BDT[i]<<std::endl;
    }
    //std::cout<<"Mjj = "<<Mjj<<std::endl;
    //std::cout<<"PtJJ = "<<PtJJ<<std::endl;
    //std::cout<<"J1_J2_dEta = "<<J1_J2_dEta<<std::endl;
    //std::cout<<"Jet1_pt = "<<Jet1_pt<<std::endl;
    //std::cout<<"Jet2_pt = "<<Jet2_pt<<std::endl;
    //std::cout<<"Jet1_eta = "<<Jet1_eta<<std::endl;
    //std::cout<<"Jet2_eta = "<<Jet2_eta<<std::endl;
    //std::cout<<"Jet1_qgl = "<<Jet1_qgl<<std::endl;
    //std::cout<<"Jet2_qgl = "<<Jet2_qgl<<std::endl;
    //std::cout<<"selLeptons_eta_0 = "<<selLeptons_eta_0<<std::endl;
    //std::cout<<"zW = "<<zW<<std::endl;
    //std::cout<<"BDT = "<<BDT[0]<<std::endl;
    newtree->Fill();
  }

  //newtree->Print();
  //newtree->AutoSave();
  newfile->cd();
  newtree->Write();
  delete oldfile;
  delete newfile;
}

int main(int argc, char *argv[]) {
    
    const char* file1 = argv[1];
    const char* file2 = argv[2];
    std::vector<const char*> weights;
    std::vector<const char*> branches;
    weights.push_back("weights/TMVA_13TeV_V27_EWK_Aug30_6Var_EWKWJets.weights.xml");
    branches.push_back("BDT_V27_EWK_Aug30_6Var");

   addBDTToTree_ewk(file1,file2,weights,branches);
}
