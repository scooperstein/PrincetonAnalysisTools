void addBDTToTree(const char* oldfilename="testingtree_mh125_evenbkg_newmvas.root",const char* newfilename="testingtree_mh125_evenbkg_allmvas.root", std::vector<const char*> weights=std::vector<const char*>(),std::vector<const char*> newbranch=std::vector<const char*>(),int Nm1Var=0, bool refill=false) {
  
#include <string>


  using std::max;
  using std::min;
  
  int nBDTs = (int) weights.size();
  std::vector<TMVA::Reader*>  tmvaReader_BDTs; 

  //Float_t   H_mass, H_pt, V_pt, hJets_btagCSV_0, hJets_btagCSV_1, HVdPhi, nAddJet_f, nAddLep_f, hJets_pt_0, hJets_pt_1, met_pt, Top1_mass_fromLepton, lepMetDPhi, selLeptons_pt_0, selLeptons_eta_0, selLeptons_relIso03_0, isWenu_f, isWmunu_f;
  Float_t   H_mass,HJ1_HJ2_dR,H_pt,V_pt,hJets_btagCSV_1,Top1_mass_fromLepton_regPT_w4MET,HVdPhi,nAddJet_f,lepMetDPhi,hJets_pt_0,hJets_pt_1,selLeptons_pt_0,selLeptons_relIso_0,hJets_btagCSV_0,selLeptons_eta_0,hJets_eta_0,hJets_eta_1,jjWPtBalance,AddJets252p9_puid_leadJet_btagCSV,AddJets252p9_puid_leadJet_pt,V_mt,met_pt,softActivityVH_njets5_f,HVdEta,HVdEta_4MET,JJEtaBal;
  Int_t     hJetInd1,hJetInd2,lepInd,softActivityVH_njets5;
  //Float_t  *Jet_btagCSV,*Jet_pt_reg,*selLeptons_pt,*selLeptons_eta,*Jet_eta;
  Float_t Jet_btagCSV[100],Jet_pt_reg[100],selLeptons_pt[10],selLeptons_eta[10],Jet_eta[100];

  //Int_t isWenu;
  Float_t BDT[nBDTs];

  for (int i=0; i<nBDTs; i++) {
      /*tmvaReader_BDTs.push_back(new TMVA::Reader("V:!Color:Silent"));
      tmvaReader_BDTs[i]->AddVariable("H_mass",&H_mass);
      tmvaReader_BDTs[i]->AddVariable("H_pt", &H_pt);
      tmvaReader_BDTs[i]->AddVariable("V_pt", &V_pt);
      tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd1]", &hJets_btagCSV_0);
      tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd2]", &hJets_btagCSV_1);
      tmvaReader_BDTs[i]->AddVariable("HVdPhi", &HVdPhi);
      tmvaReader_BDTs[i]->AddVariable("nAddJet_f", &nAddJet_f);
      tmvaReader_BDTs[i]->AddVariable("nAddLep_f", &nAddLep_f);
      tmvaReader_BDTs[i]->AddVariable("Jet_pt[hJetInd1]", &hJets_pt_0);
      tmvaReader_BDTs[i]->AddVariable("Jet_pt[hJetInd2]", &hJets_pt_1);
      tmvaReader_BDTs[i]->AddVariable("met_pt", &met_pt);
      tmvaReader_BDTs[i]->AddVariable("Top1_mass_fromLepton", &Top1_mass_fromLepton);
      tmvaReader_BDTs[i]->AddVariable("lepMetDPhi", &lepMetDPhi);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_pt[lepInd]", &selLeptons_pt_0);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_eta[lepInd]", &selLeptons_eta_0);
      tmvaReader_BDTs[i]->AddVariable("selLeptons_relIso03[lepInd]", &selLeptons_relIso03_0);
      tmvaReader_BDTs[i]->AddVariable("isWenu", &isWenu_f);
      tmvaReader_BDTs[i]->AddVariable("isWmunu_f", &isWmunu_f);*/

      tmvaReader_BDTs.push_back(new TMVA::Reader("V:!Color:Silent"));
      if (Nm1Var!=1) {
        tmvaReader_BDTs[i]->AddVariable("H_mass", &H_mass);
      }
      //if (Nm1Var!=2) {
      //  tmvaReader_BDTs[i]->AddVariable("HJ1_HJ2_dR", &HJ1_HJ2_dR);
      //}
      if (Nm1Var!=3) {
        tmvaReader_BDTs[i]->AddVariable("H_pt",   &H_pt);
      }
      if (Nm1Var!=4) {
        tmvaReader_BDTs[i]->AddVariable("V_pt",   &V_pt);
      }
      if (Nm1Var!=5) {
        tmvaReader_BDTs[i]->AddVariable("Jet_btagCSV[hJetInd2]", &hJets_btagCSV_1);
      }
      if (Nm1Var!=6) {
        tmvaReader_BDTs[i]->AddVariable("Top1_mass_fromLepton_regPT_w4MET", &Top1_mass_fromLepton_regPT_w4MET);
        //tmvaReader_BDTs[i]->AddVariable("Top1_mass_fromLepton_regPT_w4MET"", &Top1_mass_fromLepton_regPT_w4MET");
      }
      if (Nm1Var!=7) {
        tmvaReader_BDTs[i]->AddVariable("HVdPhi", &HVdPhi);
      }
      if (Nm1Var!=8) {
        tmvaReader_BDTs[i]->AddVariable("nAddJet_f", &nAddJet_f);
      }
      if (Nm1Var!=9) {
        tmvaReader_BDTs[i]->AddVariable("lepMetDPhi", &lepMetDPhi);
      }
      if (Nm1Var!=10) {
        tmvaReader_BDTs[i]->AddVariable("softActivityVH_njets5", &softActivityVH_njets5_f);
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
      if (Nm1Var!=21) {
        tmvaReader_BDTs[i]->AddVariable("V_mt", &V_mt);
      }
      if (Nm1Var!=22) {
        tmvaReader_BDTs[i]->AddVariable("met_pt", &met_pt);
      }
      //tmvaReader_BDTs[i]->AddVariable("HVdEta_4MET", &HVdEta_4MET);
      //tmvaReader_BDTs[i]->AddVariable("JJEtaBal", &JJEtaBal);
      tmvaReader_BDTs[i]->BookMVA("BDT", weights[i]);
      cout<<"finished tmvaReader_BDT "<<i<<endl;
  }

  TFile *oldfile = new TFile(oldfilename);
  
  TTree *oldtree = (TTree*)oldfile->Get("tree");

  oldtree->SetBranchAddress("H_mass", &H_mass);
  oldtree->SetBranchAddress("HJ1_HJ2_dR", &HJ1_HJ2_dR);
  oldtree->SetBranchAddress("H_pt", &H_pt);
  oldtree->SetBranchAddress("V_pt", &V_pt);
  oldtree->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);
  oldtree->SetBranchAddress("Top1_mass_fromLepton_regPT_w4MET", &Top1_mass_fromLepton_regPT_w4MET);
  oldtree->SetBranchAddress("HVdPhi", &HVdPhi);
  oldtree->SetBranchAddress("nAddJet_f", &nAddJet_f);
  oldtree->SetBranchAddress("lepMetDPhi", &lepMetDPhi);
  oldtree->SetBranchAddress("softActivityVH_njets5", &softActivityVH_njets5);
  oldtree->SetBranchAddress("Jet_pt_reg", Jet_pt_reg);
  oldtree->SetBranchAddress("selLeptons_relIso_0", &selLeptons_relIso_0);
  oldtree->SetBranchAddress("selLeptons_eta", selLeptons_eta);
  oldtree->SetBranchAddress("selLeptons_pt", &selLeptons_pt);
  oldtree->SetBranchAddress("Jet_eta", Jet_eta);
  oldtree->SetBranchAddress("jjWPtBalance", &jjWPtBalance);
  oldtree->SetBranchAddress("AddJets252p9_puid_leadJet_btagCSV", &AddJets252p9_puid_leadJet_btagCSV);
  //oldtree->SetBranchAddress("AddJets252p9_puid_leadJet_pt", &AddJets252p9_puid_leadJet_pt);
  oldtree->SetBranchAddress("V_mt", &V_mt);
  oldtree->SetBranchAddress("met_pt", &met_pt);
  oldtree->SetBranchAddress("hJetInd1", &hJetInd1);
  oldtree->SetBranchAddress("hJetInd2", &hJetInd2);
  oldtree->SetBranchAddress("lepInd", &lepInd);
  oldtree->SetBranchAddress("HVdEta", &HVdEta);
  oldtree->SetBranchAddress("JJEtaBal", &JJEtaBal);
  Long64_t nentries = oldtree->GetEntries();

  TFile *newfile = new TFile(newfilename, "recreate");
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
    oldtree->GetEntry(z);
    //isWenu_f = (float) isWenu;
    hJets_btagCSV_0 = Jet_btagCSV[hJetInd1];
    hJets_btagCSV_1 = Jet_btagCSV[hJetInd2];
    hJets_pt_0 = Jet_pt_reg[hJetInd1];
    hJets_pt_1 = Jet_pt_reg[hJetInd2];
    hJets_eta_0 = Jet_eta[hJetInd1];
    hJets_eta_1 = Jet_eta[hJetInd2];
    selLeptons_pt_0 = selLeptons_pt[lepInd];
    selLeptons_eta_0 = selLeptons_eta[lepInd];
    softActivityVH_njets5_f = (float) softActivityVH_njets5;
    
    //BDT  = tmvaReader_BDT->EvaluateMVA("Gradient");
    for (int i=0; i<nBDTs; i++) {       
        BDT[i]  = tmvaReader_BDTs[i]->EvaluateMVA("BDT");       

        if(z%50000==1) cout<<z<<"  "<<BDT[i]<<std::endl;
    }
    /*std::cout<<"H_mass = "<<H_mass<<std::endl;
    std::cout<<"HJ1_HJ2_dR = "<<HJ1_HJ2_dR<<std::endl;
    std::cout<<"H_pt = "<<H_pt<<std::endl;
    std::cout<<"V_pt = "<<V_pt<<std::endl;
    std::cout<<"hJets_btagCSV_0 = "<<hJets_btagCSV_0<<std::endl;
    std::cout<<"hJets_btagCSV_1 = "<<hJets_btagCSV_1<<std::endl;
    std::cout<<"Top1_mass_fromLepton_regPT_w4MET" = "<<Top1_mass_fromLepton_regPT_w4MET<<std::endl;
    std::cout<<"HVdPhi = "<<HVdPhi<<std::endl;
    std::cout<<"nAddJet_f = "<<nAddJet_f<<std::endl;
    std::cout<<"lepMetDPhi = "<<lepMetDPhi<<std::endl;
    std::cout<<"softActivityVH_nJets5 = "<<softActivityVH_njets5<<std::endl;
    std::cout<<"hJets_pt_0 = "<<hJets_pt_0<<std::endl;
    std::cout<<"hJets_pt_1 = "<<hJets_pt_1<<std::endl;
    std::cout<<"selLeptons_relIso_0 = "<<selLeptons_relIso_0<<std::endl;
    std::cout<<"selLeptons_eta_0 = "<<selLeptons_eta_0<<std::endl;
    std::cout<<"selLeptons_pt_0 = "<<selLeptons_pt_0<<std::endl;
    std::cout<<"hJets_eta_0 = "<<hJets_eta_0<<std::endl;
    std::cout<<"hJets_eta_1 = "<<hJets_eta_1<<std::endl;
    std::cout<<"jjWPtBalance = "<<jjWPtBalance<<std::endl;
    std::cout<<"AddJets252p9_puid_leadJet_btagCSV = "<<AddJets252p9_puid_leadJet_btagCSV<<std::endl;
    std::cout<<"V_mt = "<<V_mt<<std::endl;
    std::cout<<"met_pt = "<<met_pt<<std::endl;
    std::cout<<"hJetInd1 = "<<hJetInd1<<std::endl;
    std::cout<<"hJetInd2 = "<<hJetInd2<<std::endl;
    std::cout<<"lepInd = "<<lepInd<<std::endl;
    */    
    newtree->Fill();
  }

  //newtree->Print();
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
