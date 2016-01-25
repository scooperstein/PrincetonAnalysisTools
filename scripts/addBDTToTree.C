void addBDTToTree(const char* oldfilename="testingtree_mh125_evenbkg_newmvas.root",const char* newfilename="testingtree_mh125_evenbkg_allmvas.root", std::vector<const char*> weights=std::vector<const char*>(),std::vector<const char*> newbranch=std::vector<const char*>(),bool refill=false) {
  
#include <string>


  using std::max;
  using std::min;
  
  int nBDTs = (int) weights.size();
  std::vector<TMVA::Reader*>  tmvaReader_BDTs; 

  Float_t   H_mass, H_pt, V_pt, hJets_btagCSV_0, hJets_btagCSV_1, HVdPhi, nAddJet_f, nAddLep_f, hJets_pt_0, hJets_pt_1, met_pt, Top1_mass_fromLepton, lepMetDPhi, selLeptons_pt_0, selLeptons_eta_0, selLeptons_relIso03_0, isWenu_f, isWmunu_f;

  Int_t isWenu;
  Float_t BDT[nBDTs];

  for (int i=0; i<nBDTs; i++) {
      tmvaReader_BDTs.push_back(new TMVA::Reader("V:!Color:Silent"));
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
      tmvaReader_BDTs[i]->AddVariable("isWmunu_f", &isWmunu_f);

      tmvaReader_BDTs[i]->BookMVA("BDT", weights[i]);
      cout<<"finished tmvaReader_BDT "<<i<<endl;
  }

  TFile *oldfile = new TFile(oldfilename);
  
  TTree *oldtree = (TTree*)oldfile->Get("tree");

  oldtree->SetBranchAddress("H_mass", &H_mass);
  oldtree->SetBranchAddress("H_pt", &H_pt);
  oldtree->SetBranchAddress("V_pt", &V_pt);
  oldtree->SetBranchAddress("hJets_btagCSV_0", &hJets_btagCSV_0);
  oldtree->SetBranchAddress("hJets_btagCSV_1", &hJets_btagCSV_1);
  oldtree->SetBranchAddress("HVdPhi", &HVdPhi);
  oldtree->SetBranchAddress("nAddJet_f", &nAddJet_f);
  oldtree->SetBranchAddress("nAddLep_f", &nAddLep_f);
  oldtree->SetBranchAddress("hJets_pt_0", &hJets_pt_0);
  oldtree->SetBranchAddress("hJets_pt_1", &hJets_pt_1);
  oldtree->SetBranchAddress("met_pt", &met_pt);
  oldtree->SetBranchAddress("Top1_mass_fromLepton", &Top1_mass_fromLepton);
  oldtree->SetBranchAddress("lepMetDPhi", &lepMetDPhi);
  oldtree->SetBranchAddress("selLeptons_pt_0", &selLeptons_pt_0);
  oldtree->SetBranchAddress("selLeptons_eta_0", &selLeptons_eta_0);
  oldtree->SetBranchAddress("selLeptons_relIso03_0", &selLeptons_relIso03_0);
  oldtree->SetBranchAddress("isWenu", &isWenu);
  oldtree->SetBranchAddress("isWmunu_f", &isWmunu_f);
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

  cout<<"starting loop over "<<nentries<<" events"<<endl;
  for (Long64_t z=0; z<nentries; z++) {
    oldtree->GetEntry(z);
    isWenu_f = (float) isWenu;

    //BDT  = tmvaReader_BDT->EvaluateMVA("Gradient");
    for (int i=0; i<nBDTs; i++) {       
        BDT[i]  = tmvaReader_BDTs[i]->EvaluateMVA("BDT");       

        if(z%50000==1) cout<<z<<"  "<<BDT[i]<<std::endl;
    }
    /*std::cout<<"H_mass = "<<H_mass<<std::endl;
    std::cout<<"H_pt = "<<H_pt<<std::endl;
    std::cout<<"V_pt = "<<V_pt<<std::endl;
    std::cout<<"hJets_btagCSV_0 = "<<hJets_btagCSV_0<<std::endl;
    std::cout<<"hJets_btagCSV_1 = "<<hJets_btagCSV_1<<std::endl;
    std::cout<<"HVdPhi = "<<HVdPhi<<std::endl;
    std::cout<<"nAddJet_f = "<<nAddJet_f<<std::endl;
    std::cout<<"nAddLep_f = "<<nAddLep_f<<std::endl;
    std::cout<<"hJets_pt_0 = "<<hJets_pt_0<<std::endl;
    std::cout<<"hJets_pt_1 = "<<hJets_pt_1<<std::endl;
    std::cout<<"met_pt = "<<met_pt<<std::endl;
    std::cout<<"Top1_mass_fromLepton = "<<Top1_mass_fromLepton<<std::endl;
    std::cout<<"lepMetDPhi = "<<lepMetDPhi<<std::endl;
    std::cout<<"selLeptons_pt_0 = "<<selLeptons_pt_0<<std::endl;
    std::cout<<"selLeptons_eta_0 = "<<selLeptons_eta_0<<std::endl;
    std::cout<<"selLeptons_relIso03_0 = "<<selLeptons_relIso03_0<<std::endl;
    std::cout<<"isWenu_f = "<<isWenu_f<<std::endl;
    std::cout<<"isWmunu_f = "<<isWmunu_f<<std::endl;
    */
    newtree->Fill();
  }

  newtree->Print();
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
