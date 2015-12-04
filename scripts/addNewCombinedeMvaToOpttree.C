void addNewCombinedeMvaToOpttree(const char* oldfilename="testingtree_mh125_evenbkg_newmvas.root",const char* newfilename="testingtree_mh125_evenbkg_allmvas.root", const char* weights="weights/TMVA_vbf_dijet_diphosherpa_BDT_evenb_fixwt_Gradient.weights.xml",const char* newbranch="combined_bdt",bool refill=false) {
  
#include <string>


  using std::max;
  using std::min;
  
  TMVA::Reader  *tmvaReader_BDT; 

  Float_t   H_mass, H_pt, V_pt, hJets_btagCSV_0, hJets_btagCSV_1, HVdPhi, nAddJet_f, nAddLep_f, hJets_pt_0, hJets_pt_1, met_pt, Top1_mass_fromLepton, lepMetDPhi, selLeptons_pt_0, selLeptons_eta_0, selLeptons_relIso03_0, isWenu_f, isWmunu_f;

  Int_t isWenu;
  Float_t BDT;

  tmvaReader_BDT = new TMVA::Reader("V:!Color:Silent");
  //tmvaReader_BDT->AddVariable("H_mass",&H_mass);
  tmvaReader_BDT->AddVariable("H_pt", &H_pt);
  tmvaReader_BDT->AddVariable("V_pt", &V_pt);
  tmvaReader_BDT->AddVariable("Jet_btagCSV[hJetInd1]", &hJets_btagCSV_0);
  tmvaReader_BDT->AddVariable("Jet_btagCSV[hJetInd2]", &hJets_btagCSV_1);
  tmvaReader_BDT->AddVariable("HVdPhi", &HVdPhi);
  tmvaReader_BDT->AddVariable("nAddJet_f", &nAddJet_f);
  tmvaReader_BDT->AddVariable("nAddLep_f", &nAddLep_f);
  tmvaReader_BDT->AddVariable("Jet_pt[hJetInd1]", &hJets_pt_0);
  tmvaReader_BDT->AddVariable("Jet_pt[hJetInd2]", &hJets_pt_1);
  tmvaReader_BDT->AddVariable("met_pt", &met_pt);
  tmvaReader_BDT->AddVariable("Top1_mass_fromLepton", &Top1_mass_fromLepton);
  tmvaReader_BDT->AddVariable("lepMetDPhi", &lepMetDPhi);
  tmvaReader_BDT->AddVariable("selLeptons_pt[lepInd]", &selLeptons_pt_0);
  tmvaReader_BDT->AddVariable("selLeptons_eta[lepInd]", &selLeptons_eta_0);
  tmvaReader_BDT->AddVariable("selLeptons_relIso03[lepInd]", &selLeptons_relIso03_0);
  tmvaReader_BDT->AddVariable("isWenu", &isWenu_f);
  tmvaReader_BDT->AddVariable("isWmunu_f", &isWmunu_f);

  tmvaReader_BDT->BookMVA("BDT", weights);
  cout<<"finished tmvaReader_BDT"<<endl;

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

  TBranch *bBDT=0;
  if(refill){
    oldtree->SetBranchAddress(newbranch,  &BDT);
  } else {  
    char nametype[100];
    strcpy(nametype,newbranch);
    strcat(nametype,"/F");
    
    bBDT = newtree->Branch(newbranch,&BDT);
  }

  cout<<"starting loop over "<<nentries<<" events"<<endl;
  for (Long64_t z=0; z<nentries; z++) {
    oldtree->GetEntry(z);
    isWenu_f = (float) isWenu;

    //BDT  = tmvaReader_BDT->EvaluateMVA("Gradient");       
    BDT  = tmvaReader_BDT->EvaluateMVA("BDT");       

    if(z%50000==1) cout<<z<<"  "<<BDT<<std::endl;
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
