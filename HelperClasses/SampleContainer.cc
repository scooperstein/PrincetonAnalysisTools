// code template from https://github.com/h2gglobe/h2gglobe

#include "SampleContainer.h"
#include <utility>
#include <iostream>

inline SampleContainer::~SampleContainer() 
{}

//float SampleContainer::defaultextw=1.;

//SampleContainer::SampleContainer(const float * extw) :
//	extweight(extw)
inline SampleContainer::SampleContainer() 
{
	/*intWeight = 1;
	itype = 0;
	ind = 0;
	histoplotit = 1;
	filesshortnam = "";
	ntot = 0;
	nred = 0;
	lumi = 0;
	lumireal = 1;
	hasLumiSelection = false;
	hasEventList = false;
	pileup = "";
	forceVersion = 0;*/
    sampleName = "";
    sampleNum = 0;
    sampleChain = new TChain("Events");
    //sampleChain = new TChain("tree");
    xsec = 0;
    kFactor= 1;
    scale = 1;
    files.clear();
    processedEvents=0.;
    intWeight=1;
    nProFromFile=false;
    doJetFlavorSplit = false;
    procEff = 1.;
    CountWeightedLHEWeightScale = new TH1F("CountWeightedLHEWeightScale","CountWeightedLHEWeightScale",6,-0.5,5.5);
    CountWeightedLHEWeightPdf = new TH1F("CountWeightedLHEWeightPdf","CountWeightedLHEWeightPdf",103,-0.5,102.5);
    CountWeighted = new TH1F("CountWeighted","CountWeighted",1,0.,2.0);
    CountFullWeighted = new TH1F("CountFullWeighted","CountFullWeighted",1,0.,2.0);
    lepFlav = -1;
}

inline void SampleContainer::AddFile(const char* fname,int isBatch, int doSkim) {
    files.push_back(fname);
    std::cout<<"in SC "<<fname<<std::endl;
    std::cout<<"isBatch "<<isBatch<<std::endl;
    
    if( isBatch==1 ) return;
    
    sampleChain->Add(fname);
    //std::cout<<nProFromFile<<std::endl; 
    if(nProFromFile) {
        if (doSkim == 0 && files.size() > 1) return; // skimmed files already have the summed count histograms
        TFile *file = TFile::Open(fname);
        if (file->IsZombie()) return;
        //TH1F* counter = (TH1F*)file.Get("Count");
        //processedEvents+=counter->GetBinContent(1);
        //TH1F* counterPos = (TH1F*)file->Get("CountPosWeight");
        //TH1F* counterNeg = (TH1F*)file->Get("CountNegWeight");
        /*TH1F* counter = (TH1F*)file->Get("CountWeighted");
        TH1F* counterFullWeight = (TH1F*)file->Get("CountFullWeighted");
        //int nEffective = counterPos->GetBinContent(1) - counterNeg->GetBinContent(1); 
        float nEffective = counter->GetBinContent(1);
        if(sampleNum==49 or sampleNum==491) {
            //special prescription for WJets_BGenFilter sample weighting
            //TH1F* counterFullWeight = (TH1F*)file->Get("CountFullWeighted");
            nEffective = counterFullWeight->GetBinContent(1); 
        }
        CountWeighted->Add(counter);
        CountFullWeighted->Add(counterFullWeight);
        std::cout<<"pe = "<<processedEvents<<std::endl;
        processedEvents += nEffective;
        std::cout<<"pe = "<<processedEvents<<std::endl;
        TH1F* CountWeightedLHEWeightScale_thisfile = (TH1F*)file->Get("CountWeightedLHEWeightScale");
        TH1F* CountWeightedLHEWeightPdf_thisfile = (TH1F*)file->Get("CountWeightedLHEWeightPdf");
        std::cout<<"lhe = "<<CountWeightedLHEWeightPdf->GetBinContent(1)<<std::endl;
        CountWeightedLHEWeightScale->Add(CountWeightedLHEWeightScale_thisfile);
        CountWeightedLHEWeightPdf->Add(CountWeightedLHEWeightPdf_thisfile);
        std::cout<<"lhe = "<<CountWeightedLHEWeightPdf->GetBinContent(1)<<std::endl;*/

        // totally different setup for grabbing event count in nanoAOD
        TTree *Runs = (TTree*) file->Get("Runs");
        Long64_t genEventCount = 0;
        Runs->SetBranchAddress("genEventCount",&genEventCount);
        Runs->GetEntry(0);
        std::cout<<fname<<" genEventCount: "<<genEventCount<<std::endl;
        CountWeighted->SetBinContent(1,CountWeighted->GetBinContent(1)+genEventCount);
        processedEvents += genEventCount;
        file->Close();
    }
    
}

inline void SampleContainer::ComputeWeight(float intL) {
    if(sampleNum==0) { //this is data
        intWeight = 1; 
    } else {
        //std::cout << "Computing Weight for type - " << sampleNum << ", Using " << processedEvents << " Processed Events" << std::endl;
        intWeight = (kFactor*scale*xsec*intL)/(processedEvents*procEff);
    }
}

/* 
// ----------------------------------------------------------------------------------------------------------------------
void SampleContainer::addGoodLumi(int run, int lumi1, int lumi2 )
{
	hasLumiSelection = true;
	goodLumis[run].push_back( std::make_pair(lumi1,lumi2) );
}


// ----------------------------------------------------------------------------------------------------------------------
void SampleContainer::addEventToList(int run, int lumi, int event )
{
	hasEventList = true;
	eventList[run].push_back( std::make_pair(lumi,event) );
}*/
