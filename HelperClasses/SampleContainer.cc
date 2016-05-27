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
    sampleChain = new TChain("tree");
    xsec = 0;
    kFactor= 1;
    scale = 1;
    files.clear();
    processedEvents=0;
    intWeight=1;
    nProFromFile=false;
    doJetFlavorSplit = false;
    procEff = 1.;
}

inline void SampleContainer::AddFile(const char* fname,int isBatch) {
    files.push_back(fname);
    std::cout<<"in SC "<<fname<<std::endl;
    std::cout<<"isBatch "<<isBatch<<std::endl;
    
    if( isBatch==1 ) return;
    
    sampleChain->Add(fname);
    //std::cout<<nProFromFile<<std::endl; 
    if(nProFromFile) {
        TFile *file = TFile::Open(fname);
        //TH1F* counter = (TH1F*)file.Get("Count");
        //processedEvents+=counter->GetBinContent(1);
        TH1F* counterPos = (TH1F*)file->Get("CountPosWeight");
        TH1F* counterNeg = (TH1F*)file->Get("CountNegWeight");
        int nEffective = counterPos->GetBinContent(1) - counterNeg->GetBinContent(1); 
        if(sampleNum==49) {
            //special prescription for WJets_BGenFilter sample weighting
            TH1F* counterFullWeight = (TH1F*)file->Get("CountFullWeighted");
            nEffective = counterFullWeight->GetBinContent(1); 
        }
        processedEvents += nEffective;
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
