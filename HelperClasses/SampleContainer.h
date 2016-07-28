// code template from https://github.com/h2gglobe/h2gglobe

#ifndef SAMPLECONTAINER
#define SAMPLECONTAINER

#include <string>
#include <map>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

class SampleContainer {
	
//    static float defaultextw;
  
 public:	
    //SampleContainer(const float * extw=0);
    SampleContainer();
    ~SampleContainer();
  
    TChain *sampleChain;
    std::vector<std::string> files;
    std::string sampleName;
    int sampleNum;
    float xsec;
    float kFactor; 
    float scale;
    float processedEvents;
    float intWeight;
    bool nProFromFile;
    bool doJetFlavorSplit; // split events by jet parton flavor
    float procEff; // process a fraction of the total sample events   
    TH1F *CountWeightedLHEWeightScale;
    TH1F *CountWeightedLHEWeightPdf;
    TH1F *CountWeighted;
    int lepFlav; // if 0 only select Wmunu events, if 1 only select Wenu events
 
    void AddFile(const char* fname, int isBatch=0);
    void ComputeWeight(float);

  /** adds a lumi section range to 'goodLumis' (typically
      taken from a 'json' file containing the list of 
      certified luminosity sections).
      @param run the run to be added 
      @param lumi1 the first lumi section to be added
      @param lumi2 the last lumi secetion to be added */
  //void addGoodLumi(int run, int lumi1, int lumi2 );

  //void addEventToList(int run, int lumi, int event );
  
  //bool isdata() const { return itype == 0; };
  //float weight() const { return ( (extweight!=0 && *extweight > 0 && ! isdata()) ? (*extweight)*intWeight : intWeight); };
 
 
  /*int itype;
  int ind;
  int histoplotit;
  std::string filesshortnam;
  long long int ntot;
  int nred;
  float lumi; 
  int forceVersion;
  float lumireal;
  bool hasLumiSelection, hasEventList;
  std::map<int, std::vector<std::pair<int,int> > > goodLumis;
  std::map<int, std::vector<std::pair<int,int> > > eventList;
  */

  //std::string pileup;
  
 //private:
  //const float * extweight;
  //float intWeight;


};

#endif
