#ifndef SFCONTAINER
#define SFCONTAINER

#include <string>
#include <vector>
#include <TH2F.h>

class SFContainer {
	
 public:	
    SFContainer();
    ~SFContainer();
  
    std::string jsonFile; // file with scale factor inputs
    std::vector<std::string> branches; // branches needed as input
    TH2F *scaleMap; // scale factors binned in 2D (pT and eta for now)  
    std::string name;
    std::string binning;
    std::string branchname;

    void  AddBranch(std::string branch);
    float getScaleFactor(float pt, float eta, float &sf_err);
    
 //private:


};

#endif
