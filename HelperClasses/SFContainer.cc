#include "SFContainer.h"
//#include <utility>
//#include <iostream>

inline SFContainer::~SFContainer() 
{}

inline SFContainer::SFContainer() 
{
    jsonFile = "";
    name = "";
    binning = "";
    branchname = "weight_SF";
    branches = std::vector<std::string>();
    scaleMap = new TH2F("scaleMap","scaleMap",10,0,500,10,-2.5,2.5);
}


inline void SFContainer::AddBranch(std::string branch_) {
    branches.push_back(branch_);
}

inline float SFContainer::getScaleFactor(float pt, float eta, float &sf_err) {
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
    int biny = scaleMap->GetYaxis()->FindBin(eta);
    if (binx != -1 && binx != scaleMap->GetNbinsX() && biny != -1 && biny != scaleMap->GetNbinsY()) {
        sfactor = scaleMap->GetBinContent(binx, biny);
        sf_err = scaleMap->GetBinError(binx, biny);
    }
    return sfactor;
}
