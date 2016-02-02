#include "SFContainer.h"
#include <utility>
#include <iostream>

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
    //std::cout<<"called getScaleFactor"<<std::endl;
    //std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
    int biny = scaleMap->GetYaxis()->FindBin(eta);
    //std::cout<<binx<<": ,"<<biny<<std::endl;
    if (binx != 0 && binx != scaleMap->GetNbinsX()+1 && biny != 0 && biny != scaleMap->GetNbinsY()+1) {
        sfactor = scaleMap->GetBinContent(binx, biny);
        sf_err = scaleMap->GetBinError(binx, biny);
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}
