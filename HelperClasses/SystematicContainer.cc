// code template from https://github.com/h2gglobe/h2gglobe

#include "SystematicContainer.h"
#include <utility>
#include <iostream>

inline SystematicContainer::~SystematicContainer() 
{}

inline SystematicContainer::SystematicContainer() 
{
    name = "";
    apply=true;
    scaleVar = "";
}

inline void SystematicContainer::TurnOff(){
    apply=false;
}

inline void SystematicContainer::AddBranchName(std::string brnchName){
    branchesToEdit.push_back(brnchName);
}

inline void SystematicContainer::AddScale(float newScale){
    scales.push_back(newScale);
}

inline void SystematicContainer::AddSmear(float newSmear){
    smears.push_back(newSmear);
}

