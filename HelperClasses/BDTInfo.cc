#include "BDTInfo.h"

BDTInfo::BDTInfo(std::string _bdtname, std::string _xmlFile) {
    bdtname = _bdtname;
    inputNames = std::vector<std::string>();
    localVarNames = std::vector<std::string>();
    //method = _method;
    xmlFile = _xmlFile;
    reader = new TMVA::Reader( "!Color:Silent" );    

}

BDTInfo::BDTInfo() {
    BDTInfo("", "");
}

void BDTInfo::AddVariable(std::string varName, std::string localVarName) {
    inputNames.push_back(varName);
    localVarNames.push_back(localVarName);
}

void BDTInfo::AddSpectatorVariable(std::string varName, std::string localVarName) {
    inputSpectatorNames.push_back(varName);
    localSpectatorVarNames.push_back(localVarName);
}

