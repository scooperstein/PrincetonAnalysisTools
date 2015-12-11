#include "BDTInfo.h"

BDTInfo::BDTInfo(std::string _bdtmethod, std::string _bdtname, std::string _xmlFile) {
    bdtname = _bdtname;
    bdtmethod = _bdtmethod;
    bdtVars = std::vector<BDTVariable>();
    //method = _method;
    xmlFile = _xmlFile;
    reader = new TMVA::Reader( "!Color:Silent" );    

}

BDTInfo::BDTInfo() {
    BDTInfo("", "", "");
}

void BDTInfo::AddVariable(std::string varName, std::string localVarName, bool isExisting, bool isSpec) {
    bdtVars.push_back(BDTVariable(varName, localVarName, isExisting, isSpec));
}


