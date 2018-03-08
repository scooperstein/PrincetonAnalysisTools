#include "BDTInfo.h"

BDTInfo::BDTInfo(std::string _methodName, std::string _bdtname, std::string _xmlFile) {
    bdtname = _bdtname;
    methodName = _methodName;
    bdtVars = std::vector<BDTVariable>();
    //method = _method;
    xmlFile = _xmlFile;
    reader = new TMVA::Reader( "!Color:!Silent:Error" );    
    reader->SetVerbose(kTRUE);
}

BDTInfo::BDTInfo(BDTInfo& _bdtInfo) {
    bdtname = _bdtInfo.bdtname;
    methodName = _bdtInfo.methodName;
    bdtVars = std::vector<BDTVariable>();
    for(unsigned int iVar=0; iVar<_bdtInfo.bdtVars.size(); iVar++){
        bdtVars.push_back(_bdtInfo.bdtVars[iVar]);
    }
    xmlFile = _bdtInfo.xmlFile;
    reader = new TMVA::Reader( "!Color:!Silent:Error" );    
    reader->SetVerbose(kTRUE);
}

BDTInfo::BDTInfo() {
    BDTInfo("", "", "");
}

void BDTInfo::AddVariable(std::string varName, std::string localVarName, bool isExisting, bool isSpec) {
    bdtVars.push_back(BDTVariable(varName, localVarName, isExisting, isSpec));
}

void BDTInfo::BookMVA(){
    reader->BookMVA(methodName,xmlFile);
}
