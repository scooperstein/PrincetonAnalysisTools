#include "BDTVariable.h"

BDTVariable::BDTVariable(std::string varName_, std::string localVarName_, bool isExisting_, bool isSpec_) {
    varName = varName_;
    localVarName = localVarName_;
    isExisting = isExisting_;
    //varType = varType_;
    isSpec = isSpec_;
}

// default constructer for initialization
BDTVariable::BDTVariable() {
    BDTVariable("", "",false, false); 
} 
