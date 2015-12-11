#ifndef BDTVariable_h
#define BDTVariable_h

#include<iostream>

class BDTVariable{

public:
  std::string varName;
  std::string localVarName;
  //int varType; 
  bool isSpec;
  bool isExisting; 
  
  BDTVariable(std::string, std::string, bool, bool);
  BDTVariable();
};

#endif
