#ifndef BDTInfo_h
#define BDTInfo_h

#include "TMVA/Reader.h"
#include "BDTVariable.h"

class BDTInfo{

public:
  std::string bdtname;
  std::string bdtmethod;
  std::vector<BDTVariable> bdtVars;
  //std::string method;
  std::string xmlFile;
  TMVA::Reader *reader;

  BDTInfo(std::string, std::string, std::string);
  BDTInfo();
  void AddVariable(std::string, std::string, bool, bool);

};

#endif
