#ifndef BDTInfo_h
#define BDTInfo_h

#include "TMVA/Reader.h"

class BDTInfo{

public:
  std::string bdtname;
  std::vector<std::string> inputNames;
  std::vector<std::string> localVarNames;
  std::vector<std::string> inputSpectatorNames;
  std::vector<std::string> localSpectatorVarNames;
  //std::string method;
  std::string xmlFile;
  TMVA::Reader *reader;

  BDTInfo(std::string, std::string);
  BDTInfo();
  void AddVariable(std::string, std::string);
  void AddSpectatorVariable(std::string, std::string);

};

#endif
