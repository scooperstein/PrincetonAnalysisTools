#ifndef branchInfo_h
#define branchInfo_h

#include <map>

class BranchInfo{

public:
    std::string name;
    int type;                 // uint, int, float, double, bool, etc
    int length;               // for arrays only    
    bool onlyMC;              // only get branch for MC
    std::string prov;         // origin of branch:  existing or new
    float val;                // for settings tree... val is static
    std::string lengthBranch; // name of the branch e.g. nJets

    BranchInfo(std::string _name, int _type, int _length=-1, int _onlyMC=0, std::string _prov="existing", float _val=-999, std::string _lengthBranch="");

};

inline BranchInfo::BranchInfo(std::string _name, int _type, int _length, int _onlyMC, std::string _prov, float _val, std::string _lengthBranch) :
    name(_name),
    type(_type),
    length(_length),
    onlyMC(_onlyMC),
    prov(_prov),
    val(_val),
    lengthBranch(_lengthBranch)
{
}


#endif
