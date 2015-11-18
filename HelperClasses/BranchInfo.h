#ifndef branchInfo_h
#define branchInfo_h

#include <map>

class BranchInfo{

public:
    std::string name;
    int type;           // uint, int, float, double, bool, etc
    int length;         // for arrays only    
    std::string prov;   // origin of branch:  existing or new
    float val;          // for settings tree... val is static
    bool onlyMC;        // only get branch for MC

    BranchInfo(std::string _name, int _type, int _length=-1, int _onlyMC=0, std::string _prov="existing", float _val=-999);

};

inline BranchInfo::BranchInfo(std::string _name, int _type, int _length, int _onlyMC, std::string _prov, float _val) :
    name(_name),
    type(_type),
    length(_length),
    prov(_prov),
    val(_val),
    onlyMC(_onlyMC)
{
}

//BranchInfo::BranchInfo(std::string _name, int _type, int _length, std::string _prov){
//    BranchInfo(_name, _type, _length, _prov, -999.);
//}
//
//BranchInfo::BranchInfo(std::string _name, int _type, int _length){
//    BranchInfo(_name, _type, _length, "existing");
//}
//
//BranchInfo::BranchInfo(std::string _name, int _type, std::string _prov){
//    BranchInfo(_name, _type, -1, _prov, -999);
//}
//
//BranchInfo::BranchInfo(std::string _name, int _type){
//    BranchInfo(_name, _type, -1, "existing");
//}
//
//BranchInfo::BranchInfo() {
//}


#endif
