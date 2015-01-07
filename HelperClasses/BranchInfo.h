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

    BranchInfo(std::string _name, int _type, int _length=-1, std::string _prov="existing", float _val=-999);

};

BranchInfo::BranchInfo(std::string _name, int _type, int _length, std::string _prov, float _val) :
    name(_name),
    type(_type),
    length(_length),
    prov(_prov),
    val(_val)
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
