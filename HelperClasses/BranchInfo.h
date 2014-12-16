#ifndef branchInfo_h
#define branchInfo_h

#include <map>

class BranchInfo{

public:
    std::string name;
    int type;   // uint, int, float, double, bool, etc
    int length; // for arrays only    
    std::string prov; // origin of branch:  existing or new

    BranchInfo(std::string _name, int _type, int _length, std::string _prov);
    BranchInfo(std::string _name, int _type, std::string _prov);
    BranchInfo(std::string _name, int _type, int _length);
    BranchInfo(std::string _name, int _type);
    BranchInfo();

};

BranchInfo::BranchInfo(std::string _name, int _type, int _length, std::string _prov){
    name=_name;
    type=_type;
    length=_length;
    prov=_prov;
}

BranchInfo::BranchInfo(std::string _name, int _type, int _length){
    BranchInfo(_name, _type, _length, "existing");
}

BranchInfo::BranchInfo(std::string _name, int _type, std::string _prov){
    BranchInfo(_name, _type, -1, _prov);
}

BranchInfo::BranchInfo(std::string _name, int _type){
    BranchInfo(_name, _type, -1, "existing");
}

BranchInfo::BranchInfo() {
}


#endif
