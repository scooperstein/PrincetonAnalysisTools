//
//  Written by Stephane Cooperstein 
//
//  Analyzer for creating workspaces used in limits
//  Inheriting from Analysis Manager
//

#include <iostream>
#include "../AnalysisManager.h"

#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;


class WorkspaceAnalysis : public AnalysisManager {
    
    public:
    WorkspaceAnalysis();
    virtual ~WorkspaceAnalysis();

    virtual void InitAnalysis();
    virtual bool Preselection();
    virtual bool Analyze();
    virtual void FinishEvent();
    virtual void TermAnalysis();

    std::vector<std::string> catTypes;
    std::map<std::string, std::vector<TH1F*> > hists1D;
    std::map<std::string, std::vector<TH2F*> > hists2D;
    TFile *histout; // output root file for histograms

};

