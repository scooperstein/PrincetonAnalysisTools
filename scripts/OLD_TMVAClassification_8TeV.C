// @(#)root/tmva $Id: TMVAClassification.C 37399 2010-12-08 15:22:07Z evt $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void OLD_TMVAClassification_8TeV(int cutsType = 1, int cutsType2 = 0, int sigType = 3, int WjetsType = 1, int TTbarType = 1, int doVV = 0, int nTrees = 1000, int nEventsMin = 150, int maxDepth = 3, float adaboostbeta = 0.2, int iteration = 0, TString myMethodList = "")
{

  std::cout<<"\n--------------------------------------------------\n"
	   <<"This iteration is with the following configuration:\n"
 	   <<"nTrees = "<<nTrees<<"\nnEventsMin = "<<nEventsMin
 	   <<"\nmaxDepth = "<<maxDepth<<"\nadaBoostBeta = "<<adaboostbeta
 	   <<"\n--------------------------------------------------\n";

   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0; // 1
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0; // 1
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0; // 1
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // 1 k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // 1 Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;   //1
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins
   char* CutsName = "";
   if(cutsType == 0)
     CutsName = "noMjj";
   else if(cutsType == 1)
     CutsName = "withMjj";
   else if(cutsType == 2)
     CutsName = "newCuts";
   else if(cutsType == 3)
     CutsName = "newCuts2";
   else if(cutsType == 4)
     CutsName = "newCuts3";
   else if(cutsType == 5)
     CutsName = "newCuts4";
   else if(cutsType == 6)
     CutsName = "newCuts5";
   else if(cutsType == 7)
     CutsName = "newCuts6";

   char* CutsName2 = "";
   if(cutsType2 == 0)
     CutsName2 = "";
   else if(cutsType2 == 1)
     CutsName2 = "LowPt";
   else if(cutsType2 == 2)
     CutsName2 = "HighPt";

   char* sigName = "";
   if(sigType == -1)
     sigName = "VV";
   else if(sigType == 0)
     sigName = "H110";
   else if(sigType == 1)
     sigName = "H115";
   else if(sigType == 2)
     sigName = "H120";
   else if(sigType == 3)
     sigName = "H125";
   else if(sigType == 4)
     sigName = "H130";
   else if(sigType == 5)
     sigName = "H135";
   else if(sigType == 6)
     sigName = "H140";
   else if(sigType == 7)
     sigName = "H145";
   else if(sigType == 8)
     sigName = "H150";

   char* WjetsName = "";
   if(WjetsType == 1)
     WjetsName = "0b1b2bWjets";
   else if(WjetsType == 2)
     WjetsName = "LFHFWjets";

   char* TTbarName = "";
   if(TTbarType == 1)
     TTbarName = "NewTTbar";
   else if(TTbarType == 2)
     TTbarName = "OldTTbar";

   char* VVName = "";
   if(doVV == 1)
     VVName = "VV";

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   char outfileName[80];
   if(iteration==0) sprintf(outfileName,"TMVA_8TeV_%sSig_%s%s%sBkg_%s%s.root",sigName,WjetsName,TTbarName,VVName,CutsName,CutsName2);
   else sprintf(outfileName,"TMVA_8TeV_%sSig_%s%s%sBkg_%s%s_%i.root",sigName,WjetsName,TTbarName,VVName,CutsName,CutsName2,iteration);
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   char factoryName[80];
   if(iteration==0) sprintf(factoryName,"TMVA_8TeV_%sSig_%s%s%sBkg_%s%s",sigName,WjetsName,TTbarName,VVName,CutsName,CutsName2);
   else sprintf(factoryName,"TMVA_8TeV_%sSig_%s%s%sBkg_%s%s_%i",sigName,WjetsName,TTbarName,VVName,CutsName,CutsName2,iteration);
   TMVA::Factory *factory = new TMVA::Factory( factoryName, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   if(cutsType >= 1)
     factory->AddVariable( "H.massCorr", "H.massCorr", "", 'F' );
   factory->AddVariable( "H.ptCorr", "H.ptCorr", "", 'F' );
   factory->AddVariable( "V.pt", "V.pt", "", 'F' );
///   factory->AddVariable( "maxbtag:=max(hJet_csv[0],hJet_csv[1])", "max(hJet_csv[0],hJet_csv[1])", "", 'F' );
///   factory->AddVariable( "minbtag:=min(hJet_csv[0],hJet_csv[1])", "min(hJet_csv[0],hJet_csv[1])", "", 'F' );
   factory->AddVariable( "hJet_csvReshapedNew[0]", "hJet_csvReshapedNew[0]", "", 'F' );
   factory->AddVariable( "hJet_csvReshapedNew[1]", "hJet_csvReshapedNew[1]", "", 'F' );
   factory->AddVariable( "HVdPhi", "HVdPhi", "", 'F' );
   factory->AddVariable( "H.dEta", "H.dEta", "", 'F');
////   factory->AddVariable( "NAddJet:=Sum$(aJet_pt>20 && abs(aJet_eta)<2.5)", "Sum$(aJet_pt>20 && abs(aJet_eta)<2.5)", "", 'I' );
   factory->AddVariable( "NAddJet:=Sum$(aJet_pt>20 && abs(aJet_eta)<4.5)", "Sum$(aJet_pt>20 && abs(aJet_eta)<4.5)", "", 'I' );
   factory->AddVariable( "H.dR", "H.dR", "", 'F');

//   factory->AddVariable( "hJet_cosTheta[0]", "hJet_cosTheta[0]", "", 'F');
//   factory->AddVariable( "hJet_cosTheta[1]", "hJet_cosTheta[1]", "", 'F');
   factory->AddVariable( "PullAngle:=abs(deltaPullAngle)", "abs(deltaPullAngle)", "", 'F');

   factory->AddVariable( "hJet_ptCorr[0]", "hJet_ptCorr[0]", "", 'F');
   factory->AddVariable( "hJet_ptCorr[1]", "hJet_ptCorr[1]", "", 'F');
   if(cutsType == 2)
   {
     //factory->AddVariable( "vLepton_charge[0]", "vLepton_charge[0]", "", 'I');
     factory->AddVariable( "TopM_regres", "TopM_regres", "", 'F');
     factory->AddVariable( "TopPt_regres", "TopPt_regres", "", 'F');
   }
   else if(cutsType == 3)
   {
     factory->AddVariable( "NAddJetFwd:=Sum$(aJet_pt>20 && abs(aJet_eta)>2.5 && abs(aJet_eta)<4.5)", "Sum$(aJet_pt>20 && abs(aJet_eta)>2.5 && abs(aJet_eta)<4.5)", "", 'I');
   }
   else if(cutsType == 4)
   {
     //factory->AddVariable( "NAddJetBtag:=Sum$(aJet_pt>20 && abs(aJet_eta)<2.5 && aJet_csvReshapedNew>0.4)", "Sum$(aJet_pt>20 && abs(aJet_eta)<2.5 && aJet_csvReshapedNew>0.4)", "", 'I');
     //factory->AddVariable( "mTWH", "mTWH", "", 'F');
     //factory->AddVariable( "pTWH", "pTWH", "", 'F');
     //factory->AddVariable("mTWH := evalMTWH(H.massCorr, H.ptCorr, H.phi, 80.385, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","evalMTWH(H.massCorr, H.ptCorr, H.phi, 80.385, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","",'F');
     //factory->AddVariable("pTWH := evalPTWH(H.ptCorr, H.phi, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","evalPTWH(H.ptCorr, H.phi, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","",'F');
   }
   else if(cutsType == 5)
   {
     //factory->AddVariable( "FatHmass:=FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "FatH.filteredmass * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "", 'F');
     //factory->AddVariable( "FatHpt:=FatH.filteredpt * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "FatH.filteredpt * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "", 'F');
     //factory->AddVariable( "FatHdeltaM:=(H.massCorr - FatH.filteredmass) * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "(H.massCorr - FatH.filteredmass) * (nfathFilterJets>0 && FatH.FatHiggsFlag==1)", "", 'F');
     factory->AddVariable("maxAddCSV := MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)", "MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)", "", 'F');
     factory->AddVariable("mindRAddJetH := MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)", "MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)", "", 'F');
   }
   else if(cutsType == 6)
   {
     //factory->AddVariable("cosThetaHbb := abs(evalCosThetaHbb(hJet_ptCorr[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptCorr[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", "abs(evalCosThetaHbb(hJet_ptCorr[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptCorr[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", "", 'F');
     factory->AddVariable("cosThetaHbb := abs(hJet_cosTheta[0])","abs(hJet_cosTheta[0])","",'F');
   }
   else if(cutsType == 7)
   {
     factory->AddVariable( "TopM_regres", "TopM_regres", "", 'F');
     factory->AddVariable( "TopPt_regres", "TopPt_regres", "", 'F');
     factory->AddVariable( "NAddJetFwd:=Sum$(aJet_pt>20 && abs(aJet_eta)>2.5 && abs(aJet_eta)<4.5)", "Sum$(aJet_pt>20 && abs(aJet_eta)>2.5 && abs(aJet_eta)<4.5)", "", 'I');
     //factory->AddVariable( "mTWH", "mTWH", "", 'F');
     //factory->AddVariable( "pTWH", "pTWH", "", 'F');
     //factory->AddVariable("mTWH := evalMTWH(H.massCorr, H.ptCorr, H.phi, 80.385, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","evalMTWH(H.massCorr, H.ptCorr, H.phi, 80.385, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","",'F');
     //factory->AddVariable("pTWH := evalPTWH(H.ptCorr, H.phi, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","evalPTWH(H.ptCorr, H.phi, vLepton_pt[0], vLepton_phi[0], METtype1corr.et, METtype1corr.phi)","",'F');
     factory->AddVariable("maxAddCSV := MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)", "MaxIf$(aJet_csvReshapedNew,aJet_pt>20 && abs(aJet_eta)<2.5)", "", 'F');
     factory->AddVariable("mindRAddJetH := MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)", "MinIf$(evalDeltaR(aJet_eta,aJet_phi,H.eta,H.phi),aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_puJetIdL>0)", "", 'F');
     //factory->AddVariable("cosThetaHbb := abs(evalCosThetaHbb(hJet_ptCorr[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptCorr[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", "abs(evalCosThetaHbb(hJet_ptCorr[0], hJet_pt[0], hJet_eta[0], hJet_phi[0], hJet_e[0], hJet_ptCorr[1], hJet_pt[1], hJet_eta[1], hJet_phi[1], hJet_e[1]))", "", 'F');
     factory->AddVariable("cosThetaHbb := abs(hJet_cosTheta[0])","abs(hJet_cosTheta[0])","",'F');
   }

   factory->AddSpectator( "Vtype",  "Vtype", "", 'I' );
   factory->AddSpectator( "Addleptons:=Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )",  "Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )", "", 'I' );
   factory->AddSpectator( "METtype1corr.et", "METtype1corr.et", "", 'F' );
//   factory->AddSpectator( "hJet_pt[0]", "hJet_pt[0]", "", 'F');
//   factory->AddSpectator( "hJet_pt[1]", "hJet_pt[1]", "", 'F');


   char *higgsName = (char*) "125";
   if(sigType == 0)
     higgsName = (char*) "110";
   else if(sigType == 1)
     higgsName = (char*) "115";
   else if(sigType == 2)
     higgsName = (char*) "120";
   else if(sigType == 3)
     higgsName = (char*) "125";
   else if(sigType == 4)
     higgsName = (char*) "130";
   else if(sigType == 5)
     higgsName = (char*) "135";
   else if(sigType == 6)
     higgsName = (char*) "140";
   else if(sigType == 7)
     higgsName = (char*) "145";
   else if(sigType == 8)
     higgsName = (char*) "150";

   TFile fSig_Training(Form("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_WH_%s_lumiWeighted.root",higgsName));
   TTree *sig_Training = (TTree*)fSig_Training.Get("tree");

   TFile fBkgWjetLF_100to180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_oddEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjetLF_100to180_Training = (TTree*) (fBkgWjetLF_100to180_Training.Get("tree")); 
   TFile fBkgWjetHF_100to180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_oddEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjetHF_100to180_Training = (TTree*) (fBkgWjetHF_100to180_Training.Get("tree"));
   TFile fBkgWjet0b_100to180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/0b_oddEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet0b_100to180_Training = (TTree*) (fBkgWjet0b_100to180_Training.Get("tree"));
   TFile fBkgWjet1b_100to180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/1b_oddEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet1b_100to180_Training = (TTree*) (fBkgWjet1b_100to180_Training.Get("tree"));
   TFile fBkgWjet2b_100to180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/2b_oddEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet2b_100to180_Training = (TTree*) (fBkgWjet2b_100to180_Training.Get("tree"));

   TFile fBkgWjetLF_180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_oddEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjetLF_180_Training = (TTree*) (fBkgWjetLF_180_Training.Get("tree")); 
   TFile fBkgWjetHF_180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_oddEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjetHF_180_Training = (TTree*) (fBkgWjetHF_180_Training.Get("tree"));
   TFile fBkgWjet0b_180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/0b_oddEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet0b_180_Training = (TTree*) (fBkgWjet0b_180_Training.Get("tree"));
   TFile fBkgWjet1b_180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/1b_oddEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet1b_180_Training = (TTree*) (fBkgWjet1b_180_Training.Get("tree"));
   TFile fBkgWjet2b_180_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/2b_oddEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet2b_180_Training = (TTree*) (fBkgWjet2b_180_Training.Get("tree"));
   
   TFile fBkgTTbar_OLD_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/TTbar_old_lumiWeighted.root");
   TTree* bkgTTbar_OLD_Training = (TTree*) (fBkgTTbar_OLD_Training.Get("tree"));
   /*TFile fBkgTTbar_SemiLept_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/oddEvents_DiJetPt_TTJets2_SemiLeptMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_SemiLept_Training = (TTree*) (fBkgTTbar_SemiLept_Training.Get("tree"));
   TFile fBkgTTbar_FullLept_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/oddEvents_DiJetPt_TTJets2_FullLeptMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_FullLept_Training = (TTree*) (fBkgTTbar_FullLept_Training.Get("tree"));
   TFile fBkgTTbar_Hadronic_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/oddEvents_DiJetPt_TTJets2_HadronicMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_Hadronic_Training = (TTree*) (fBkgTTbar_Hadronic_Training.Get("tree"));
   */
   TFile fBkgTTbar_SemiLept_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_TTbar_SemiLept_lumiWeighted.root");
   TTree* bkgTTbar_SemiLept_Training = (TTree*) (fBkgTTbar_SemiLept_Training.Get("tree"));
   TFile fBkgTTbar_FullLept_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_TTbar_FullLept_lumiWeighted.root");
   TTree* bkgTTbar_FullLept_Training = (TTree*) (fBkgTTbar_FullLept_Training.Get("tree"));
   TFile fBkgTTbar_Hadronic_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_TTbar_Hadronic_lumiWeighted.root");
   TTree* bkgTTbar_Hadronic_Training = (TTree*) (fBkgTTbar_Hadronic_Training.Get("tree"));

   TFile fBkgWW_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_WW_lumiWeighted.root");
   TTree* bkgWW_Training = (TTree*) (fBkgWW_Training.Get("tree"));
   TFile fBkgWWLF_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_oddEvents_WW_lumiWeighted.root");
   TTree* bkgWWLF_Training = (TTree*) (fBkgWWLF_Training.Get("tree"));
   TFile fBkgWWHF_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_oddEvents_WW_lumiWeighted.root");
   TTree* bkgWWHF_Training = (TTree*) (fBkgWWHF_Training.Get("tree"));
   TFile fBkgWZ_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/oddEvents_WZ_lumiWeighted.root");
   TTree* bkgWZ_Training = (TTree*) (fBkgWZ_Training.Get("tree"));
   TFile fBkgWZLF_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_oddEvents_WZ_lumiWeighted.root");
   TTree* bkgWZLF_Training = (TTree*) (fBkgWZLF_Training.Get("tree"));
   TFile fBkgWZHF_Training("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_oddEvents_WZ_lumiWeighted.root");
   TTree* bkgWZHF_Training = (TTree*) (fBkgWZHF_Training.Get("tree"));


   TFile fSig_Test(Form("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_WH_%s_lumiWeighted.root",higgsName));
   TTree *sig_Test = (TTree*)fSig_Test.Get("tree");

   TFile fBkgWjetLF_100to180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_evenEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjetLF_100to180_Test = (TTree*) (fBkgWjetLF_100to180_Test.Get("tree")); 
   TFile fBkgWjetHF_100to180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_evenEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjetHF_100to180_Test = (TTree*) (fBkgWjetHF_100to180_Test.Get("tree"));
   TFile fBkgWjet0b_100to180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/0b_evenEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet0b_100to180_Test = (TTree*) (fBkgWjet0b_100to180_Test.Get("tree"));
   TFile fBkgWjet1b_100to180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/1b_evenEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet1b_100to180_Test = (TTree*) (fBkgWjet1b_100to180_Test.Get("tree"));
   TFile fBkgWjet2b_100to180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/2b_evenEvents_WJetsPt100To180Madgraph_lumiWeighted.root");
   TTree* bkgWjet2b_100to180_Test = (TTree*) (fBkgWjet2b_100to180_Test.Get("tree"));

   TFile fBkgWjetLF_180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_evenEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjetLF_180_Test = (TTree*) (fBkgWjetLF_180_Test.Get("tree")); 
   TFile fBkgWjetHF_180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_evenEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjetHF_180_Test = (TTree*) (fBkgWjetHF_180_Test.Get("tree"));
   TFile fBkgWjet0b_180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/0b_evenEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet0b_180_Test = (TTree*) (fBkgWjet0b_180_Test.Get("tree"));
   TFile fBkgWjet1b_180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/1b_evenEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet1b_180_Test = (TTree*) (fBkgWjet1b_180_Test.Get("tree"));
   TFile fBkgWjet2b_180_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/2b_evenEvents_WJetsPt180Madgraph_lumiWeighted.root");
   TTree* bkgWjet2b_180_Test = (TTree*) (fBkgWjet2b_180_Test.Get("tree"));
   
   TFile fBkgTTbar_OLD_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/TTbar_old_lumiWeighted.root");
   TTree* bkgTTbar_OLD_Test = (TTree*) (fBkgTTbar_OLD_Test.Get("tree"));
   /*TFile fBkgTTbar_SemiLept_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/evenEvents_DiJetPt_TTJets2_SemiLeptMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_SemiLept_Test = (TTree*) (fBkgTTbar_SemiLept_Test.Get("tree"));
   TFile fBkgTTbar_FullLept_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/evenEvents_DiJetPt_TTJets2_FullLeptMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_FullLept_Test = (TTree*) (fBkgTTbar_FullLept_Test.Get("tree"));
   TFile fBkgTTbar_Hadronic_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/evenEvents_DiJetPt_TTJets2_HadronicMGDecays_8TeV-madgraph-part_lumiWeight.root");
   TTree* bkgTTbar_Hadronic_Test = (TTree*) (fBkgTTbar_Hadronic_Test.Get("tree"));
   */
   TFile fBkgTTbar_SemiLept_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_TTbar_SemiLept_lumiWeighted.root");
   TTree* bkgTTbar_SemiLept_Test = (TTree*) (fBkgTTbar_SemiLept_Test.Get("tree"));
   TFile fBkgTTbar_FullLept_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_TTbar_FullLept_lumiWeighted.root");
   TTree* bkgTTbar_FullLept_Test = (TTree*) (fBkgTTbar_FullLept_Test.Get("tree"));
   TFile fBkgTTbar_Hadronic_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_TTbar_Hadronic_lumiWeighted.root");
   TTree* bkgTTbar_Hadronic_Test = (TTree*) (fBkgTTbar_Hadronic_Test.Get("tree"));

   TFile fBkgWW_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_WW_lumiWeighted.root");
   TTree* bkgWW_Test = (TTree*) (fBkgWW_Test.Get("tree"));
   TFile fBkgWWLF_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_evenEvents_WW_lumiWeighted.root");
   TTree* bkgWWLF_Test = (TTree*) (fBkgWWLF_Test.Get("tree"));
   TFile fBkgWWHF_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_evenEvents_WW_lumiWeighted.root");
   TTree* bkgWWHF_Test = (TTree*) (fBkgWWHF_Test.Get("tree"));
   TFile fBkgWZ_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v16/nominal/evenEvents_WZ_lumiWeighted.root");
   TTree* bkgWZ_Test = (TTree*) (fBkgWZ_Test.Get("tree"));
   TFile fBkgWZLF_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/LF_evenEvents_WZ_lumiWeighted.root");
   TTree* bkgWZLF_Test = (TTree*) (fBkgWZLF_Test.Get("tree"));
   TFile fBkgWZHF_Test("/eos/uscms/store/user/sbc01/v16_transfer/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_genSplit_varsAddedSummed_v16/HF_evenEvents_WZ_lumiWeighted.root");
   TTree* bkgWZHF_Test = (TTree*) (fBkgWZHF_Test.Get("tree"));


   Double_t sigWeight = 2.0 * 20.0 * (1.0)/(7684.8); // 125 GeV
   if(sigType == 0)
     sigWeight = 2.0 * 20.0 * (1.0)/(3877.2); // 110 GeV
   else if(sigType == 1)
     sigWeight = 2.0 * 20.0 * (1.0)/(4788.9); // 115 GeV
   else if(sigType == 2)
     sigWeight = 2.0 * 20.0 * (1.0)/(5986.3); // 120 GeV
   else if(sigType == 3)
     sigWeight = 2.0 * 20.0 * (1.0)/(7684.8); // 125 GeV
   else if(sigType == 4)
     sigWeight = 2.0 * 20.0 * (1.0)/(10285.6); // 130 GeV
   else if(sigType == 5)
     sigWeight = 2.0 * 20.0 * (1.0)/(14323.0); // 135 GeV
   else if(sigType == 6)
     sigWeight = 2.0 * 20.0 * (1.0)/(16203.8); // 140 GeV
   else if(sigType == 7)
     sigWeight = 2.0 * 20.0 * (1.0)/(18360.7); // 145 GeV
   else if(sigType == 8)
     sigWeight = 2.0 * 20.0 * (1.0)/(20700.3); // 150 GeV
     
   Double_t bkgWeightWjetLF_100to180 = 2.0 * 20.0 * (0.95)/(42.73);
   Double_t bkgWeightWjetHF_100to180 = 2.0 * 20.0 * (1.9)/(42.73);
   Double_t bkgWeightWjet0b_100to180 = 2.0 * 20.0 * (0.95)/(42.73);
   Double_t bkgWeightWjet1b_100to180 = 2.0 * 20.0 * (3.0)/(42.73);
   Double_t bkgWeightWjet2b_100to180 = 2.0 * 20.0 * (0.95)/(42.73);
   Double_t bkgWeightWjetLF_180 = 2.0 * 20.0 * (0.95)/(265.2);
   Double_t bkgWeightWjetHF_180 = 2.0 * 20.0 * (1.9)/(265.2);
   Double_t bkgWeightWjet0b_180 = 2.0 * 20.0 * (0.95)/(265.2);
   Double_t bkgWeightWjet1b_180 = 2.0 * 20.0 * (3.0)/(265.2);
   Double_t bkgWeightWjet2b_180 = 2.0 * 20.0 * (0.95)/(265.2);
   Double_t bkgWeightTTbar_OLD = 1.0 * 20.0 * (1.0)/(112.8);
   Double_t bkgWeightTTbar_SemiLept = 2.0 * 20.0 * (0.95)/(226.9);
   Double_t bkgWeightTTbar_FullLept = 2.0 * 20.0 * (0.95)/(257.1);
   Double_t bkgWeightTTbar_Hadronic = 2.0 * 20.0 * (0.95)/(206.4);
   Double_t bkgWeightWW = 2.0 * 20.0 * (1.0)/(176.0);
   Double_t bkgWeightWZ = 2.0 * 20.0 * (1.0)/(295.2);
   
   if(sigType >= 0)
   {
     factory->AddSignalTree(sig_Training,sigWeight,TMVA::Types::kTraining);
     factory->AddSignalTree(sig_Test,sigWeight,TMVA::Types::kTesting);
   }
   else if(sigType == -1)
   {
     factory->AddSignalTree(bkgWW_Training,bkgWeightWW,TMVA::Types::kTraining);
     factory->AddSignalTree(bkgWZ_Training,bkgWeightWZ,TMVA::Types::kTraining);
     factory->AddSignalTree(bkgWW_Test,bkgWeightWW,TMVA::Types::kTesting);
     factory->AddSignalTree(bkgWZ_Test,bkgWeightWZ,TMVA::Types::kTesting);
   }
   
   if(WjetsType == 1)
   {
     factory->AddBackgroundTree(bkgWjet0b_100to180_Training,bkgWeightWjet0b_100to180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet1b_100to180_Training,bkgWeightWjet1b_100to180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet2b_100to180_Training,bkgWeightWjet2b_100to180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet0b_180_Training,bkgWeightWjet0b_180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet1b_180_Training,bkgWeightWjet1b_180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet2b_180_Training,bkgWeightWjet2b_180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjet0b_100to180_Test,bkgWeightWjet0b_100to180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjet1b_100to180_Test,bkgWeightWjet1b_100to180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjet2b_100to180_Test,bkgWeightWjet2b_100to180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjet0b_180_Test,bkgWeightWjet0b_180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjet1b_180_Test,bkgWeightWjet1b_180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjet2b_180_Test,bkgWeightWjet2b_180,TMVA::Types::kTesting);

     //factory->AddBackgroundTree(bkgWWLF_Training,bkgWeightWW,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWZLF_Training,bkgWeightWZ,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWWLF_Test,bkgWeightWW,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWZLF_Test,bkgWeightWZ,TMVA::Types::kTesting);
   }
   else if(WjetsType == 2)
   {
     factory->AddBackgroundTree(bkgWjetLF_100to180_Training,bkgWeightWjetLF_100to180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjetHF_100to180_Training,bkgWeightWjetHF_100to180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjetLF_180_Training,bkgWeightWjetLF_180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjetHF_180_Training,bkgWeightWjetHF_180,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWjetLF_100to180_Test,bkgWeightWjetLF_100to180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjetHF_100to180_Test,bkgWeightWjetHF_100to180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjetLF_180_Test,bkgWeightWjetLF_180,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWjetHF_180_Test,bkgWeightWjetHF_180,TMVA::Types::kTesting);

     //factory->AddBackgroundTree(bkgWWLF_Training,bkgWeightWW,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWZLF_Training,bkgWeightWZ,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWWLF_Test,bkgWeightWW,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWZLF_Test,bkgWeightWZ,TMVA::Types::kTesting);
   }
 
   if(TTbarType == 1)
   {
     factory->AddBackgroundTree(bkgTTbar_SemiLept_Training,bkgWeightTTbar_SemiLept,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgTTbar_FullLept_Training,bkgWeightTTbar_FullLept,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgTTbar_Hadronic_Training,bkgWeightTTbar_Hadronic,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgTTbar_SemiLept_Test,bkgWeightTTbar_SemiLept,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgTTbar_FullLept_Test,bkgWeightTTbar_FullLept,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgTTbar_Hadronic_Test,bkgWeightTTbar_Hadronic,TMVA::Types::kTesting);
   }
   else if(TTbarType == 2)
   {
     factory->AddBackgroundTree(bkgTTbar_OLD_Training,bkgWeightTTbar_OLD,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgTTbar_SemiLept_Test,bkgWeightTTbar_SemiLept,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgTTbar_FullLept_Test,bkgWeightTTbar_FullLept,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgTTbar_Hadronic_Test,bkgWeightTTbar_Hadronic,TMVA::Types::kTesting);
   }

   if(doVV == 1)
   {
     factory->AddBackgroundTree(bkgWW_Training,bkgWeightWW,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWZ_Training,bkgWeightWZ,TMVA::Types::kTraining);
     factory->AddBackgroundTree(bkgWW_Test,bkgWeightWW,TMVA::Types::kTesting);
     factory->AddBackgroundTree(bkgWZ_Test,bkgWeightWZ,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWWHF_Training,bkgWeightWW,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWZHF_Training,bkgWeightWZ,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWWHF_Test,bkgWeightWW,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWZHF_Test,bkgWeightWZ,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWjetHF_100to180_Training,bkgWeightWjetHF_100to180,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWjetHF_180_Training,bkgWeightWjetHF_180,TMVA::Types::kTraining);
     //factory->AddBackgroundTree(bkgWjetHF_100to180_Test,bkgWeightWjetHF_100to180,TMVA::Types::kTesting);
     //factory->AddBackgroundTree(bkgWjetHF_180_Test,bkgWeightWjetHF_180,TMVA::Types::kTesting);
   }

   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
///   factory->SetBackgroundWeightExpression( "weightMC_dij" );

// try e+mu

   TCut mycutsSig;
   TCut mycutsBkg;

   if (cutsType2 == 0)
   {
     mycutsSig = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>120. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
     mycutsBkg = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>120. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
   }
   else if (cutsType2 == 1)
   {
     mycutsSig = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>120. && V.pt<170. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
     mycutsBkg = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>120. && V.pt<170. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
   }
   else if (cutsType2 == 2)
   {
     mycutsSig = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>170. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
     mycutsBkg = "H.massCorr>0. && H.massCorr<256. && H.ptCorr>120. && ((Vtype==2 && vLepton_pt[0]>20.) || (Vtype==3 && vLepton_pt[0]>30. && vLepton_wp80[0]>0 && METtype1corr.et>45.)) && V.pt>170. && Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && METtype1corr.et<1000. && V.pt<1000 && max(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4  && min(hJet_csvReshapedNew[0],hJet_csvReshapedNew[1])>0.4 && hJet_csvReshapedNew[0]>0. && hJet_csvReshapedNew[1]>0. && hJet_ptCorr[0]>30. && hJet_ptCorr[1]>30. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0 && hbhe && abs(deltaPullAngle)<5.";
   }

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
//       factory->PrepareTrainingAndTestTree( mycutsSig, mycutsBkg, "SplitMode=random:!V" );
//factory->PrepareTrainingAndTestTree(mycutsSig,mycutsBkg,"NormMode=None:!V");
  factory->PrepareTrainingAndTestTree(mycutsSig,mycutsBkg,"!V");
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
//   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=300:NEventsMax=800:MaxVIterations=150" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.1333:nActiveCells=50:nSampl=500:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   char parString[150];
   //sprintf(parString,"!H:!V:NTrees=%i:nEventsMin=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=MisClassificationError:nCuts=-1:PruneMethod=NoPruning"
   sprintf(parString,"!H:!V:NTrees=%i:nEventsMin=%i:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=%.2f:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning",nTrees,nEventsMin,maxDepth,adaboostbeta);
   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT", parString );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}

//inline double evalCosThetaHbb(double ptCor0, double pt0, double eta0, double phi0, double e0,double ptCor1, double pt1, double eta1, double phi1, double e1,bool rescaleEnergy=true)
//{
//  TLorentzVector j0, j1;
//  j0.SetPtEtaPhiE(ptCor0, eta0, phi0, (rescaleEnergy ? e0 * ptCor0 / pt0 : e0));
//  j1.SetPtEtaPhiE(ptCor1, eta1, phi1, (rescaleEnergy ? e1 * ptCor1 / pt1 : e1));
//  TLorentzVector jsum = j0+j1;
//  TVector3 boostsum = jsum.BoostVector();
//  j0.Boost(-boostsum);
//  j1.Boost(-boostsum);
//  TVector3 boost;
//  if(int(pt0) % 2 == 0)
//    boost = j0.BoostVector();
//  else
//    boost = j1.BoostVector();
//  
//  double cosTheta = boost.Dot(boostsum) / (boost.Mag() * boostsum.Mag());
//  return cosTheta;
//}

inline double evalMTWH(double m0, double pt0, double phi0, double m1, double pt1, double phi1, double pt2, double phi2)
{
  double ptW = TMath::Sqrt(pow(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2),2) + pow(pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2),2));
  double phiW = TMath::ATan((pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2))/(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2)));

  double dPhi = fabs(phi0-phiW);
  if(dPhi > 3.14159)
    dPhi = 2.0*3.14159 - dPhi;

  return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+ptW*ptW)) - pt0*pt1*TMath::Cos(dPhi)));
}

inline double evalPTWH(double pt0, double phi0, double pt1, double phi1, double pt2, double phi2)
{
  double ptW = TMath::Sqrt(pow(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2),2) + pow(pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2),2));
  double phiW = TMath::ATan((pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2))/(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2)));

  return TMath::Sqrt(pow(pt0*TMath::Cos(phi0)+ptW*TMath::Cos(phiW),2) + pow(pt0*TMath::Sin(phi0)+ptW*TMath::Sin(phiW),2));
}

inline double evalDeltaR(double eta0, double phi0, double eta1, double phi1)
{
  double dEta = fabs(eta0-eta1);
  double dPhi = fabs(phi0-phi1);
  
  if(dPhi > 3.14159)
    dPhi = 2.0*3.14159 - dPhi;

  return TMath::Sqrt(TMath::Power(dEta,2)+TMath::Power(dPhi,2));
}
