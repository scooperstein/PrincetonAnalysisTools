#cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/StatTools/Jan21
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#cmsenv
#for bdtvar in noHVdPhi
#for bdtvar in top14Var noJet2CSV noMbb noDRjj noHPt noWPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noLepPt noLepEta noMET noNSA5Jet noPtBal 
#for bdtvar in top16Var noJet2CSV noMbb noDRjj noHPt noWPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noLepPt noLepEta noAddJetCSV noWMt noMET noNSA5Jet noPtBal 
#for bdtvar in noWMt noMET noNSA5Jet 
#for bdtvar in top10
#for bdtvar in noJet2Pt noLepPt noLepIso noJet1CSV noLepEta noJet1Eta noJet2Eta noPtBal noAddJetCSV noWMt noMET noNSA5Jet 
#for bdtvar in noJet2CSV noMbb noDRjj noHPt noWPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noJet1Pt noJet2Pt noLepPt noJet1CSV noLepEta noJet1Eta noJet2Eta noPtBal noAddJetCSV noWMt noMET noNSA5Jet 
#for bdtvar in allVar noJet2CSV noMbb noDRjj noHPt noWPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noJet1Pt noJet2Pt noLepPt noLepIso noJet1CSV noLepEta noJet1Eta noJet2Eta noPtBal noAddJetCSV noWMt noMET noNSA5Jet 
#for bdtvar in allVar noJet2CSV noMbb noHPt noVPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noMET noNSA5Jet 
#for bdtvar in noJet2CSV noMbb noHPt noVPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noWMt noMET noNSA5Jet 
#for bdtvar in noLepPt
#for bdtvar in allVar noJet2CSV noMbb 
#for bdtvar in noHPt noVPt noTopMass
#for bdtvar in noHVdPhi noNAddJet noLepMetDPhi
#for bdtvar in noLepPt noMET noNSA5Jet 
#for bdtvar in allVar noJet2CSV noHPt noAddJetCSV noWMt noMET
#for bdtvar in allVar noAddJetCSV
#for bdtvar in allVar noLepIso
#for bdtvar in allVar top10Var noJet2CSV noMbb noDRjj noHPt noWPt noTopMass noHVdPhi noNAddJet noLepMetDPhi noJet1Pt noJet2Pt noLepPt noLepIso noJet1CSV noLepEta noJet1Eta noJet2Eta noPtBal noAddJetCSV noWMt noMET noNSA5Jet 
#for bdtvar in 200_3 200_4 200_5
#for bdtvar in 150_3 150_4 150_5 200_3 200_4 200_5
#for bdtvar in 250_3 250_4 250_5 300_3 300_4 300_5
#for bdtvar in 350_3 350_4 350_5 400_3 400_4 400_5
#for bdtvar in 450_3 450_4 450_5 500_3 500_4 500_5
#for bdtvar in 550_3 550_4 550_5 600_3 600_4 600_5
#for bdtvar in 650_3 650_4 650_5 700_3 700_4 700_5
#for bdtvar in 750_3 750_4 750_5 800_3 800_4 800_5
#for bdtvar in 150_3 150_4 150_5 200_3 200_4 200_5 250_3 250_4 250_5 300_3 300_4 300_5 350_3 350_4 350_5 400_3 400_4 400_5 450_3 450_4 450_5 500_3 500_4 500_5 550_3 550_4 550_5 600_3 600_4 600_5 650_3 650_4 650_5 700_3 700_4 700_5 750_3 750_4 750_5 800_3 800_4 800_5
#for bdtvar in 150_3 200_3 250_3 300_3 350_3 400_3 450_3 500_3 550_3 600_3 650_3 700_3 750_3 800_3 150_4 200_4 250_4 300_4 350_4 400_4 450_4 500_4 550_4 600_4 650_4 700_4 750_4 800_4 150_5 200_5 250_5 300_5 350_5 400_5 450_5 500_5 550_5 600_5 650_5 700_5 750_5 800_5
for bdtvar in 400_5
#for bdtvar in 150_3 200_3 250_3 300_3 350_3 400_3 450_3 500_3 550_3 600_3 650_3 700_3 750_3
#for bdtvar in 300_3
#for bdtvar in $1
#for vptcut in 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180 185 190 195 200 205 210 215 220 225 230 235 240 245 250
#for vptcut in 255 260 265 270 275 280 285 290 295 300
#for bdtvar in wFullTopMass wMWj
#for bdtvar in 450_3 wHVdEta wJJEtaBal
#for bdtvar in wHVdEta_4MET wJJEtaBal_and_HVdEta_4MET
#for bdtvar in TMBestCSV
do 
    echo $bdtvar
    #echo $vptcut
    mkdir -p V24_Sep20_fullTTStats_${bdtvar}
    cd V24_Sep20_fullTTStats_${bdtvar}
    #mkdir -p May6_VPtSplit_$vptcut
    #cd May6_VPtSplit_$vptcut
    #cd 10Nm1_${bdtvar}_fixedBinning
    #bash ../runStats.sh ../../../VHbbAnalysis/V21_Wlnu_April6_SR_forOpt/output_mc_wMay2BDTs.root blah ../../../VHbbAnalysis/V21_Wlnu_April6_SR_forOpt/output_data_wMay2BDTs.root blah 1 180 BDT_May3_$bdtvar
    ## an example call for signal region datacards
    #bash ../runStats.sh ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_mc.root blah ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_data.root blah 1 100 BDT_V24_Sep20_fullTTStats_$bdtvar ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_ttpowheg.root
    ## an example call for control region datacards
    #bash ../runStats.sh bsl ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_mc.root blah dsf ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_data.root 1 100 BDT_V24_Sep20_$bdtvar ../../../VHbbAnalysis/V24_Wln_SR_Sep15/output_ttpowheg.root
    #bash ../runStats.sh ../../../VHbbAnalysis/V21_Wlnu_May30_SR_splitWJetSFs/output_mc.root blah ../../../VHbbAnalysis/V21_Wlnu_May30_SR_splitWJetSFs/output_data.root blah 1 $vptcut BDT_May6_450_3
    bash ../test.sh
    #cd ..
    #cp vhbb_W*.txt hists_W*.root binStats_W*.txt April28_$bdtvar
    #
    # Calc limits
    #cd April28_$bdtvar
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/vhbb_tt*.txt .
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/hists_tt*.root .
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/vhbb_wlf*.txt .
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/hists_wlf*.root .
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/vhbb_whf*.txt .
    #cp ../Sep15_CR_noSF_12p9fb_oldHLTEleEffs_v2/hists_whf*.root .
    combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_Wln_13TeV.txt
    #combineCards.py WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_SR_13TeV.txt 
    #combineCards.py WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_SR_13TeV.txt 
    combine -M ProfileLikelihood --significance vhbb_Wln_13TeV.txt -t -1 --expectSignal=1| grep Significance
    #combine -M ProfileLikelihood --significance vhbb_WmnLowPt_13TeV.txt -t -1 --expectSignal=1 | grep Significance
    #combine -M ProfileLikelihood --significance vhbb_WenLowPt_13TeV.txt -t -1 --expectSignal=1 | grep Significance
    #combine -M ProfileLikelihood --significance vhbb_WmnHighPt_13TeV.txt -t -1 --expectSignal=1 | grep Significance
    #combine -M ProfileLikelihood --significance vhbb_WenHighPt_13TeV.txt -t -1 --expectSignal=1 | grep Significance
    #combine -M ProfileLikelihood --significance vhbb_SR_13TeV.txt -t -1 --expectSignal=1 -S 0 | grep Significance
    #combine -M MaxLikelihoodFit --expectSignal=1 -t -1 vhbb_Wln_13TeV.txt --justFit | grep "Best fit r"
    #combine -M Asymptotic --run expected -v 2 --X-rtd ADDNLL_RECURSIVE=0 --rMin=-4 --rMax=12  -m 125 vhbb_SR_13TeV.txt | grep "Expected 50.0%"
    #combine -M Asymptotic -t -1 vhbb_Wln_13TeV.txt
    cd ..
done
