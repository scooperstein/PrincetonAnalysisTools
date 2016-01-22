python ../splitSamples.py ../../PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_Jan8_wBDT.root WmnLowPt 'isWmunu&&V_pt>100&&V_pt<180'
python ../splitSamples.py ../../PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_Jan8_wBDT.root WmnHighPt 'isWmunu&&V_pt>180'
python ../splitSamples.py ../../PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_Jan8_wBDT.root WenLowPt 'isWenu&&V_pt>100&&V_pt<180'
python ../splitSamples.py ../../PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_Jan8_wBDT.root WenHighPt 'isWenu&&V_pt>=180'
python ../printYields.py vhbb_WmnHighPt_13TeV.txt WmnHighPt ../systematics_Wmn.txt binStats_WmnHighPt.txt
python ../printYields.py vhbb_WmnLowPt_13TeV.txt WmnLowPt ../systematics_Wmn.txt  binStats_WmnLowPt.txt
python ../printYields.py vhbb_WenHighPt_13TeV.txt WenHighPt ../systematics_Wen.txt binStats_WenHighPt.txt
python ../printYields.py vhbb_WenLowPt_13TeV.txt WenLowPt ../systematics_Wen.txt  binStats_WenLowPt.txt
combineCards.py WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt > vhbb_Wln_13TeV.txt
