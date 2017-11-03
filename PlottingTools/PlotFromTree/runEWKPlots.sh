#python plotFromTree_2.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree_2.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_mc.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -o output_mc.root -w 'weight*(1.0/1.0)*(12.89/12.89)**(2.28/9.24)'
#python plotFromTree_2.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_ttpowheg.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w 'weight*(1.0/1.0)*(12.89/12.89)**(2.28/9.24)'
#python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_March13_data_reminiaod/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_March13_NLOWJets/output_wjetsnlo_2.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlynlowjets.dat -o output_nlowjets.root -w 'weight*bTagWeightMoriondCMVA*(1./CS_SF)'
python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10_data/output_data.root -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_April6_v2_data/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_CR_March13_TTMadNoSys/output_ttmadgraph_3.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttmad.dat -o output_ttpow.root -w 'weight*(1./CS_SF)*bTagWeightMoriondCMVA*VPtCorrFactorSplit3'
python plotFromTree_2.py -f /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk.dat -o output_ewk.root -w 'weight*corrQGL_norm*(1+(sampleIndex==-100))'
#python plotFromTree_2.py -f /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_EWK_SR_July5/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm*1.021*(1 + isWenu*(-1+0.9668))*(1 + isWmunu*(-1+0.9958))'
#python plotFromTree_2.py -f root://cmseos.fnal.gov://store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm*1.021*(1 + isWenu*(-1+0.982929*0.997822))*(1 + isWmunu*(-1+ 0.98666*0.998233))'
python plotFromTree_2.py -f root://cmseos.fnal.gov://store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm'
#python plotFromTree_2.py -f /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm*1.021'
#python plotFromTree_2.py -f /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_EWK_SR_May3/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk.dat -I inputfiles_ewk.dat -o output_ewk.root -w 'weight*(1./corrQGL)*(1./Lep_SF)'
#hadd -f output.root output_ewk.root output_data.root
hadd -f output.root output_ewk.root output_data.root output_wjets.root
cp output.root sr.root
#hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root
#hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root output_nlowjets.root
#python plotFromTree_2.py -i ../../V24_Wln_CR_Oct13_v2/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree_2.py -i ../../V24_Wln_CR_Oct13_v2/output_mc.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -o output_mc.root 
#python plotFromTree_2.py -i ../../V24_Wln_CR_Oct13_v2/output_ttpowheg.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttpowheg.dat -o output_ttpow.root 
#python plotFromTree_2.py -i ../../V24_Wln_CR_Oct13_v2/output_signal.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlysignal.dat -o output_signal.root
#hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root

python plotFromTree_2.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10_data/output_data.root -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk_cr.dat -I inputfiles_onlydata.dat -o output_data.root
python plotFromTree_2.py -f /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk_cr.dat -I inputfiles_ewk.dat -o output_ewk.root -w 'weight*corrQGL_norm*(1+(sampleIndex==-100))'
#python plotFromTree_2.py -f root://cmseos.fnal.gov://store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk_cr.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm*1.021*(1 + isWenu*(-1+0.982929*0.997822))*(1 + isWmunu*(-1+ 0.98666*0.998233))'
python plotFromTree_2.py -f root://cmseos.fnal.gov://store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct10/haddjobs/ -P plotvariables_ewk.dat -C categories_ewk.dat -S selection_ewk_cr.dat -I inputfiles_ewk_wjets.dat -o output_wjets.root -w 'weight*corrQGL_norm'
hadd -f output.root output_ewk.root output_data.root output_wjets.root
cp output.root cr.root
python addQCDFromData.py cr.root qcd.root
hadd -f output.root sr.root qcd.root
