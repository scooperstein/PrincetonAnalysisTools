#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_data.root -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_mc.root -o output_mc.root -w '(1./CS_SF)*weight*(2.28/9.24)'
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_ttpowheg.root -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w '(1./CS_SF)*weight*(2.28/9.24)'
python plotFromTree.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_data.root -I inputfiles_onlydata.dat -o output_data.root
python plotFromTree.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_mc.root -o output_mc.root -w '(bTagWeightICHEP/bTagWeight)*(1.0/CS_SF)*weight*(12.89/12.89)'
python plotFromTree.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_ttpowheg.root -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w '(bTagWeightICHEP/bTagWeight)*(1.0/CS_SF)*weight*(12.89/12.89)'
python plotFromTree.py -i /eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V24_Wlnu_SRVV_Nov23/output_signal.root -I inputfiles_onlysignal.dat -o output_signal.root -w '(bTagWeightICHEP/bTagWeight)*(1.0/CS_SF)*weight*(12.89/12.89)'
hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root
