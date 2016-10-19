#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_mc.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -o output_mc.root -w 'weight*(1.0/1.0)*(12.89/12.89)**(2.28/9.24)'
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July1.0/output_ttpowheg.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w 'weight*(1.0/1.0)*(12.89/12.89)**(2.28/9.24)'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_CR_Sep23_noMW/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_CR_Sep23_noMW/output_mc.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -o output_mc.root -w 'weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))*(bTagWeight/bTagWeight)'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_CR_Sep23_noMW/output_ttpowheg.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w 'weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))*(bTagWeight/bTagWeight)'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_CR_Sep23_noMW/output_signal.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlysignal.dat -o output_signal.root -w 'weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))*(bTagWeight/bTagWeight)'
hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root
#python plotFromTree.py -i ../../V24_Wln_CR_Sep23_noMW/output_data.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree.py -i ../../V24_Wln_CR_Sep23_noMW/output_mc.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -o output_mc.root 
#python plotFromTree.py -i ../../V24_Wln_CR_Sep23_noMW/output_ttpowheg.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlyttpowheg.dat -o output_ttpow.root 
#python plotFromTree.py -i ../../V24_Wln_CR_Sep23_noMW/output_signal.root -P plotvariables_CS.dat -C categories_CS.dat -S selection_CS.dat -I inputfiles_onlysignal.dat -o output_signal.root
#hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root
