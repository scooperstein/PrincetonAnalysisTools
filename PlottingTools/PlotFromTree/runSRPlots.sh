#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July19/output_data.root -I inputfiles_onlydata.dat -o output_data.root
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July19/output_mc.root -o output_mc.root -w '(1./CS_SF)*weight*(2.28/9.24)'
#python plotFromTree.py -i ../../VHbbAnalysis/V23_Wlnu_CR_July19/output_ttpowheg.root -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w '(1./CS_SF)*weight*(2.28/9.24)'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_SR_Sep15/output_data.root -I inputfiles_onlydata.dat -o output_data.root
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_SR_Sep15/output_mc.root -o output_mc.root -w '(1./CS_SF)*weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_SR_Sep15/output_ttpowheg.root -I inputfiles_onlyttpowheg.dat -o output_ttpow.root -w '(1./CS_SF)*weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))'
python plotFromTree.py -i ../../VHbbAnalysis/V24_Wln_SR_Sep15/output_signal.root -I inputfiles_onlysignal.dat -o output_signal.root -w '(1./CS_SF)*weight*(12.89/12.89)*(1 + isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))'
hadd -f output.root output_mc.root output_data.root output_ttpow.root output_signal.root
