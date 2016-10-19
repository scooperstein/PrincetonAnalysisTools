python ../../splitSamples.py -i ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_mc.root -c WmnHighPt -p "(sampleIndex!=0||run<=276811)&&selLeptons_relIso_0<0.06&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=10&&isWmunu&&V_pt>100&&Vtype==2" -s ../../systematics_Wmn_CS.txt -w 2.0,CS_SF_new2 -t ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_data.root -v CMS_vhbb_BDT_Wln_13TeV -tt ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_ttpowheg.root  -d 1 -n 15 -o $ORIG_DIR/hists_WmnHighPt.root -b $ORIG_DIR/binStats_WmnHighPt.txt
python ../../splitSamples.py -i ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_mc.root -c WenHighPt -p "selLeptons_pt[lepInd]>30&&(sampleIndex!=0||run<=276811)&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=10&&isWenu&&V_pt>100&&Vtype==3" -s ../../systematics_Wen_CS.txt -w 2.0,CS_SF_new2 -t ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_data.root -v CMS_vhbb_BDT_Wln_13TeV -tt ../../../VHbbAnalysis/V24_Wln_SR_Oct3_forPreapproval_v2/output_ttpowheg.root  -d 1 -n 15 -o $ORIG_DIR/hists_WenHighPt.root -b $ORIG_DIR/binStats_WenHighPt.txt
