#!bin/bash
## Wrapper script to make datacards from skimmed ntuples
## $1 ntuple for signal region
## $2 ntuple for control regions
## $3 data ntuple for signal region
## $4 data ntuple for control regions
## $5 if 1 do electron and muon channel combined, if 0 do separately (default)
## $6 W pT split point between low and high pT channels
## $7 BDT name in ntuples
## $8 ntuple for ttbar powheg

ptSplit=180
if [ ! -z $6 ]; then
    ptSplit=$6
fi

bdtname="CMS_vhbb_BDT_Wln_13TeV"
if [ ! -z $7 ]; then
    bdtname=$7
fi

echo "ptSplit = "
echo "$ptSplit"
echo "bdtname = "
echo $bdtname

## FIXME take out application of CS_SF,8.26 weight, just using it for optimization studies. Keeping it in will ruin the MaxLikelihoodFits!!! Also set back when going back to scale factor fits
# Signal Region Workspaces
#python ../../splitSamples.py -i $1 -c WmnLowPt -p "run<=276811&&selLeptons_relIso_0<0.06&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWmunu&&V_pt>100&&V_pt<${ptSplit}&&Vtype==2" -s ../../systematics_Wmn_CS.txt -w 2.0 -t $3 -v $7 -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WmnHighPt -p "run<=276811&&selLeptons_relIso_0<0.06&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWmunu&&V_pt>${ptSplit}&&Vtype==2" -s ../../systematics_Wmn_CS.txt -w 2.0 -t $3 -v $7 -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WenLowPt -p "run<=276811&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWenu&&V_pt>100&&V_pt<${ptSplit}&&Vtype==3" -s ../../systematics_Wen_CS.txt -w 2.0 -t $3 -v $7 -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WenHighPt -p "run<=276811&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWenu&&V_pt>=${ptSplit}&&Vtype==3" -s ../../systematics_Wen_CS.txt -w 2.0 -t $3 -v $7 -tt $8 -d 1 -n 15
# Signal Region Workspaces
#python ../../splitSamples.py -i $1 -c WmnLowPt -p "run<=276811&&selLeptons_relIso_0<0.06&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWmunu&&V_pt>100&&V_pt<${ptSplit}&&Vtype==2" -s ../../systematics_Wmn_CS.txt -w 2.0 -t $3 -v BDT_V24_Sep20_fullTTStats_VPt100To150_500_4 -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WmnHighPt -p "run<=276811&&selLeptons_relIso_0<0.06&&(evt%2==0||sampleIndex==0)&&H_mass>60&&H_mass<120&&cutFlow>=9&&isWmunu&&V_pt>${ptSplit}&&Vtype==2" -s ../../systematics_Wmn_CS.txt -w 2.0 -t $3 -v BDT_V24_Sep20_fullTTStats_500_4_VVSig -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WenLowPt -p "run<=276811&&(evt%2==0||sampleIndex==0)&&H_mass>90&&H_mass<150&&cutFlow>=9&&isWenu&&V_pt>100&&V_pt<${ptSplit}&&Vtype==3" -s ../../systematics_Wen_CS.txt -w 2.0 -t $3 -v BDT_V24_Sep20_fullTTStats_VPt100To150_500_4 -tt $8 -d 1 -n 15
#python ../../splitSamples.py -i $1 -c WenHighPt -p "run<=276811&&(evt%2==0||sampleIndex==0)&&H_mass>60&&H_mass<120&&cutFlow>=9&&isWenu&&V_pt>=${ptSplit}&&Vtype==3" -s ../../systematics_Wen_CS.txt -w 2.0 -t $3 -v BDT_V24_Sep20_fullTTStats_500_4_VVSig -tt $8 -d 1 -n 15

### Control Regions Workspaces
#python ../../splitSamples.py -i $2 -c ttWmn -p "selLeptons_relIso_0<0.06&&run<=276811&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==1" -s ../../systematics_Wmn_CS.txt -v HJ1_HJ2_dR --xlow 0.0 --xhigh 3.5 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 
#python ../../splitSamples.py -i $2 -c ttWen -p "run<=276811&&isWenu&&V_pt>100&&Vtype==3&&controlSample==1" -s ../../systematics_Wen_CS.txt  -v HJ1_HJ2_dR --xlow 0.0 --xhigh 3.5 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 
#python ../../splitSamples.py -i $2 -c whfWmn -p "abs(HJ1_HJ2_dPhi)<1.0&&selLeptons_relIso_0<0.06&&run<=276811&&(H_mass<90||H_mass>150)&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wmn_CS.txt -v HJ1_HJ2_dR --xlow 0.0 --xhigh 3.5 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 
#python ../../splitSamples.py -i $2 -c whfWen -p "abs(HJ1_HJ2_dPhi)<1.0&&run<=276811&&(H_mass<90||H_mass>150)&&isWenu&&V_pt>100&&Vtype==3&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wen_CS.txt -v HJ1_HJ2_dR --xlow 0.0 --xhigh 3.5 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 
#python ../../splitSamples.py -i $2 -c wlfWmn -p "selLeptons_relIso_0<0.06&&run<=276811&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wmn_CS.txt -v V_pt --xlow 100 --xhigh 350 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 
#python ../../splitSamples.py -i $2 -c wlfWen -p "run<=276811&&isWenu&&V_pt>100&&Vtype==3&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wen_CS.txt -v V_pt --xlow 100 --xhigh 350 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 

#python ../../splitSamples.py -i $2 -c ttWmnLowPt -p "V_pt<150&&Jet_btagCSV[hJetInd1]>0.935&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==1" -s ../../systematics_Wmn_CS_EF_LowPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559 
#python ../../splitSamples.py -i $2 -c ttWmnHighPt -p "V_pt>150&&Jet_btagCSV[hJetInd1]>0.935&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==1" -s ../../systematics_Wmn_CS_EF_HighPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559 
#python ../../splitSamples.py -i $2 -c ttWenLowPt -p "V_pt<150&&Jet_btagCSV[hJetInd1]>0.935&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWenu&&V_pt>100&&Vtype==3&&controlSample==1" -s ../../systematics_Wen_CS_EF_LowPt.txt  -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c ttWenHighPt -p "V_pt>150&&Jet_btagCSV[hJetInd1]>0.935&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWenu&&V_pt>100&&Vtype==3&&controlSample==1" -s ../../systematics_Wen_CS_EF_HighPt.txt  -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c whfWmnLowPt -p "V_pt<150&&Jet_btagCSV[hJetInd1]>0.935&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&(H_mass<60||H_mass>120)&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wmn_CS_EF_LowPt.txt -v "Jet_btagCSV[hJetInd2]" --xlow 0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c whfWmnHighPt -p "V_pt>150&&Jet_btagCSV[hJetInd1]>0.935&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&(H_mass<60||H_mass>120)&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wmn_CS_EF_HighPt.txt -v "Jet_btagCSV[hJetInd2]" --xlow 0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c whfWenLowPt -p "V_pt<150&&Jet_btagCSV[hJetInd1]>0.935&&(sampleIndex!=0||(run>276811&&run<=278808))&&(H_mass<60||H_mass>120)&&isWenu&&V_pt>100&&Vtype==3&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wen_CS_EF_LowPt.txt -v "Jet_btagCSV[hJetInd2]" --xlow 0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c whfWenHighPt -p "V_pt>150&&Jet_btagCSV[hJetInd1]>0.935&&(sampleIndex!=0||(run>276811&&run<=278808))&&(H_mass<60||H_mass>120)&&isWenu&&V_pt>100&&Vtype==3&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wen_CS_EF_HighPt.txt -v "Jet_btagCSV[hJetInd2]" --xlow 0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c wlfWmnLowPt -p "V_pt<150&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wmn_CS_EF_LowPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c wlfWmnHighPt -p "V_pt>150&&selLeptons_relIso_0<0.06&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wmn_CS_EF_HighPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c wlfWenLowPt -p "V_pt<150&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWenu&&V_pt>100&&Vtype==3&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wen_CS_EF_LowPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c wlfWenHighPt -p "V_pt>150&&(sampleIndex!=0||(run>276811&&run<=278808))&&isWenu&&V_pt>100&&Vtype==3&&controlSample==3&&Jet_btagCSV[hJetInd1]<0.7" -s ../../systematics_Wen_CS_EF_HighPt.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 -w 0.559
#python ../../splitSamples.py -i $2 -c whfWmn -p "selLeptons_relIso_0<0.06&&run<=276811&&(H_mass<90||H_mass>150)&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wmn_CS.txt -v "(H_mass>120) + 2*(Jet_btagCSV[hJetInd2]>0.5)" --xlow 0 --xhigh 4 -t $4 -tt $8 -d 1 -n 4  -tol 0.75 -r 0
#python ../../splitSamples.py -i $2 -c whfWen -p "run<=276811&&(H_mass<90||H_mass>150)&&isWenu&&V_pt>100&&Vtype==3&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wen_CS.txt -v "(H_mass>120) + 2*(Jet_btagCSV[hJetInd2]>0.5)" --xlow 0 --xhigh 4 -t $4 -tt $8  -d 1 -n 4  -tol 0.75 -r 0


#python ../../splitSamples.py -i $2 -c whfWmn -p "selLeptons_relIso_0<0.06&&run<=276811&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wmn_CS_2015.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8 -d 1 -n 15  -tol 0.75
#python ../../splitSamples.py -i $2 -c whfWen -p "run<=276811&&isWenu&&V_pt>100&&Vtype==3&&controlSample==2&&nAddJets252p9_puid<1" -s ../../systematics_Wen_CS_2015.txt -v Jet_btagCSV[hJetInd2] --xlow 0.0 --xhigh 1.0 -t $4 -tt $8  -d 1 -n 15  -tol 0.75 

export SFs_Wmn="TT_Wmn,Wj0b_Wmn,Wj1b_Wmn,Wj2b_Wmn"
#export SFs_Wmn="TT_Wmn,Wj0b_Wmn,Wj1b_Wj2b_Wmn"
export SFs_Wen="TT_Wen,Wj0b_Wen,Wj1b_Wen,Wj2b_Wen"
#export SFs_Wen="TT_Wen,Wj0b_Wen,Wj1b_Wj2b_Wen"
echo "arg 3 = ${5}"
if [ $5 -eq '1' ]
then
    export SFs_Wmn="TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln"
    #export SFs_Wmn="TT_Wln,Wj0b_Wln,Wj1b_Wj2b_Wln"
    export SFs_Wen=$SFs_Wmn
fi

### FIXME this is just for the SR optimization studies
#export SFs_Wmn=""
#export SFs_Wen=""

# Signal Region Datacards
#python ../../printYields.py -c WmnHighPt -i hists_WmnHighPt.root -s ../../systematics_Wmn.txt -b binStats_WmnHighPt.txt -r $SFs_Wmn
#python ../../printYields.py -c WmnLowPt -i hists_WmnLowPt.root -s ../../systematics_Wmn.txt -b binStats_WmnLowPt.txt -r $SFs_Wmn 
#python ../../printYields.py -c WenHighPt -i hists_WenHighPt.root -s ../../systematics_Wen.txt -b binStats_WenHighPt.txt -r $SFs_Wen
#python ../../printYields.py -c WenLowPt -i hists_WenLowPt.root -s ../../systematics_Wen.txt -b  binStats_WenLowPt.txt -r $SFs_Wen

## FIXME: Put rate params back here
#python ../../printYields.py -c WmnHighPt -i hists_WmnHighPt.root -s ../../systematics_Wmn.txt -b binStats_WmnHighPt.txt
#python ../../printYields.py -c WmnLowPt -i hists_WmnLowPt.root -s ../../systematics_Wmn.txt -b binStats_WmnLowPt.txt
#python ../../printYields.py -c WenHighPt -i hists_WenHighPt.root -s ../../systematics_Wen.txt -b binStats_WenHighPt.txt
#python ../../printYields.py -c WenLowPt -i hists_WenLowPt.root -s ../../systematics_Wen.txt -b  binStats_WenLowPt.txt

# Control Region Datacards
#python ../../printYields.py -c ttWmn -i hists_ttWmn.root -s ../../systematics_Wmn_CS.txt -b binStats_ttWmn.txt -r $SFs_Wmn 
#python ../../printYields.py -c ttWen -i hists_ttWen.root -s ../../systematics_Wen_CS.txt -b binStats_ttWen.txt -r $SFs_Wen 
#python ../../printYields.py -c whfWmn -i hists_whfWmn.root -s ../../systematics_Wmn_CS.txt -b binStats_whfWmn.txt -r $SFs_Wmn 
#python ../../printYields.py -c whfWen -i hists_whfWen.root -s ../../systematics_Wen_CS.txt -b binStats_whfWen.txt -r $SFs_Wen 
#python ../../printYields.py -c wlfWmn -i hists_wlfWmn.root -s ../../systematics_Wmn_CS.txt -b binStats_wlfWmn.txt -r $SFs_Wmn 
#python ../../printYields.py -c wlfWen -i hists_wlfWen.root -s ../../systematics_Wen_CS.txt -b binStats_wlfWen.txt -r $SFs_Wen 

## Should make a separate script for actually running combine since it needs root5, but I'll leave some useful commands here
#combineCards.py WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_Wen_13TeV.txt
#combineCards.py WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt > vhbb_Wmn_13TeV.txt
#combineCards.py WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt > vhbb_Wln_13TeV.txt
#combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt > vhbb_CS_13TeV.txt
#combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt  WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_Wln_13TeV.txt
#combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_Wln_13TeV.txt
#combineCards.py WmnLowPt=vhbb_WmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt WenLowPt=vhbb_WenLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt > vhbb_SR_13TeV.txt
#combine -M MaxLikelihoodFit vhbb_Wln_13TeV.txt --saveShapes --saveWithUncertainties -v 3 --expectSignal=0 --rMin 0 --rMax 0
#combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 vhbb_Wln_13TeV.txt
#combineCards.py ttWenHighPt=vhbb_ttWenHighPt_13TeV.txt ttWenLowPt=vhbb_ttWenLowPt_13TeV.txt ttWmnHighPt=vhbb_ttWmnHighPt_13TeV.txt ttWmnLowPt=vhbb_ttWenLowPt_13TeV.txt whfWenHighPt=vhbb_whfWenHighPt_13TeV.txt whfWenLowPt=vhbb_whfWenLowPt_13TeV.txt whfWmnHighPt=vhbb_whfWmnHighPt_13TeV.txt whfWmnLowPt=vhbb_whfWmnLowPt_13TeV.txt wlfWenHighPt=vhbb_wlfWenHighPt_13TeV.txt wlfWenLowPt=vhbb_wlfWenLowPt_13TeV.txt wlfWmnHighPt=vhbb_wlfWmnHighPt_13TeV.txt wlfWmnLowPt=vhbb_wlfWmnLowPt_13TeV.txt > vhbb_CS_13TeV.txt
#combineCards.py ttWenHighPt=vhbb_ttWenHighPt_13TeV.txt ttWenLowPt=vhbb_ttWenLowPt_13TeV.txt ttWmnHighPt=vhbb_ttWmnHighPt_13TeV.txt ttWmnLowPt=vhbb_ttWenLowPt_13TeV.txt whfWenHighPt=vhbb_whfWenHighPt_13TeV.txt whfWenLowPt=vhbb_whfWenLowPt_13TeV.txt whfWmnHighPt=vhbb_whfWmnHighPt_13TeV.txt whfWmnLowPt=vhbb_whfWmnLowPt_13TeV.txt wlfWenHighPt=vhbb_wlfWenHighPt_13TeV.txt wlfWenLowPt=vhbb_wlfWenLowPt_13TeV.txt wlfWmnHighPt=vhbb_wlfWmnHighPt_13TeV.txt wlfWmnLowPt=vhbb_wlfWmnLowPt_13TeV.txt WmnHighPt=vhbb_WmnHighPt_13TeV.txt WmnLowPt=vhbb_WmnLowPt_13TeV.txt WenHighPt=vhbb_WenHighPt_13TeV.txt WenLowPt=vhbb_WenLowPt_13TeV.txt > vhbb_Wln_13TeV.txt
