import ROOT
import sys
import math

filename=sys.argv[1]

tfile=ROOT.TFile(filename)
tree=tfile.Get("tree")

sampleTypes=[-125,21]
genCuts="(genElePt>20&&abs(genEleEta)<2.5&&b2Pt_gen>25&&b1Pt_gen>25&&abs(b2Eta_gen)<2.5&&abs(b1Eta_gen)<2.5&& genHPt>100)"

addCuts=["1==1",
#"L1_EG30isoer==1",
#"L1_EG35er==1",
#"L1_EG40==1",
#"L1_HTT175==1",
#"L1_ETM70==1",
"(L1_EG30isoer==1||L1_EG40==1)",
"(L1_EG30isoer==1||L1_EG40==1||L1_ETM70==1)",
"(L1_EG30isoer==1||L1_EG40==1||L1_HTT175==1)",
"(L1_EG30isoer==1||L1_EG40==1||L1_HTT175==1||L1_ETM70==1)",
"HLT_Ele32_WP75_NoER==1",
"HLT_Ele40_WP85==1",
"HLT_Ele40_WP85_NoER==1",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1)",
"(HLT_Ele20_WP85_DijetPT50==1)",
"(HLT_Ele20_WP85_DijetPT100==1)",
"(HLT_Ele20_WP85_DijetPT150==1)",
"(HLT_Ele30_WP85_DijetPT50==1)",
"(HLT_Ele30_WP85_DijetPT100==1)",
"(HLT_Ele30_WP85_DijetPT150==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele20_WP85_DijetPT50==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele20_WP85_DijetPT100==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele20_WP85_DijetPT150==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele30_WP85_DijetPT50==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele30_WP85_DijetPT100==1)",
"(HLT_Ele32_WP75_NoER==1||HLT_Ele40_WP85==1||HLT_Ele30_WP85_DijetPT150==1)",
]
#"HLT_Ele32_WP75==1",
#"HLT_Ele35_WP85_NoER==1",
#"HLT_Ele35_WP85_OLDL1_NoER==1",
#"HLT_Ele35_WP85_NoER==1",
#"HLT_Ele35_WP85==1",
#"HLT_Ele40_WP85_OLDL1_NoER==1",
#"HLT_Ele40_WP85_OLDL1==1",
#"HLT_40CSV0p5_ETA2p4_BESTPT==1",
#"HLT_Ele20_WP85_40CSV0p5_ETA2p4_BESTPT==1",
#"(HLT_Ele20_WP85_40CSV0p5_ETA2p4_BESTPT==1||HLT_Ele32_WP75==1||HLT_Ele40_WP85)",
#"(HLT_Ele20_WP85_40CSV0p5_ETA2p4_BESTPT==1||HLT_Ele32_WP75==1||HLT_Ele40_WP85||HLT_Ele20_WP85_DijetPT50==1)",
#"(HLT_Ele20_WP85_40CSV0p5_ETA2p4_BESTPT==1||HLT_Ele32_WP75==1||HLT_Ele20_WP85_DijetPT50==1)"]

for addCut in addCuts:
    print addCut,
    for sampleType in sampleTypes:
        basecut="sampleIndex=="+str(sampleType)

        if sampleType<0:
            basecut=basecut+"&&"+genCuts

        cut=basecut+"&&"+addCut
        #entries=tree.GetEntries(cut)
        #temp=ROOT.TH1F("temp"+str(sampleType)+"_"+addCut,"",2,0,2)
        temp=ROOT.TH1F("temp","",2,0,2)
        tree.Draw("L1_EG30>>temp",cut,"goff")
        entries=temp.Integral()

        tree.Draw("L1_EG30>>temp","("+cut+")*weight","goff")
        weighted=temp.Integral()
        #print sampleType,entries,math.sqrt(entries),weighted,
        print sampleType,entries,weighted,weighted/math.sqrt(entries),
        del temp
    print
