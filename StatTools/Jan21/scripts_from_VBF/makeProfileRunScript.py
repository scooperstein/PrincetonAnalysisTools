##
## Helper script to generate many profile plots.
##
## Author: Stephane Cooperstein
##

vs = {}
vs["H_mass"] = (100,0.,300.,"M_bb")
vs["H_pt"] = (80,100.,180.,"p_{T}(jj) [GeV]")
vs["V_pt"] = (80,100.,180.,"p_{T}(W) [GeV]")
vs["nAddJet_f"] = (8,2.,10.,"N Additional Jets")
vs["lepMetDPhi"] = (20,2.0,3.15,"#Delta#phi(lep,MET)")
vs["HVdPhi"] = (20,2.5,3.15,"#Delta#phi(jj,W)")
vs["Jet_btagCSV[hJetInd2]"] = (25,0,1,"Sub-leading Jet CSV")
vs["Top1_mass_fromLepton_regPT_w4MET"] = (100,100,400,"Top Mass [GeV]")
vs["softActivityVH_njets5"] = (10,0,10,"N Soft Activity Jets with p_{T} > 5 GeV")
vs["V_mt"] = (100,0,200,"m_{T}(W)[GeV]")
vs["met_pt"] = (100,0,300,"Missing E_{T} [GeV]")
vs["CMS_vhbb_BDT_Wln_13TeV"] = (20,-1,1,"BDT Score")

otextfile = open("profile_runscript.sh","w")
otextfile.write("g++ mbb_profile.C -g -o mbb_profile `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA -lRooFit -lRooFitCore -lRooStats -lMinuit -lHtml\n")
vars_done = []
for var1 in vs.keys():
    for var2 in vs.keys():
        #if var1 != var2 and var2 not in vars_done:
        if var1 != var2:
            otextfile.write("./mbb_profile \"%s\" \"%s\" %i %f %f %i %f %f \"%s\" \"%s\"\n" % (var1,var2,vs[var1][0],vs[var1][1],vs[var1][2],vs[var2][0],vs[var2][1],vs[var2][2],vs[var1][3],vs[var2][3])) 
    vars_done.append(var1)
