#! /usr/bin/env python
#
#   Plotting from flat tree
#   Based on a program written by Matteo Sani (UCSD)
#   
#   Augmented slightly by Chris Palmer (Princeton)
#


import ROOT
import getopt, sys, array

def inputfileTree(mysamples):

    inputfiletree = ROOT.TTree("inputfiles", "inputfiles provenance information");

    nfiles = array.array('i', [1]) #junk
    nindfiles = array.array('i', [len(mysamples)])  #junk
    intlumi = array.array('f', [12.])  #junk
    sampleIndex = array.array('i', nindfiles[0]*[0])
    infoind = array.array('i', nfiles[0]*[1])  #junk
    histoindfromfiles = array.array('i', nfiles[0]*[1])  #junk
    inshortnames = ROOT.TClonesArray("TObjString", len(mysamples))
    infilenames = ROOT.TClonesArray("TObjString", len(mysamples)) #junk
        
    inputfiletree.Branch("nfiles", nfiles, "nfiles/I");
    inputfiletree.Branch("nindfiles", nindfiles, "nindfiles/I");
    inputfiletree.Branch("intlumi",  intlumi, "intlumi/F");
    inputfiletree.Branch("sampleIndex", sampleIndex, "sampleIndex[nindfiles]/I");
    inputfiletree.Branch("histoind", histoindfromfiles, "histoindfromfiles[nfiles]/I");
    inputfiletree.Branch("infoind", infoind, "infoind[nindfiles]/I");
    inputfiletree.Branch("inshortnames", "TClonesArray", ROOT.AddressOf(inshortnames), 32000, 0)
    inputfiletree.Branch("infilenames", "TClonesArray", ROOT.AddressOf(infilenames), 32000, 0)

    for i,s in enumerate(mysamples):
        sampleIndex[i] = s[0]
        temp = ROOT.TObjString()
        temp.SetString(s[1])
        inshortnames[i] = temp
    
    inputfiletree.Fill()

    return inputfiletree

def plotvariableTree(allHistos):

    plotvartree = ROOT.TTree("plotvariables","plotvariables provenance information");

    Nvar = array.array('i', [len(allHistos.name)])
    histoncat = array.array('i', Nvar[0]*[0])
    typplotall = array.array('i', [0])  #junk
    doplot = array.array('i', Nvar[0]*[0]) #junk
    h2d = array.array('i', Nvar[0]*[0]) #junk
    typplot = array.array('i', Nvar[0]*[0]) #junk
    histoncatindtonames = array.array('i', Nvar[0]*[0]) #junk
    nbinsx = array.array('i', Nvar[0]*[0])
    nbinsy = array.array('i', Nvar[0]*[0])
    lowlim = array.array('f', Nvar[0]*[0])
    highlim = array.array('f', Nvar[0]*[0])
    lowlim2 = array.array('f', Nvar[0]*[0])
    highlim2 = array.array('f', Nvar[0]*[0])
    xaxislabels = ROOT.TClonesArray("TObjString", Nvar[0])
    yaxislabels = ROOT.TClonesArray("TObjString", Nvar[0])
    plotvarnames = ROOT.TClonesArray("TObjString", Nvar[0])
        
    plotvartree.Branch("Nvar", Nvar, "Nvar/I");
    plotvartree.Branch("typplotall", typplotall, "typplotall/I");
    plotvartree.Branch("doplot", doplot, "doplot[Nvar]/I");
    plotvartree.Branch("h2d", h2d, "h2d[Nvar]/I");
    plotvartree.Branch("typplot", typplot, "typplot[Nvar]/I");
    plotvartree.Branch("histoncat", histoncat, "histoncat[Nvar]/I");
    plotvartree.Branch("histoncatindtonames", histoncatindtonames, "histoncatindtonames[Nvar]/I");
    plotvartree.Branch("nbinsx", nbinsx, "nbinsx[Nvar]/I");
    plotvartree.Branch("nbinsy", nbinsy, "nbinsy[Nvar]/I");
    plotvartree.Branch("lowlim", lowlim, "lowlim[Nvar]/F");
    plotvartree.Branch("highlim", highlim, "highlim[Nvar]/F");
    plotvartree.Branch("lowlim2", lowlim2, "lowlim2[Nvar]/F");
    plotvartree.Branch("highlim2", highlim2, "highlim2[Nvar]/F");
    plotvartree.Branch("xaxislabels", "TClonesArray", ROOT.AddressOf(xaxislabels), 32000, 0)
    plotvartree.Branch("yaxislabels", "TClonesArray", ROOT.AddressOf(yaxislabels), 32000, 0)
    plotvartree.Branch("plotvarnames", "TClonesArray", ROOT.AddressOf(plotvarnames), 32000, 0)
 
    for i in xrange(Nvar[0]):
        h2d[i] = 0
        doplot[i] = 1
        typplot[i] = 0
        histoncat[i] = allHistos.ncat[i]
        histoncatindtonames[i] = -1
        nbinsx[i]   = allHistos.xbins[i]
        nbinsy[i]   = allHistos.ybins[i]
        lowlim[i]   = allHistos.xmin[i]
        highlim[i]  = allHistos.xmax[i]
        lowlim2[i]  = allHistos.ymin[i]
        highlim2[i] = allHistos.ymax[i]

        temp = ROOT.TObjString()
        temp.SetString(allHistos.xaxis[i])
        xaxislabels[i] = temp

        temp.SetString(allHistos.yaxis[i])
        yaxislabels[i] = temp

        temp.SetString(allHistos.name[i])
        plotvarnames[i] = temp

    plotvartree.Fill()

    return plotvartree
  
class histoContainer:
    def __init__(self):
        self.histo_type = []
        self.ncat = []
        self.xbins = []
        self.ybins = []
        self.xmin, self.xmax = [], []
        self.ymin, self.ymax = [], []
        self.name = []
        self.xaxis = []
        self.yaxis = []
        self.vars = []
        self.catTypes =[]

def parseSelection(filename):    
    file = open(filename)
    lines = file.readlines()
    file.close()

    return lines[0].split("\n")[0]

def parseCategories(filename):
    categories = dict()
    file = open(filename)
    lines = file.readlines()
    file.close()
    
    for l in lines:
        if ("#" not in l):
            catname = ""
            catdef = []
            items = l.split()
            for item in items:
                if (item.find("catName=") is not -1):
                    catname = str(item.split("=")[1])
                elif (item.find("catDef:") is not -1):
                    temp = (str(item.split(":")[1])).split(",")
                    for t in temp:
                        catdef.append(t)
                #else:
                #    print "Cannot parse categories"
            categories[catname] = catdef 

    return categories

def parseInputfiles(filename):
    samples = []
    file = open(filename)
    lines = file.readlines()
    file.close()
    for l in lines:
        if ("#" not in l and "type" in l):
            items = l.split()
            sampleIndex = -1
            name = ""
            for item in items:
                if (item.find("type=") is not -1):
                    sampleIndex = int(item.split("=")[1])
                elif (item.find("name=") is not -1):
                    name = str(item.split("=")[1])
                #else:
                #    print "Item",item,"not parsed."
            samples.append((sampleIndex, name))

    for i in range(len(samples)-1,-1,-1):
        if samples.count(samples[i]) > 1:
            samples.pop(i)

    return samples

def parsePlotvariables(filename, samples):
    h = histoContainer()

    file = open(filename)
    lines = file.readlines()
    file.close()

    for l in lines:
      if ("#"!=l[0] and "htyp" in l):
        if(deepDebug): print "# at the front",int("#"==l[0])
        items = l.split()
        for item in items:
            if ("default" in item or item == ""):
                continue
            elif (item.find("htyp=") is not -1):
                h.histo_type.append(int(item.split("=")[1]))
            elif item.find("ncat=") is not -1:
                h.ncat.append(int(item.split("=")[1]))
            elif item.find("xbins=") is not -1:
                h.xbins.append(int(item.split("=")[1]))
            elif item.find("ybins=") is not -1:
                h.ybins.append(int(item.split("=")[1]))
            elif item.find("xmin=") is not -1:
                h.xmin.append(float(item.split("=")[1]))
            elif item.find("xmax=") is not -1:
                h.xmax.append(float(item.split("=")[1]))
            elif item.find("ymin=") is not -1:
                h.ymin.append(float(item.split("=")[1]))
            elif item.find("ymax=") is not -1:
                h.ymax.append(float(item.split("=")[1]))
            elif item.find("name=") is not -1:
                h.name.append(str(item.split("=")[1]))
            elif item.find("xaxis=") is not -1:
                dummy = str(item.split("=")[1])
                dummy = dummy.replace("@"," ")
                h.xaxis.append(dummy)
            elif item.find("yaxis=") is not -1:
                dummy = str(item.split("=")[1])
                dummy = dummy.replace("@"," ")
                h.yaxis.append(dummy)
            elif item.find("name=") is not -1:
                h.name.append(str(item.split("=")[1]))
            elif item.find("var=") is not -1:
                h.vars.append(str(item.split("=")[1]))
            elif item.find("catType=") is not -1:
                h.catTypes.append(str(item.split("=")[1]))
            #else:
            #    print "Item",item,"not parsed."
    return h

def usage():
    use="""
plotFromOptree [options]

-h, --help               print this message
-w, --weight             weight branchname (default "weight")
-i, --rootinputfile      inputfile name
-o, --rootoutputfile     outputfile name (default "output.root")
-t, --treename           inputtree name (default "tree")
-S, --selectionfile      selection .dat (default "selection.dat")
-I, --inputfile          inputfiles .dat (default "inputfiles.dat")
-P, --plotvariables      plotvariables .dat (default "plotvariables.dat")
-C, --categories         categories .dat (default "categories.dat")
-d, --debug              display debugging messages
-D, --deepDebug          display more debugging messages
"""
    print use

if __name__ == "__main__":  

    # default parameters
    weightName = "weight"
    rootInputFile = ""
    rootOutputFile = "output.root"
    treeName = "tree"
    selectionFile = "selection.dat"
    inputfilesFile = "inputfiles.dat"
    plotvariablesFile = "plotvariables.dat"
    categoriesFile = "categories.dat"
    debug=False
    deepDebug=False
    inputDirectory = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hw:i:o:t:S:I:P:C:f:dD", ["help", "weight", "rootinputfile", "roototputfile", "treename", "selectionfile", "inputfile", "plotvariables", "categories","debug","deepDebug","inputDirectory"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        print opt,arg
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-w", "--weight"):
            weightName = arg
        elif opt in ("-i", "--rootinputfile"):
            rootInputFile = arg
        elif opt in ("-o", "--rootoutputfile"):
            rootOutputFile = arg
        elif opt in ("-t", "--treename"):
            treeName = arg
        elif opt in ("-S", "--selectionFile"):
            selectionFile = arg
        elif opt in ("-I", "--inputfile"):
            inputfilesFile = arg
        elif opt in ("-P", "--plotvariables"):
            plotvariablesFile = arg
        elif opt in ("-C", "--categories"):
            categoriesFile = arg
        elif opt in ("-d", "--debug"):
            debug=True
        elif opt in ("-D", "--deepDebug"):
            debug=True
            deepDebug=True
        elif opt in ("-f", "--inputDirectory"):
            inputDirectory = arg
        else:
            print "Unrecongnized option " + opt

    if (rootInputFile == "" and inputDirectory == ""):
        print "You have to select one input file or input directory."
        usage()
        sys.exit(3)

    if (weightName == ""):
        weightName = "1"
    
    print "Parsing all files"
    print "inputDirectory = ",inputDirectory
    finalHistos = []
    if(debug): print "Parsing selection file..."
    selection_cut = parseSelection(selectionFile)
    if(debug): print "Parsing inputfiles..."
    samples = parseInputfiles(inputfilesFile)
    if(debug): print "Parsing plotvariables..."
    allHistos = parsePlotvariables(plotvariablesFile, samples)
    if(debug): print "Writing inputfiles Tree..."
    inputfiletree = inputfileTree(samples)
    if(debug): print "Writing plotvariables Tree..."
    plotvariabletree = plotvariableTree(allHistos)
    if(debug): print "Parsing categories file..."
    categories = parseCategories(categoriesFile)

    if (inputDirectory == ""):
        file = ROOT.TFile.Open(rootInputFile)
        tree = file.Get(treeName)


    sampleNameMap = {}
    sampleNameMap["WH125p"] = "WplusH125_powheg" 
    sampleNameMap["WH125m"] = "WminusH125_powheg" 
    sampleNameMap["ZH125"]  = "ZH125_powheg"
    sampleNameMap["WZ_udcsg"] = "WZ_fil"
    sampleNameMap["WZ_b"] = "WZ_fil"
    sampleNameMap["WZ_bb"] = "WZ_fil"
    sampleNameMap["ZZ_udcsg"] = "ZZ_fil"
    sampleNameMap["ZZ_b"] = "ZZ_fil"
    sampleNameMap["ZZ_bb"] = "ZZ_fil"
    sampleNameMap["WW_udcsg"] = "WW_fil"
    sampleNameMap["WW_b"] = "WW_fil"
    sampleNameMap["WW_bb"] = "WW_fil"
    sampleNameMap["TT_powheg"] = "TT_powheg"
    sampleNameMap["TToLeptons_s"] = "TToLeptons_s"
    sampleNameMap["TToLeptons_t"] = ["TToLeptons_t_powheg","TBarToLeptons_t_powheg"]    
    #sampleNameMap["TBarToLeptons_t_powheg"] = "TBarToLeptons_t_powheg" 
    #sampleNameMap["TToLeptons_t_powheg"] = "TToLeptons_t_powheg" 
    sampleNameMap["T_tW"] = "T_tW"
    sampleNameMap["Tbar_tW"] = "Tbar_tW"
    sampleNameMap["W_udcsg"] = "WJets_madgraph"
    sampleNameMap["W_b"] = "WJets_madgraph"
    sampleNameMap["W_bb"] = "WJets_madgraph"
    sampleNameMap["W_udcsg_HT100To200"] = "WJets-HT100To200"
    sampleNameMap["W_b_HT100To200"] = "WJets-HT100To200"
    sampleNameMap["W_bb_HT100To200"] = "WJets-HT100To200"
    sampleNameMap["W_udcsg_HT200To400"] = "WJets-HT200To400"
    sampleNameMap["W_b_HT200To400"] = "WJets-HT200To400"
    sampleNameMap["W_bb_HT200To400"] = "WJets-HT200To400"
    sampleNameMap["W_udcsg_HT400To600"] = "WJets-HT400To600"
    sampleNameMap["W_b_HT400To600"] = "WJets-HT400To600"
    sampleNameMap["W_bb_HT400To600"] = "WJets-HT400To600"
    sampleNameMap["W_udcsg_HT600To800"] = "WJets-HT600To800"
    sampleNameMap["W_b_HT600To800"] = "WJets-HT600To800"
    sampleNameMap["W_bb_HT600To800"] = "WJets-HT600To800"
    sampleNameMap["W_udcsg_HT800To1200"] = "WJets-HT800To1200"
    sampleNameMap["W_b_HT800To1200"] = "WJets-HT800To1200"
    sampleNameMap["W_bb_HT800To1200"] = "WJets-HT800To1200"
    sampleNameMap["W_udcsg_HT1200To2500"] = "WJets-HT1200To2500"
    sampleNameMap["W_b_HT1200To2500"] = "WJets-HT1200To2500"
    sampleNameMap["W_bb_HT1200To2500"] = "WJets-HT1200To2500"
    sampleNameMap["W_udcsg_HT2500ToInf"] = "WJets-HT2500ToInf"
    sampleNameMap["W_b_HT2500ToInf"] = "WJets-HT2500ToInf"
    sampleNameMap["W_bb_HT2500ToInf"] = "WJets-HT2500ToInf"
    sampleNameMap["W_udcsg_WBJets_Pt100To200"] = "WBJets-Pt100To200"
    sampleNameMap["W_b_WBJets_Pt100To200"] = "WBJets-Pt100To200"
    sampleNameMap["W_bb_WBJets_Pt100To200"] = "WBJets-Pt100To200"
    sampleNameMap["W_udcsg_WBJets_Pt200ToInf"] = "WBJets-Pt200ToInf"
    sampleNameMap["W_b_WBJets_Pt200ToInf"] = "WBJets-Pt200ToInf"
    sampleNameMap["W_bb_WBJets_Pt200ToInf"] = "WBJets-Pt200ToInf"
    sampleNameMap["W_udcsg_BGenFilter_Pt100To200"] = "WJets_BGenFilter-Pt100To200"
    sampleNameMap["W_b_BGenFilter_Pt100To200"] = "WJets_BGenFilter-Pt100To200"
    sampleNameMap["W_bb_BGenFilter_Pt100To200"] = "WJets_BGenFilter-Pt100To200"
    sampleNameMap["W_udcsg_BGenFilter_Pt200ToInf"] = "WJets_BGenFilter-Pt200ToInf"
    sampleNameMap["W_b_BGenFilter_Pt200ToInf"] = "WJets_BGenFilter-Pt200ToInf"
    sampleNameMap["W_bb_BGenFilter_Pt200ToInf"] = "WJets_BGenFilter-Pt200ToInf"
    sampleNameMap["Z_udcsg"] = "DYToLL_madgraph"
    sampleNameMap["Z_b"] = "DYToLL_madgraph"
    sampleNameMap["Z_bb"] = "DYToLL_madgraph"
    sampleNameMap["Z_udcsg_HT100To200"] = "DYToLL_HT100to200"
    sampleNameMap["Z_b_HT100To200"] = "DYToLL_HT100to200"
    sampleNameMap["Z_bb_HT100To200"] = "DYToLL_HT100to200"
    sampleNameMap["Z_udcsg_HT200To400"] = "DYToLL_HT200to400"
    sampleNameMap["Z_b_HT200To400"] = "DYToLL_HT200to400"
    sampleNameMap["Z_bb_HT200To400"] = "DYToLL_HT200to400"
    sampleNameMap["Z_udcsg_HT400To600"] = "DYToLL_HT400to600"
    sampleNameMap["Z_b_HT400To600"] = "DYToLL_HT400to600"
    sampleNameMap["Z_bb_HT400To600"] = "DYToLL_HT400to600"
    ####sampleNameMap["QCD_HT100To200"] = "QCD_HT100to200"
    sampleNameMap["QCD_HT100To200"] = "IntEWKWJets"
    sampleNameMap["QCD_HT200To300"] = "QCD_HT200to300"
    sampleNameMap["QCD_HT300To500"] = "QCD_HT300to500"
    sampleNameMap["QCD_HT500To700"] = "QCD_HT500to700"
    sampleNameMap["QCD_HT700To1000"] = "QCD_HT700to1000"
    sampleNameMap["QCD_HT1000To1500"] = "QCD_HT1000to1500"
    sampleNameMap["QCD_HT1500To2000"] = "QCD_HT1500to2000"
    sampleNameMap["QCD_HT2000ToInf"] = "QCD_HT2000toInf"
    sampleNameMap["WJets_0J"] = "WJets_0J"
    sampleNameMap["WJets_1J"] = "WJets_1J"
    sampleNameMap["WJets_2J"] = "WJets_2J"
    sampleNameMap["EWKWJets"] = "EWKWJets"
    sampleNameMap["EWKWJets_herwig"] = "EWKWJets_herwig"
    sampleNameMap["IntEWKWJets"] = "IntEWKWJets"
    sampleNameMap["WZ_fil"] = "WZ_fil"
    sampleNameMap["WW_fil"] = "WW_fil"
    sampleNameMap["ZZ_fil"] = "ZZ_fil"
    sampleNameMap["DYToLL_madgraph"] = "DYToLL_madgraph"
    sampleNameMap["DYToLL_HT100to200"] = "DYToLL_HT100to200"
    sampleNameMap["DYToLL_HT200to400"] = "DYToLL_HT200to400"
    sampleNameMap["DYToLL_HT400to600"] = "DYToLL_HT400to600"
    sampleNameMap["DYToLL_HT600to800"] = "DYToLL_HT600to800"
    sampleNameMap["DYToLL_HT800to1200"] = "DYToLL_HT800to1200"
    sampleNameMap["DYToLL_HT1200to2500"] = "DYToLL_HT1200to2500"
    sampleNameMap["DYToLL_HT2500toInf"] = "DYToLL_HT2500toInf"
    sampleNameMap["ZJets_0J"] = "ZJets_0J"
    sampleNameMap["ZJets_1J"] = "ZJets_1J"
    sampleNameMap["ZJets_2J"] = "ZJets_2J"
    sampleNameMap["WJets_inc"] = "WJets_madgraph"
    sampleNameMap["WJets_HT100To200"] = "WJets-HT100To200"
    sampleNameMap["WJets_HT200To400"] = "WJets-HT200To400"
    sampleNameMap["WJets_HT400To600"] = "WJets-HT400To600"
    sampleNameMap["WJets_HT600To800"] = "WJets-HT600To800"
    sampleNameMap["WJets_HT800To1200"] = "WJets-HT800To1200"
    sampleNameMap["WJets_HT1200To2500"] = "WJets-HT1200To2500"
    sampleNameMap["WJets_HT2500ToInf"] = "WJets-HT2500ToInf"
    ##sampleNameMap["QCD_HT100To200"] = "ZJets_0J"

    print "Looping over samples and histograms"
    print "inputDirectory = ",inputDirectory
    output = ROOT.TFile.Open(rootOutputFile, "recreate")
    for ns, s in enumerate(samples):
        #if(debug): print "Processing sample: ", s[1]
        print "Processing sample: ", s[1]
        if s[1] not in sampleNameMap:
            sampleNameMap[s[1]] = s[1]
        if (inputDirectory != ""):
            if isinstance(sampleNameMap[s[1]], (list, tuple)):
                tree = ROOT.TChain(treeName)
                for item in sampleNameMap[s[1]]:
                    print "Adding to tree: " + inputDirectory + "/sum_" + item + "_weighted.root"
                    #tree.Add(inputDirectory + "/sum_" + item + "_weighted3.root")
                    #tree.Add(inputDirectory + "/sum_" + item + "_weighted2.root")
                    #tree.Add(inputDirectory + "/sum_" + item + "_3.root")
                    tree.Add(inputDirectory + "/sum_" + item + ".root")
            else:
                print inputDirectory + "/sum_" + sampleNameMap[s[1]] + "_weighted2.root"
                #file = ROOT.TFile(inputDirectory + "/sum_" + sampleNameMap[s[1]] + "_weighted3.root")
                #file = ROOT.TFile(inputDirectory + "/sum_" + sampleNameMap[s[1]] + "_weighted2.root")
                #file = ROOT.TFile(inputDirectory + "/sum_" + sampleNameMap[s[1]] + "_3.root")
                file = ROOT.TFile.Open(inputDirectory + "/sum_" + sampleNameMap[s[1]] + ".root")
                tree = file.Get(treeName)
        for nv, v in enumerate(allHistos.vars):
            if(deepDebug): print v
            temp_hist = ""
            if (allHistos.histo_type[nv] == 0):
                temp_hist = "temp_hist("+ str(allHistos.xbins[nv])+","+str(allHistos.xmin[nv])+","+str(allHistos.xmax[nv])+")"

            for nc in xrange(allHistos.ncat[nv]):
                if(deepDebug): print "cat ",nc
                cut = "sampleIndex==" + str(s[0])

                if (allHistos.ncat[nv] > 1 and allHistos.ncat[nv] > len(categories[allHistos.catTypes[nv]])):
                    print categories
                    print "Wrong category definition !!!"
                    sys.exit(-1)

                category_cut = ""
                if (allHistos.ncat[nv] > 1):
                    category_cut = categories[allHistos.catTypes[nv]][nc]

                if (category_cut != ""):
                    cut = cut + " && " + category_cut
                    
                if (selection_cut != ""):
                    cut = cut + " && " + selection_cut

                cut = "(" + cut + ")* " + weightName
                    
                try:
                    if(deepDebug): tree.Print("*"+v+"*")
                    tree.Draw(v+" >> "+temp_hist, cut, "goff")
            
                    final_h = ROOT.gDirectory.Get("temp_hist")
                    final_name = allHistos.name[nv] + "_cat"+str(nc)+"_"+s[1]
                    final_h.SetName(final_name)
                    final_h.SetTitle(final_name)
                    final_h.GetXaxis().SetTitle(allHistos.xaxis[nv])
                    final_h.GetYaxis().SetTitle(allHistos.yaxis[nv])
                    finalHistos.append(final_h)
                    output.cd()
                    final_h.Write()
                except:
                    print "Failed to draw"
                    print "v",v
                    print "cut",cut
                    print "temp_hist",temp_hist
                    sys.exit(4)
        if (inputDirectory != ""):
            tree.Reset()
            file.Close()            
    #output = ROOT.TFile(rootOutputFile, "recreate")
    #for h in finalHistos:
    #    h.Write()
    
    output.cd()    
    inputfiletree.Write()
    plotvariabletree.Write()

    output.Close()
