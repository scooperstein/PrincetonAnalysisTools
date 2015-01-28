#~ /usr/bin/python
import sys
import ROOT
ROOT.gSystem.Load("SampleContainer_cc.so")
ROOT.gSystem.Load("AnalysisManager_cc.so")
ROOT.gSystem.Load("VHbbAnalysis_cc.so")
ROOT.gSystem.Load("BDTInfo_h.so")
ROOT.gSystem.Load("VHbbTrigger_h.so")

debug=0

def ReadTextFile(filename, filetype):
    if debug > 100:
         print "filetype is ", filetype
         print "filename is ", filename

    textfile=open(filename, 'r')

    filelines=textfile.readlines()
    filelines=CleanLines(filelines)

    if debug > 1000:
        print filelines
    if filetype is "cfg":

        settings=MakeConfigMap(filelines)

        if settings.has_key("analysis"):
            am=ROOT.__getattr__(settings["analysis"])()
            
        if settings.has_key("samples"):
            aminitialized=0
            samples=ReadTextFile(settings["samples"], "samplefile")
            for name in samples:
                sample=samples[name]
                #print sample, sampledic[sample]
                samplecon = ROOT.SampleContainer()
                if sample.has_key("type"):
                    samplecon.sampleNum         = sample["type"]
                if sample.has_key("xsec"):
                    samplecon.xsec              = sample["xsec"]
                if sample.has_key("kfac"):
                    samplecon.kFactor           = sample["kfac"]
                if sample.has_key("scale"):
                    samplecon.scale             = sample["scale"]
                if sample.has_key("npro"):
                    samplecon.processedEvents   = sample["npro"]
                else:
                    samplecon.nProFromFile      = True
                for filename in sample["files"]:
                    # AnalysisManager needs to be initialized
                    # with one file at the beginning
                    if aminitialized == 0:
                        am.Initialize(filename)
                        aminitialized=1
                        # FIXME can this go elsewhere?
                        if settings.has_key("outputname"):
                            am.outputTreeName=settings["outputname"]
                    samplecon.AddFile(filename)
                am.AddSample(samplecon)

        else:
            print "There are no samples in the config file."
            sys.exit(0)
 
        if settings.has_key("earlybranches"):
            branches=ReadTextFile(settings["earlybranches"], "branchlist")
            for branch in branches:
                am.SetupBranch(branch,branches[branch][0], branches[branch][1], "early")
        else:
            print "There are no existing branches in the config file."

        if settings.has_key("existingbranches"):
            branches=ReadTextFile(settings["existingbranches"], "branchlist")
            for branch in branches:
                am.SetupBranch(branch,branches[branch][0], branches[branch][1])
        else:
            print "There are no existing branches in the config file."

        am.ConfigureOutputTree()
            
        if settings.has_key("newbranches"):
            branches=ReadTextFile(settings["newbranches"], "branchlist")
            for branch in branches:
                am.SetupNewBranch(branch,branches[branch][0], branches[branch][1])
        else:
            print "There are no new branches in the config file."
        
        if settings.has_key("settings"):
            branches=ReadTextFile(settings["settings"], "branchlist")
            for branch in branches:
                am.SetupNewBranch(branch,branches[branch][0], branches[branch][1], True, "settings", branches[branch][2])
        else:
            print "There are no settings branches in the config file."
             
        if settings.has_key("bdtsettings"):
            print "Adding a BDT configuration..."
            bdtInfo=ReadTextFile(settings["bdtsettings"], "bdt")
            print "read the BDT settings text file for BDT %s" % bdtInfo.bdtname
            # now set up any of the branches if they don't exist yet (must be floats for BDT)
            for varname in bdtInfo.localVarNames:
                am.SetupNewBranch(varname, 2)
            for varname in bdtInfo.localSpectatorVarNames:
                am.SetupNewBranch(varname, 2)       
            am.AddBDT(bdtInfo)
            print "added BDT to analysis manager"
        if settings.has_key("reg1settings"):
            print "Adding a Jet 1 Energy Regresion..."
            reg1 = ReadTextFile(settings["reg1settings"], "bdt") 
            for varname in reg1.localVarNames:
                am.SetupNewBranch(varname, 2)
            for varname in reg1.localSpectatorVarNames:
                am.SetupNewBranch(varname, 2) 
            am.SetJet1EnergyRegression(reg1) 
        if settings.has_key("reg2settings"):
            print "Adding a Jet 2 Energy Regresion..."
            reg2 = ReadTextFile(settings["reg2settings"], "bdt") 
            for varname in reg2.localVarNames:
                am.SetupNewBranch(varname, 2)
            for varname in reg2.localSpectatorVarNames:
                am.SetupNewBranch(varname, 2) 
            am.SetJet2EnergyRegression(reg2) 
        return am    
    elif filetype is "samplefile":
        samples=MakeSampleMap(filelines)
        return samples
    elif filetype is "branchlist":
        branches=MakeBranchMap(filelines)
        return branches
    elif filetype is "bdt":
        bdtInfo=SetupBDT(filelines)
        return bdtInfo
    else:
        print "Unknown filetype ", filetype


def CleanLines(lines):
    cleanedLines=[]
    for line in lines:
        # skip if it is 100% whitespace
        if line.isspace():
            continue
        # remove comments
        if line.find("#") is -1:        # no comment
            cleanedLines.append(line.split("\n")[0])
        elif line.find("#") is not 0:   # there is an inline comment
            cleanedLines.append(line[:line.find("#")])
        else:                           # whole line is commented
            if debug > 1000:
                print "this line is commented and will be removed"
                print line
    return cleanedLines
    

def MakeConfigMap(lines):
    if debug > 0: print "In MakeConfigMap"
    settings={}

    for line in lines:
        for item in line.split():
            if len(item.split("=")) == 2:
                name,value=item.split("=")
                settings[name]=value

            else:
                print "not in setting=value format\n\t",item

    return settings


def MakeSampleMap(lines):
    if debug > 100: print "Reading samples"
    samples={}

    for line in lines:
        #print line
        samplename=""
        samplepaths=[]
        sampletype=-1
        samplexsec=-1
        samplekfac=1
        samplescale=1
        sample={}

        for item in line.split():
            name,value=item.split("=")
            if name.find("name") is 0:
                sample["name"]=str(value)
            if name.find("file") is 0:
                samplepaths.append(str(value))
            if name.find("dir") is 0:
                from os import listdir
                from os.path import isfile, join
                onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
                for rootfile in onlyfiles:
                    #print rootfile
                    if rootfile.find(".root") != -1:
                        samplepaths.append(str(value)+"/"+str(rootfile))
            if name.find("type") is 0:
                sample["type"]=int(value)
            if name.find("xsec") is 0:
                sample["xsec"]=float(value)
            if name.find("kfac") is 0:
                sample["kfac"]=float(value)
            if name.find("scale") is 0:
                sample["scale"]=float(value)
            if name.find("npro") is 0:
                sample["npro"]=int(value)
       
        sample["files"]=samplepaths 
        if sample.has_key("name"):
            samples[sample["name"]]=sample
        else:
            print "sample name is empty",samplename,"not filling"
    
    return samples


def MakeBranchMap(lines):
    branches={}

    for line in lines:
        branchname=""
        branchtype=-1
        arraylength=-1
        val=-999
        
        for item in line.split():
            name,value = item.split("=")
            if name.find("name") is 0:
                branchname=value
            if name.find("type") is 0:
                branchtype=int(value)
            if name.find("max") is 0:
                arraylength=int(value)
            if name.find("val") is 0:
                val=float(value)

        branches[branchname]= [branchtype,arraylength,val]

    return branches            

def SetupBDT(lines):
    bdtname = ""
    xmlFile = ""
    inputNames = []
    localVarNames = []
    vars = {} 

    for line in lines:
        inputName = ""
        localVarName = ""
        type = -1
        order = -1
        for item in line.split():
            name,value = item.split("=")
            if name.find("bdtname") is 0:
                bdtname=value
                break
            #if name.find("method") is 0:
            #    method=value.replace('@', ' ')
            #    print "method set to %s" % method
            #    break
            if name.find("xmlFile") is 0:
                xmlFile=value
                break
            if name.find("name") is 0:
                inputName=value
            if name.find("lname") is 0:
                localVarName=value
            if name.find("type") is 0:
                type=int(value)
            if name.find("order") is 0:
                order=int(value)
        vars[order] = (inputName,localVarName,type)
   
    bdt = ROOT.BDTInfo(bdtname, xmlFile)
    
    keys = vars.keys()
    keys.sort()
 
    for key in keys:
        if (key == -1): continue
        name, lname, type = vars[key]
        if (type == 1):
            print "adding variable %s (%s) " % (name,lname)
            bdt.AddVariable(name, lname)
        if (type == 0):
            print "adding spectator variable %s (%s) " % (name,lname)
            bdt.AddSpectatorVariable(name, lname)
    return bdt 
