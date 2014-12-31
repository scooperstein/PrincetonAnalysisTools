#~ /usr/bin/python
import sys
import ROOT
ROOT.gSystem.Load("SampleContainer_cc.so")
ROOT.gSystem.Load("AnalysisManager_cc.so")
ROOT.gSystem.Load("VHbbAnalysis_h.so")

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
            for sample in samples:
                #print sample, sampledic[sample]
                samplecon = ROOT.SampleContainer()
                samplecon.sampleNum = samples[sample][0]
                samplecon.xsec      = samples[sample][1]
                samplecon.kfactor   = samples[sample][2]
                samplecon.scale     = samples[sample][3]
                for filename in samples[sample][4]:
                    # AnalysisManager needs to be initialized
                    # with one file at the beginning
                    if aminitialized == 0:
                        am.Initialize(filename)
                        aminitialized=1
                    samplecon.AddFile(filename)
                am.AddSample(samplecon)

        else:
            print "There are no samples in the config file."
            #sys.exit(0)
 
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
        
        return am    
    elif filetype is "samplefile":
        samples=MakeSampleMap(filelines)
        return samples
    elif filetype is "branchlist":
        branches=MakeBranchMap(filelines)
        return branches
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

        for item in line.split():
            name,value=item.split("=")
            if name.find("name") is 0:
                samplename=str(value)
            if name.find("file") is 0:
                samplepaths.append(str(value))
            if name.find("dir") is 0:
                from os import listdir
                from os.path import isfile, join
                onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
                for rootfile in onlyfiles:
                    print rootfile
                    if rootfile.find(".root") != -1:
                        samplepaths.append(str(value)+"/"+str(rootfile))
            if name.find("type") is 0:
                sampletype=int(value)
            if name.find("xsec") is 0:
                samplexsec=float(value)
            if name.find("kfac") is 0:
                samplekfac=float(value)
            if name.find("scale") is 0:
                samplescale=float(value)
        
        if samplename != "":
            samples[samplename]=[sampletype,samplexsec,samplekfac,samplescale,samplepaths]
        else:
            print "sample name is empty",samplename,"not filling"
    
    return samples


def MakeBranchMap(lines):
    branches={}

    for line in lines:
        branchname=""
        branchtype=-1
        arraylength=-1
        
        for item in line.split():
            name,value = item.split("=")
            if name.find("name") is 0:
                branchname=value
            if name.find("type") is 0:
                branchtype=int(value)
            if name.find("max") is 0:
                arraylength=int(value)

        branches[branchname]= [branchtype,arraylength]

    return branches            

