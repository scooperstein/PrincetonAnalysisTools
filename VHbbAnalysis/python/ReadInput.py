#~ /usr/bin/python

debug=0

def ReadTextFile(filename, filetype):
    if debug > 1000:
         print "filetype is ", filetype
         print "filename is ", filename

    textfile=open(filename, 'r')

    filelines=textfile.readlines()
    filelines=CleanLines(filelines)

    print filelines
    if filetype is "samplefile":
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
    

def MakeSampleMap(lines):
    print "Reading samples"
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
                print "Need to put a function here that 'ls's the content of a dir and appends the list"
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

