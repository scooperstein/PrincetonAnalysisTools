import ROOT
import sys,time,math
import TdrStyles

#
#   This program takes a ROOT file containing
#   histograms with names of a particular format.
#   Various histograms with identical parameters
#   but containing different samples are plotted
#   together.
#
#   This program is based on a program written
#   by Sean Simon (UCSD).  It was re-written in 
#   python and made stand-alone by Chris Palmer
#   (Princeton).  Matteo Sani (UCSD) also 
#   contributed to the development of this 
#   software.
#


NFILES=-1
NIND=-1

plotinfos={}
samples={}

# Set Defaults
pads = []

xsize=1000
ysize=800
can=ROOT.TCanvas("can","can",10,10,xsize,ysize)

legx1=0.6
legx2=0.95
legy1=0.6
legy2=0.9
legend = ROOT.TLegend(legx1,legy1,legx2,legy2)
legend2= ROOT.TLegend(legx1,legy1,legx2,legy2)
nLeg=1
nLegEntries=0
dolegend=True

Normalize = False

textx1=0.3
textx2=0.6
texty1=0.8
texty2=0.9
plottext = ROOT.TPaveText(textx1,texty1,textx2,texty2,"brNDC");
dotext=True
text=""

Ncol=1
Nrow=1

stackmax=0
stackmin=0.01
stackminlog=0.01
maxscaleup=1.1
linewidth=2
NReBin=1
plotscale=ysize/700

cutline=ROOT.TLine()

Debug=True
DebugNew=False

dodivide=False
dosoverb=False
dosoversqrtb=False
dolog=False
dogridx=False
dogridy=False
dorebin=False
doreplot=True
doxtitle=True
doytitle=True
dooflow=False
douflow=False
doscale=True
docats=True
singalcat=-1
doline=False
linex=-2.0
domergecats=False
StaticMin=False
StaticMax=False
dotitles=False
dotdr=True
dochi2=True
#cattitles=["EB-EB-hiR9-hiR9","EB-EB-!(hiR9-hiR9)","!(EB-EB)-hiR9-hiR9","!(EB-EB)-!(hiR9-hiR9)"]
#cattitles=["Electron Barrel (Eta < 0.8)","Electron Barrel (Eta > 0.8)","Electron Endcaps"]
#cattitles=["|#eta| < 0.8","0.8 < |#eta| < 1.4442","|#eta| > 1.556"]
#cattitles = ["TT CS (Wl#nu)","W+HF CS (Wl#nu)","W+LF CS (Wl#nu)"]
cattitles = ["TT CS (W#mu#nu)","W+HF CS (W#mu#nu)","W+LF CS (W#mu#nu)","TT CS (We#nu)","W+HF CS (We#nu)", "W+LF CS (We#nu)"]
#cattitles=["W#mu#nu","We#nu"]
dointegrals=False
domean=False

allintegrals=False

datasampleIndex=0

linex=-1
dodata=True
dobkg=True
doautogroup=False

readinputfiles=True
DoFileSave=True
filesetup=True
outfilename="PI_output.root"
outfile=ROOT.TFile(outfilename,"recreate")


DoBigLegend=False
DoPopSig=True
DoData=True
DoFill=True
DoFillSig=True
DoCats=True
DoStack=True

DoMulti=False
DoWriteAll=False
DoPrintFlows=False
DoPrintBins=False
DoStats=False
DoSumw2=False
DoInt=False
DoRevInt=False




class PlotInfo:
    def __init__(self):
        self.ncat=1
        self.nbinsx=1
        self.nbinsy=1
        self.lowlim=0
        self.highlim=0
        self.lowlim2=0
        self.highlim2=0
        self.xaxislabel=""
        self.yaxislabel=""
        self.plotvarname=""
        self.plotvarcatname=""
        self.catid=0
        self.index=-1
        if Debug:
            "New PlotInfo"
        
    def Print(self):
        print "ncat",self.ncat
        print "nbinsx",self.nbinsx
        print "nbinsy",self.nbinsy
        print "lowlim",self.lowlim
        print "highlim",self.highlim
        print "lowlim2",self.lowlim2
        print "highlim2",self.highlim2
        print "xaxislabel",self.xaxislabel
        print "yaxislabel",self.yaxislabel
        print "plotvarname",self.plotvarname
        print "plotvarcatname",self.plotvarcatname
        print "catid",self.catid

    def ColsRows(self):
        if self.ncat==1:
            return 1,1
        elif self.ncat==2:
            return 2,1
        elif self.ncat==3:
            return 3,1
        elif self.ncat==4:
            return 2,2
        elif self.ncat==5 or self.ncat==6:
            return 3,2
        elif self.ncat==7 or self.ncat==8:
            return 4,2
        elif self.ncat==9:
            return 3,3
        else:
            return 4,5

        

class SampleInfo:
    def __init__(self):
        self.sampleIndex=-99999
        self.inshortnames=""
        self.plotsample=1
        self.addtoleg=1
        self.order=-1
        self.color=-1
        self.scale=1.0
        self.displayname=""
        self.configline=""

    def Print(self):
        print "sampleIndex",         self.sampleIndex
        print "inshortnames",  self.inshortnames
        print "plotsample",    self.plotsample
        print "addtoleg",      self.addtoleg
        print "order",         self.order
        print "color",         self.color
        print "configline",    self.configline
        

    def ParseConfigLines(self):
        if self.configline is not "":
            items=self.configline.split()
            for item in items:
                if item.find("sampleIndex=") is not -1:
                    continue
                elif item.find("order=") is not -1:
                    self.order=int(item.split("=")[1])
                elif item.find("plot=") is not -1:
                    self.plotsample=int(item.split("=")[1])
                elif item.find("leg=") is not -1:
                    self.addtoleg=int(item.split("=")[1])
                elif item.find("color=") is not -1:
                    self.color=int(item.split("=")[1])
                elif item.find("marker=") is not -1:
                    self.maker=int(item.split("=")[1])
                elif item.find("scale=") is not -1:
                    self.scale=float(item.split("=")[1])
                elif item.find("displayName=") is not -1:
                    self.displayname=ROOT.TString(item.split("=")[1])
                    self.displayname.ReplaceAll("@"," ")
                else:
                    print "Item",item,"not parsed."
        else:
            self.plotsample=0
            print "No configuration--no plotting"
   


NEWVALS=["datasampleIndex","linex","stackmax","stackmin","stackminlog","maxscaleup","linewidth","NReBin","plotscale","legx1","legx2","legy1","legy2","textx1","textx2","texty1","texty2","setminlog"]

PLOTPROPS=["dolog","dogridx","dogridy","doline","domergecats","dotdr","dochi2","doreplot","dointegrals","dotitles","doxtitle","doytitle","dooflow","douflow","dorebin","StaticMin","StaticMax","dolegend","dotext","dodata","dobkg","dodivide", "dosoverb", "dosoversqrtb", "Normalize" ,"Debug","DebugNew"]
def FindFunction(option):
    print option
    if option == "START":
        if len(sys.argv) <2:
            print "Please pass configuration"
            sys.exit(0)
        print "Running startup"
        Startup(sys.argv[1])
    elif option.split()[0] == "settext":
        ReSetText(float(option.split()[1]),
        float(option.split()[2]),
        float(option.split()[3]),
        float(option.split()[4]))
    elif option.split()[0] == "setleg":
        SetLeg(float(option.split()[1]),
        float(option.split()[2]),
        float(option.split()[3]),
        float(option.split()[4]))
    elif option == "ls":
        PrintPlotNames(plotinfos)
    elif option == "samples":
        PrintSamples(samples)
    elif option == "integrals":
        PrintTableOfIntegrals(samples)
    elif option == "groupintegrals":
        PrintTableOfIntegrals(samples,True)
    elif option == "vars":
        PrintAllVariables()
    elif option == "save":
        SaveHist()
    elif option.split()[0] == "linex":
        SetLine(float(option.split()[1]))
    elif option.split()[0] == "mvlegx":
        MoveLegX(float(option.split()[1]))
    elif option.split()[0] == "mvlegy":
        MoveLegY(float(option.split()[1]))
    elif option.split()[0] == "mvtextx":
        MoveTextX(float(option.split()[1]))
    elif option.split()[0] == "mvtexty":
        MoveTextY(float(option.split()[1]))
    elif option.split()[0] == "singlecat":
        SetCat(int(option.split()[1]))
    elif option.split()[0] == "text":
        NewText(option)
    elif option.split()[0] == "ploton":
        PlotSample(int(option.split()[1]))
    elif option.split()[0] == "j":
        Plot(option.split()[1])
    elif option.split()[0] == "cmds":
        ParseTextCommands(option.split()[1])
    elif option == "n":
        NextPlot()
    elif option == "p":
        PrevPlot()
    elif option=="allintegrals":
        GetSampleIntegrals()
    elif option.split()[0]=="rebin":
        SetReBin(option.split()[1])
    elif option.split()[0]=="maxval":
        SetMax(option.split()[1])
    elif option.split()[0]=="minval":
        SetMin(option.split()[1])
    elif option.split()[0]=="resize":
        Resize(option.split()[1],option.split()[2])
    elif option.split()[0] == "write":
        if len(option.split()) == 2:
            Write(str(option.split()[1]))
        else: 
            if option.split()[1]=="all":
                tempind=cur_plot.index
                for iplot in range(0,len(plotinfos)):
                    Plot(iplot)
                    Write(str(option.split()[2]))
                Plot(int(tempind))
            else:
                WriteCat(str(option.split()[1]),int(option.split()[2]))
    elif option in PLOTPROPS:
        Switch(option)
    elif len(option.split())==2:
        itemlist=option.split()
        NewValue(str(itemlist[0]),itemlist[1])
    elif option == "help":
        ListCMDs()
    else:
        print option, "not found"



def ListCMDs():
    print "START:               Read in plotvariables and inputfiles trees"
    print "ls:                  List plots"
    print "samples:             List samples"
    print "ploton sampleIndex:        Turns sample with sampleIndex on/off"
    print "resize x y           Resizes the canvas"
    print "j #:                 Standard plot of #th plot in plot list"
    print "p:                   Previous plot"
    print "n:                   Next plot"
    print "rebin #:             Rebin histograms by #"
    print "maxval #:            Set the maximum of the plot (0 turns off)"
    print "minval #:            Set the minimum of the plot (0 turns off)"
    print "setleg x1 x2 y1 y2   Moves position of the legend"
    print "settext x1 x2 y1 y2  Moves position of the text"
    print "variable num:        Change variable to num"
    print "vars:                Prints all global variables"
    print "singlecat #:         Display only # cat (-1 turns off)"
    print "mvlegx/y #:          Move legend by # in x/y"
    print "mvtextx/y #:         Move text by # in x/y"
    print "write appendix:      Write current plot with appendix"
    print "write all appendix:  Writes all plots with current settings"
    print "cmds:             Will loop over the lines of a file and"
    print "                     excute those lines as commands"
    print "variablename newval:",NEWVALS
    print "Switch on/off:      ",PLOTPROPS



### sample editing functions
def PlotSample(sampleIndex):
    if samples.has_key(sampleIndex):
        if samples[sampleIndex].plotsample == 1:
            samples[sampleIndex].plotsample=0
        else:
            samples[sampleIndex].plotsample=1
        ReDraw()
    else:
        print "no sample",sampleIndex
        PrintSamples(samples)
        

def PrintAllVariables():
    keys=globals().keys()
    for var in keys:
        if var == "samples":
            continue
        if var == "plotinfos":
            continue
        print "Name Val",var, globals()[var]


def ParseTextCommands(textfilename):
    try:
        textfile=open(str(textfilename))
    except:
        print "There is no file",textfilename
        return;

    for line in textfile:
        line=line.split("\n")[0]
        try:
            FindFunction(line)
        except:
            print "This line cannot be parsed:"
            print line



### Plot settings / Canvas Properties
def Switch(option):
    var = globals()[option]
    if Debug:
        print "var",var
    if var == True:
        globals()[option]=False
    else:
        globals()[option]=True
    ReDraw()

        
def NewValue(name,newval):
    if globals().has_key(name):
        globals()[name]=type(globals()[name])(newval)
        ReDraw()
    else:
        print "There is no",name,"to modify."

def SetReBin(nrebin):
    global dorebin
    dorebin = True
    global NReBin
    NReBin = int(nrebin)
    ReDraw()


def SetMin(minval):
    global StaticMin,stackmin
    stackmin=float(minval)
    if stackmin==0:
        StaticMin=False
    else:
        StaticMin=True
    ReDraw()


def SetMax(maxval):
    global StaticMax,stackmax
    stackmax=float(maxval)
    if stackmax==0:
        StaticMax=False
    else:
        StaticMax=True
    ReDraw()


def Resize(x,y):
    global xsize
    xsize = int(x)
    global ysize
    ysize = int(y)
    ReDraw()

def SetLeg(x1,x2,y1,y2):
    global dolegend,legx1,legx2,legy1,legy2
    legx1=x1
    legx2=x2
    legy1=y1
    legy2=y2
    dolegend=True
    ReDraw()


def ReSetText(x1,x2,y1,y2):
    global dotext,textx1,textx2,texty1,texty2
    textx1=x1
    textx2=x2
    texty1=y1
    texty2=y2
    dotext=True
    ReDraw()


def NewText(fulloption):
    global text
    fulloption=fulloption.lstrip("text")
    text=fulloption.lstrip()
    ReDraw()


#  Information

def PrintPlotNames(plotinfos):
    plotinds = plotinfos.keys()
    plotinds.sort()
    for ind in plotinds:
        print ind, plotinfos[ind].plotvarname

def PrintSamples(samples):
    sampleIndices = samples.keys()
    for sampleIndex in sampleIndices:
        print "Plotting:",samples[sampleIndex].plotsample," sampleIndex:",sampleIndex,"\t",samples[sampleIndex].inshortnames

def PrintTableOfIntegrals(samples,grouped=False):
    sampleIndices = samples.keys()
    nCats=cur_plot.ncat
    orderMap={}
    for sampleIndex in sampleIndices:
        orderMap[samples[sampleIndex].order]=sampleIndex
    ordered=orderMap.keys()
    ordered.sort()
    orderedSampleIndices=[]
    groupedSampleIndices=[]
#print "samples[sampleIndex].addtoleg",samples[sampleIndex].addtoleg
    for iOrder in ordered:
        if samples[orderMap[iOrder]].addtoleg:
            groupedSampleIndices.append(orderMap[iOrder])
        orderedSampleIndices.append(orderMap[iOrder])

    table={}
    if grouped:
        for sampleIndex in groupedSampleIndices:
            table[samples[sampleIndex].inshortnames]=[0]*nCats
    else:
        for sampleIndex in orderedSampleIndices:
            table[samples[sampleIndex].inshortnames]=[0]*nCats
    
    iGroup=0
    for sampleIndex in orderedSampleIndices:
        for iCat in range(nCats):
            try:
                histname=str(cur_plot.plotvarname)+"_cat"+str(iCat)+"_"+str(samples[sampleIndex].inshortnames)
                hist=ROOT.gROOT.FindObject(histname).Clone()
                if grouped:
                    table[samples[groupedSampleIndices[iGroup]].inshortnames][iCat]=table[samples[groupedSampleIndices[iGroup]].inshortnames][iCat]+hist.Integral()*samples[sampleIndex].scale
                else:
                    table[samples[sampleIndex].inshortnames][iCat]=table[samples[sampleIndex].inshortnames][iCat]+hist.Integral()*samples[sampleIndex].scale
            except:
                print "Can't sum for sample",samples[sampleIndex].inshortnames,"in cat",iCat
        if sampleIndex==groupedSampleIndices[iGroup]:
            iGroup=iGroup+1

    if grouped:
        orderedListForNames=groupedSampleIndices
    else:
        orderedListForNames=orderedSampleIndices

    for index in orderedListForNames:
        name=samples[index].inshortnames
        if grouped:
            print repr(samples[index].displayname).rjust(25),
        else:
            print repr(name).rjust(25),
        for iCat in range(len(table[name])):
            try:
                #print '{:.2e}'.format(table[name][iCat]),
                print '{:.5e}'.format(table[name][iCat]),
            except:
                print '------  ',
    
        print


#  Setup Functions

def Startup(configfilename):
    print "Reading in inputfiles and plotvariables from rootfile"
    print "  and reading plot features from "+configfilename
    file=open(configfilename, "r")
    configlines=file.readlines()
    file.close()
    #print configlines
    for line in configlines:
        if Debug:
            print "line",line
        if line.find("file=") is not -1:
            rootline=(line.split("file=")[1]).split()[0]
            global rootfile
            rootfile=ROOT.TFile(rootline)
            break
    SetupGlobalVariables(configlines)
    print "rootline",rootline
    if rootfile is not "":
        ReadPlotVariables(rootfile)
        if readinputfiles:
            ReadInputFiles(rootfile)
        
    else:
        print "need root file in config.dat"
        sys.exit(0)
    
    AssociateConfigToSample(configlines)

    for sampleIndex in samples:
        samples[sampleIndex].ParseConfigLines()

    #if DoFileSave:
    #    SetupSaveFile(outfilename)

def CountSamplesInLegend():
    global nLegEntries
    nLegEntries=0
    for sampleIndex in samples:
        if samples[sampleIndex].plotsample==1 and samples[sampleIndex].addtoleg==1:
            nLegEntries=nLegEntries+1
    

def AssociateConfigToSample(configlines):
    for line in configlines:
        if line.find("name=") is not -1:
            if Debug:
                print "AssociateConfigToSample",line
            items=line.split()
            for item in items:
                print item
                if item.find("sampleIndex=") is not -1:
                    sampleIndex_config=int(item.split("=")[1])
                    if samples.has_key(sampleIndex_config) == False:
                        if DebugNew:
                            print "add new sample with name and sampleIndex from line"
                        newsample=SampleInfo()
                        newsample.sampleIndex         = sampleIndex_config
                        if DebugNew:
                            print "find the new name"
                        name=""
                        for item2 in items:
                            print item2
                            if item2.find("name=") is not -1:
                                name=str(item2.split("=")[1])
                                if DebugNew:
                                    print "found name",name
                                break
                        newsample.inshortnames  = str(name)
                        samples[sampleIndex_config]=newsample
                    if samples[sampleIndex_config].configline=="":
                        samples[sampleIndex_config].configline=line
                    elif DebugNew:
                        print "Found another input line for",sampleIndex_config
                        print "line",line


    
def SetupGlobalVariables(configlines):
    vars2check=globals().keys()
    for line in configlines:
        if line.find("name=") is not -1:
            continue # this is a files line
        if line.find("text=") is not -1:
            global text
            text=line.lstrip("text=")
            text=text.rstrip("\n")
            print line
            print text
            continue 
        if Debug:
            print "SetupGlobalVariables",line
        for item in line.split():
            if Debug:
                print item
            itemlist=item.split("=")
            if itemlist[0] in vars2check:
                if Debug:
                    print globals()[itemlist[0]]
                    print itemlist,type(itemlist[1])
                if itemlist[1]=='True':
                    globals()[itemlist[0]] = True
                elif itemlist[1]=='False':
                    globals()[itemlist[0]] = False
                else:
                    globals()[itemlist[0]] = type(globals()[itemlist[0]])(itemlist[1])
                if Debug:
                    print globals()[itemlist[0]]


def ReadPlotVariables(rootfile):
    rootfile.cd()
    plotvariables=ROOT.gDirectory.Get("plotvariables")
    entries=plotvariables.GetEntriesFast()
    pvnames=ROOT.TClonesArray("TObjString",entries)
    xaxesnames=ROOT.TClonesArray("TObjString",entries)
    yaxesnames=ROOT.TClonesArray("TObjString",entries)
    
    plotvariables.SetBranchAddress("plotvarnames",ROOT.AddressOf(pvnames))
    plotvariables.SetBranchAddress("xaxislabels",ROOT.AddressOf(xaxesnames))
    plotvariables.SetBranchAddress("yaxislabels",ROOT.AddressOf(yaxesnames))
     
    ientry = plotvariables.LoadTree(0)
    plotvariables.GetEntry(0)
    
    if Debug:
        print "plotvariables.Nvar",plotvariables.Nvar
    for i in range(0,plotvariables.Nvar):
        newplotinfo=PlotInfo()
   
        newplotinfo.index=i 
        newplotinfo.ncat=plotvariables.histoncat[i]
        newplotinfo.nbinsx=plotvariables.nbinsx[i]
        newplotinfo.nbinsy=plotvariables.nbinsy[i]
        newplotinfo.lowlim=plotvariables.lowlim[i]
        newplotinfo.highlim=plotvariables.highlim[i]
        newplotinfo.lowlim2=plotvariables.lowlim2[i]
        newplotinfo.highlim2=plotvariables.highlim2[i]
        newplotinfo.xaxislabel=plotvariables.xaxislabels[i].GetString()
        newplotinfo.xaxislabel.ReplaceAll("@"," ")
        newplotinfo.yaxislabel=plotvariables.yaxislabels[i].GetString()
        newplotinfo.yaxislabel.ReplaceAll("@"," ")
        newplotinfo.plotvarname=plotvariables.plotvarnames[i].GetString()
        #newplotinfo.plotvarcatname=
        #newplotinfo.catid=0
        
        plotinfos[i]=newplotinfo
        if Debug:
            plotinfos[i].Print()
        #newplotinfo.Print()


def ReadInputFiles(rootfile):
    rootfile.cd()
    inputfiles=ROOT.gDirectory.Get("inputfiles")
    entries=inputfiles.GetEntriesFast()
    shortnames=ROOT.TClonesArray("TObjString",entries)
    inputfiles.SetBranchAddress("inshortnames",ROOT.AddressOf(shortnames))
    
    ientry = inputfiles.LoadTree(0)
    inputfiles.GetEntry(0)

    NFILES      = inputfiles.nfiles
    NIND        = inputfiles.nindfiles
    intlumi     = inputfiles.intlumi

    for ifile in range(NIND):
        newsample=SampleInfo()
        newsample.sampleIndex           =inputfiles.sampleIndex[ifile]
        newsample.inshortnames    =inputfiles.inshortnames[ifile].GetString()
        samples[inputfiles.sampleIndex[ifile]]=newsample


#def SetupSaveFile(outfilename):
#    global filesetup
#    outfilename=outfilename.split(".")
#    if not filesetup:
#        if globals().has_key(outfilename):
#            print "Output file already created but filesetup is",filesetup
#            filesetup=True
#        else:
#            try:
#                globals()[outfilename]=ROOT.TFile(outfilename,"recreate")
#                filesetup=True
#            except:
#                print "Tried to setup file... fail"
#    else:
#        print "File is already setup",filesetup,outfilename


def SaveHist():
    #global filesetup
    #if not filesetup:
    #    global outfilename
    #    SetupSaveFile(outfilename)
   
    print outfile
    outfile.cd()
    print "cd'ed"
    cur_plot.Write()

### Ploting Functions 

def Plot(num,printsuffix="",printcat=-1):
    global cur_plot
    cur_plot=plotinfos[int(num)]
    print num, cur_plot.plotvarname
    if docats:
        Ncol,Nrow=cur_plot.ColsRows()
    else:
        Ncol = 1
        Nrow = 1

    can.Clear()
    can.SetWindowSize(xsize,ysize)
    can.SetCanvasSize(xsize-4,ysize-28)
    
    if dotdr:
        TdrStyles.tdrStyle()   
 
    if domergecats:
        can.cd()
        can.cd().SetLogy(dolog)
        can.cd().SetGridx(dogridx)
        can.cd().SetGridy(dogridy)
        
    else:
        nCanvases = Ncol*Nrow
        if (not dodivide and not dosoverb and not dosoversqrtb):
            can.Divide(Ncol, Nrow)
            for ican in range(1, nCanvases+1):
                can.cd(ican)
                can.cd(ican).SetLogy(dolog)
                can.cd(ican).SetGridx(dogridx)
                can.cd(ican).SetGridy(dogridy) 
                can.cd(ican).SetTopMargin(0.01)
                can.cd(ican).SetBottomMargin(0.1)
                can.cd(ican).SetRightMargin(0.01)
                can.cd(ican).SetLeftMargin(0.07)
        else:
            global pads
            pads = []
            nCanvases = Ncol*Nrow*2
            nx = Ncol
            ny = Nrow*2
            xmargin=0.0
            ymargin=0.0 
            dy = 1./float(Nrow)
            dx = 1./float(nx)
            dy_odd  = dy*.3
            dy_even = dy*.7
            for iy in range(ny):
                if iy%2==0:
                    y1 = 1 - (iy/2)*dy - ymargin
                    y2 = 1 - (iy/2)*dy - dy_even + ymargin 
                else:
                    y1 = 1 - (iy/2)*dy - dy_even - ymargin
                    y2 = 1 - (iy/2+1)*dy + ymargin
                for ix in range(nx):
                    x1 = ix*dx + xmargin
                    x2 = x1 +dx -2*xmargin
                    if (x1 > x2):
                        continue
                    name = can.GetName()+"_"+str(ix+iy*nx)
                    pad = ROOT.TPad(name, name, x1, y1, x2, y2, ROOT.kWhite)
                    #pad.SetTopMargin(0.07)
                    #pad.SetBottomMargin(0.1)
                    #pad.SetRightMargin(0.01)
                    #pad.SetLeftMargin(0.07) 
                    pad.SetFrameBorderMode(0)
                    if iy%2==0:
                        pad.SetTopMargin(0.05);
                        pad.SetBottomMargin(0.0)
                    else:
                        pad.SetTopMargin(0.00)
                        pad.SetBottomMargin(0.3)
                    pad.SetRightMargin(0.05)
                    pad.SetLeftMargin(0.12) 
                    pad.SetFillStyle(4000)
                    pad.SetFrameFillStyle(1000)
                    #pad.SetRightMargin(0.0)
                    #pad.SetLeftMargin(0.0) 
                    pads.append(pad)
                    pads[-1].SetNumber(ix+iy*nx+1)
                    pads[-1].Draw()

            for ican in range(1, nCanvases+1):
                if (((ican-1)/Ncol)%2 == 0):
                    can.cd(ican)
                    can.cd(ican).SetLogy(dolog)
                    can.cd(ican).SetGridx(dogridx)
                    can.cd(ican).SetGridy(dogridy)
                
    stacks={}
    if doautogroup:
        stacktypes=AutoGroupNames()
    else:
        if dodata:
            if dobkg:
                stacktypes=["bkg","sig","data"]
            else:
                stacktypes=["sig","data"]
        else:
            if dobkg:
                stacktypes=["bkg","sig"]
            else:
                stacktypes=["sig"]


    first=1
    stackmaxima={}
    plottexttitle = [""]*cur_plot.ncat
    if docats:
        cats=xrange(cur_plot.ncat)
    else:
        cats=[]
        cats.append(singlecat)
    
    if domergecats:
        stackmaxima[0]=[]
        for stacktype in stacktypes:
            if DebugNew:
                print "stacktype",stacktype
                print "MakeMergeStack"
            stacks[stacktype]=MakeMergeStack(stacktype)
            if DebugNew:
                print "MakeMergeOverlayLines"
            stacks[stacktype+"lines"]=MakeMergeOverlayLines(stacktype)
            stackmaxima[0].append(stacks[stacktype].GetMaximum())          
    else:
        for icat in cats:
            stackmaxima[icat]=[]
            for stacktype in stacktypes:
                stacks[stacktype+str(icat)]=MakeStack(stacktype,icat)
                stacks[stacktype+"lines"+str(icat)]=MakeOverlayLines(stacktype, icat)
                stackmaxima[icat].append(stacks[stacktype+str(icat)].GetMaximum())
            
    if dolegend:
        CountSamplesInLegend()
        SetLegend()

    if dotext:
        SetText()

    if dointegrals:
        stackintegrals={}


    if domergecats:
        stackmaxima[0].sort()
        
        if StaticMax:
            global stackmax
        else:
            if dolog is True:
                logmax=math.log10(stackmaxima[0][-1])
                logmin=math.log10(stackminlog)
                stackmax=math.pow(10,maxscaleup*logmax - (maxscaleup-1)*logmin)
            else:
                stackmax=stackmaxima[0][-1]*maxscaleup
        if DebugNew:
            print "mergecats stackmaxima",stackmaxima
            print "mergecats stackmax",stackmax

      
        stacks["bkg"].Draw("hist")
        ROOT.gStyle.SetOptStat(0)
        
        stacks["bkg"].SetMaximum(stackmax)
            
        if StaticMin:
            stacks["bkg"].SetMinimum(stackmin)

        if doxtitle==True:
            stacks["bkg"].GetXaxis().SetTitle(str(cur_plot.xaxislabel))
        else:
            stacks["bkg"].GetXaxis().SetTitle("")
        
        if doytitle==True:
            stacks["bkg"].GetYaxis().SetTitle(str(cur_plot.yaxislabel))
        else:
            stacks["bkg"].GetYaxis().SetTitle("")
            
        
        iLeg=0.0
        if dodivide:
            stacks["bkg"].GetXaxis().SetTitleSize(0.0)
            stacks["bkg"].GetXaxis().SetTitleOffset(999)
        else:
            stacks["bkg"].GetXaxis().SetTitleSize(0.04)
        stacks["bkg"].GetYaxis().SetTitleSize(0.13)
        stacks["bkg"].GetYaxis().SetTitleOffset(0.8)
        stacks["bkg"].SetTitle("All Categories")

        dataIntegral = -1
        if Debug:
            print "finished start"
       
        if dodata: 
            lineorder = stacks["datalines"].keys()
            lineorder.sort()
            lineorder.reverse()
            if Debug:
                print "lineorder",lineorder
            for index in lineorder:
                sampleIndex=index[1]
                if dolegend and nLeg==1:
                    legend.AddEntry(stacks["datalines"][index],str(samples[sampleIndex].displayname),"ep"); 
                elif dolegend and nLeg==2:
                    if iLeg/float(nLegEntries) <0.5:
                        legend.AddEntry(stacks["datalines"][index],str(samples[sampleIndex].displayname),"ep"); 
                    else:
                        legend2.AddEntry(stacks["datalines"][index],str(samples[sampleIndex].displayname),"ep"); 
                    iLeg=iLeg+1.0
            #dataIntegral = stacks["datalines"+str(icat)][lineorder[0]].Integral()
            print dataIntegral
        if DebugNew:
            print "finished data"
        
        lineorder = stacks["siglines"].keys()
        lineorder.sort()
        lineorder.reverse()
        for index in lineorder:
            sampleIndex=index[1]
            if dolegend and nLeg==1:
                legend.AddEntry(stacks["siglines"][index],str(samples[sampleIndex].displayname),"l"); 
            elif dolegend and nLeg==2:
                if iLeg/float(nLegEntries) <0.5:
                    legend.AddEntry(stacks["siglines"][index],str(samples[sampleIndex].displayname),"l"); 
                else:
                    legend2.AddEntry(stacks["siglines"][index],str(samples[sampleIndex].displayname),"l"); 
                iLeg=iLeg+1.0
            #sigIntegral = stacks["siglines"+str(icat)][index].Integral()
            #print sigIntegral
        
        lineorder = stacks["bkglines"].keys()
        lineorder.sort()
        lineorder.reverse()
        if DebugNew:
            print "lineorder",lineorder
        first=1
        if (Normalize):
            if stacks["bkglines"+str(icat)][lineorder[0]].Integral() > 0:
                toscale=dataIntegral/stacks["bkglines"+str(icat)][lineorder[0]].Integral()
            else:
                toscale=1.0
        else:
            toscale=1.0
        for index in lineorder:
            sampleIndex=index[1]
            if DebugNew:
                print "index",index
                print "sampleIndex",sampleIndex
                print "samples[sampleIndex].color",samples[sampleIndex].color
            stacks["bkglines"][index].SetLineColor(ROOT.kBlack)
            stacks["bkglines"][index].SetLineWidth(linewidth*plotscale)
            stacks["bkglines"][index].SetFillStyle(1001)
            stacks["bkglines"][index].SetFillColor(int(samples[sampleIndex].color))
            stacks["bkglines"][index].Scale(toscale)
            stacks["bkglines"][index].Draw("histsame")
            if dolegend and nLeg==1:
                legend.AddEntry(stacks["bkglines"][index],str(samples[sampleIndex].displayname),"f"); 
            elif dolegend and nLeg==2:
                if iLeg/float(nLegEntries) <0.5:
                    legend.AddEntry(stacks["bkglines"][index],str(samples[sampleIndex].displayname),"f"); 
                else:
                    legend2.AddEntry(stacks["bkglines"][index],str(samples[sampleIndex].displayname),"f"); 
                iLeg=iLeg+1.0
        #print "bkg integral",stacks["bkglines"+str(icat)][lineorder[0]].Integral()
        
        if dointegrals:
            oflowbin = int(stacks["bkglines"][lineorder[0]].GetNbinsX()+1)
            stackintegrals["bkg"]=stacks["bkglines"][lineorder[0]].Integral(0,oflowbin)
        
        lineorder = stacks["siglines"].keys()
        lineorder.sort()
        lineorder.reverse()
        for index in lineorder:
            sampleIndex=index[1]
            stacks["siglines"][index].Draw("histsame")
            #stacks["siglines"+str(icat)][index].Draw("histsame")
       
        if dointegrals:
            oflowbin = int(stacks["siglines"][lineorder[0]].GetNbinsX()+1)
            stackintegrals["sig"]=stacks["siglines"][lineorder[0]].Integral(0,oflowbin)
       
        
        if dodata: 
            lineorder = stacks["datalines"].keys()
            lineorder.sort()
            lineorder.reverse()
            for index in lineorder:
                sampleIndex=index[1]
                stacks["datalines"][index].Draw("epsame")
                #stacks["datalines"+str(icat)][index].Draw("epsame")
            
            if dointegrals:
                oflowbin = int(stacks["datalines"][lineorder[0]].GetNbinsX()+1)
                stackintegrals["data"]=stacks["datalines"][lineorder[0]].Integral(0,oflowbin)
        
        if dolegend:
            if nLeg==1:
                legend.Draw()
            elif nLeg==2:
                legend.Draw()
                legend2.Draw()

        if dotext:
            plottext.Draw()
            if dotitles:
                plottexttitle[icat] = SetTextTitle() # have to refresh to get rid of old category title
                plottexttitle[icat].AddText(cattitles[icat])
                plottexttitle[icat].Draw()
    
        if doline:
            cutline.SetLineWidth(4*plotscale)
            cutline.SetLineColor(633)
            cutline.DrawLine(linex,0.0,linex,stackmax/1.1);
            
    else:   #  this is the typical plot, domergecats plots all cats together
        
        hist_1 = [""]*cur_plot.ncat
        dataTot = [""]*cur_plot.ncat
        mcTot = [""]*cur_plot.ncat
        mcErrBand = [""]*cur_plot.ncat
        sigTot = [""]*cur_plot.ncat
        chi2text = [""]*cur_plot.ncat
        for icat in cats:
            stackmaxima[icat].sort()
            if StaticMax:
                global stackmax
            else:
                if dolog is True:
                    logmax=math.log10(stackmaxima[icat][-1])
                    logmin=math.log10(stackminlog)
                    stackmax=math.pow(10,maxscaleup*logmax - (maxscaleup-1)*logmin)
                else:
                    stackmax=stackmaxima[icat][-1]*maxscaleup
                if stackmax==0:
                    stackmax = 1
            if DebugNew:
                print "docats stackmaxima",stackmaxima
                print "docats stackmax",stackmax

            if docats:
                if (not dodivide and not dosoverb and not dosoversqrtb):
                    if DebugNew:
                        print "don't divide"
                    can.cd(icat+1)
                else:
                    pad_id = icat%Ncol + 1 + ((icat)/Ncol) * 2 * Ncol # (icat%Ncol)+1 + ((icat/Ncol)+1)*Ncol
                    print "icat, pad_id,",icat, pad_id
                    can.cd(pad_id)
       
            # Matteo no stack is plotted to allow normalization (and stats are not plot as well)
            ROOT.gStyle.SetOptStat(0)
            initalstack=""
            if dobkg:
                initalstack="bkg"
            elif dodata:
                initalstack="data"
            else:
                initalstack="sig"

            stacks[initalstack+str(icat)].Draw("hist")
            stacks[initalstack+str(icat)].SetMaximum(stackmax)
            if DebugNew:
                print "initial stack drawn/setmax"

            if StaticMin:
                stacks[initalstack+str(icat)].SetMinimum(stackmin)

            if dolog==True:
                stacks[initalstack+str(icat)].SetMinimum(stackminlog)

            if doxtitle==True:
                stacks[initalstack+str(icat)].GetXaxis().SetTitle(str(cur_plot.xaxislabel))
            else:
                stacks[initalstack+str(icat)].GetXaxis().SetTitle("")
            
            if doytitle==True:
                stacks[initalstack+str(icat)].GetYaxis().SetTitle(str(cur_plot.yaxislabel))
            else:
                stacks[initalstack+str(icat)].GetYaxis().SetTitle("")
            
            if not dodivide:              
                stacks[initalstack+str(icat)].GetXaxis().SetTitleSize(0.045)
            else:
                stacks[initalstack+str(icat)].GetXaxis().SetTitleSize(0.0)
                stacks[initalstack+str(icat)].GetXaxis().SetTitleOffset(999)

            stacks[initalstack+str(icat)].GetYaxis().SetTitleSize(0.13)
            stacks[initalstack+str(icat)].GetYaxis().SetTitleOffset(0.8)
            
            if dotitles:
                if len(cattitles)==int(cur_plot.ncat):
                    stacks[initalstack+str(icat)].SetTitle(str(cattitles[icat]))
                else:
                    stacks[initalstack+str(icat)].SetTitle(str(cur_plot.plotvarname))
            else:
                stacks[initalstack+str(icat)].SetTitle("")
            
            if DebugNew:
                print "initial stack done"

            iLeg=0
            dataIntegral = -1
            if dodata: 
                lineorder = stacks["datalines"+str(icat)].keys()
                lineorder.sort()
                lineorder.reverse()
                for index in lineorder:
                    sampleIndex=index[1]
                    if icat==cats[0]:
                        if dolegend and nLeg==1:
                            legend.AddEntry(stacks["datalines"+str(icat)][index],str(samples[sampleIndex].displayname),"ep"); 
                        elif dolegend and nLeg==2:
                            if iLeg/float(nLegEntries) <0.5:
                                legend.AddEntry(stacks["datalines"+str(icat)][index],str(samples[sampleIndex].displayname),"ep"); 
                            else:
                                legend2.AddEntry(stacks["datalines"+str(icat)][index],str(samples[sampleIndex].displayname),"ep"); 
                        iLeg=iLeg+1.0
                dataIntegral = stacks["datalines"+str(icat)][lineorder[0]].Integral()
                print "data int ", dataIntegral
                if domean:
                    print "data mean","%.3f"%stacks["datalines"+str(icat)][lineorder[0]].GetMean()
                if dodivide:
                    dataTot[icat] =  stacks["datalines"+str(icat)][lineorder[0]].Clone("dataTot")
 
            lineorder = stacks["siglines"+str(icat)].keys()
            lineorder.sort()
            lineorder.reverse()
            for index in lineorder:
                sampleIndex=index[1]
                if icat==cats[0]:
                    if dolegend and nLeg==1:
                        legend.AddEntry(stacks["siglines"+str(icat)][index],str(samples[sampleIndex].displayname),"l"); 
                    elif dolegend and nLeg==2:
                        if iLeg/float(nLegEntries) <0.5:
                            legend.AddEntry(stacks["siglines"+str(icat)][index],str(samples[sampleIndex].displayname),"l"); 
                        else:
                            legend2.AddEntry(stacks["siglines"+str(icat)][index],str(samples[sampleIndex].displayname),"l"); 
                        iLeg=iLeg+1.0
                sigIntegral = stacks["siglines"+str(icat)][index].Integral()
                print "sig int ", sigIntegral
                if dosoverb or dosoversqrtb:
                    sigTot[icat] =  stacks["siglines"+str(icat)][lineorder[0]].Clone("sigTot")
 
            if dobkg:
                lineorder = stacks["bkglines"+str(icat)].keys()
                lineorder.sort()
                lineorder.reverse()
                if DebugNew:
                    print "lineorder",lineorder
                first=1
                if (Normalize):
                    if stacks["bkglines"+str(icat)][lineorder[0]].Integral() > 0:
                        toscale=dataIntegral/stacks["bkglines"+str(icat)][lineorder[0]].Integral()
                    else:
                        toscale=1.0
                else:
                    toscale=1.0
                for index in lineorder:
                    sampleIndex=index[1]
                    if DebugNew:
                        print "index",index
                        print "sampleIndex",sampleIndex
                        print "samples[sampleIndex].color",samples[sampleIndex].color
                    stacks["bkglines"+str(icat)][index].SetLineColor(ROOT.kBlack)
                    stacks["bkglines"+str(icat)][index].SetLineWidth(linewidth*plotscale)
                    stacks["bkglines"+str(icat)][index].SetFillStyle(1001)
                    stacks["bkglines"+str(icat)][index].SetFillColor(int(samples[sampleIndex].color))
                    stacks["bkglines"+str(icat)][index].Scale(toscale)
                    stacks["bkglines"+str(icat)][index].Draw("histsame")
                    if icat==cats[0]:
                        if dolegend and nLeg==1:
                            legend.AddEntry(stacks["bkglines"+str(icat)][index],str(samples[sampleIndex].displayname),"f"); 
                        elif dolegend and nLeg==2:
                            if iLeg/float(nLegEntries) <0.5:
                                legend.AddEntry(stacks["bkglines"+str(icat)][index],str(samples[sampleIndex].displayname),"f"); 
                            else:
                                legend2.AddEntry(stacks["bkglines"+str(icat)][index],str(samples[sampleIndex].displayname),"f"); 
                            iLeg=iLeg+1.0
                    if dodivide or dosoverb or dosoversqrtb:
                        if (index == lineorder[0]):
                            mcTot[icat] = stacks["bkglines"+str(icat)][index].Clone("mcTot")
                            if dosoversqrtb:
                                for bin in range(0,mcTot[icat].GetNbinsX()+1):
                                    mcTot[icat].SetBinContent(bin,math.sqrt(mcTot[icat].GetBinContent(bin)))
                # Draw MC statistical uncertainty 
                mcTot[icat].SetFillStyle(3013)
                mcTot[icat].SetFillColor(ROOT.kBlack)
                mcTot[icat].Draw("E2 same")
                if stacks["bkglines"+str(icat)][lineorder[0]].GetEffectiveEntries() > 0:
                    error=math.sqrt(1/float(stacks["bkglines"+str(icat)][lineorder[0]].GetEffectiveEntries()))*stacks["bkglines"+str(icat)][lineorder[0]].Integral() 
                    print "bkg int","%.2f"%stacks["bkglines"+str(icat)][lineorder[0]].Integral(),"+/-","%.2f"%error
                if domean:
                    print "bkg mean","%.3f"%stacks["bkglines"+str(icat)][lineorder[0]].GetMean()
 
                if dointegrals:
                    oflowbin = int(stacks["bkglines"+str(icat)][lineorder[0]].GetNbinsX()+1)
                    stackintegrals["bkg"+str(icat)]=stacks["bkglines"+str(icat)][lineorder[0]].Integral(0,oflowbin)
                    if DebugNew:
                        print "stack integrals"
                        print "last lineorder",lineorder[0]
                        print "bkg"+str(icat),stackintegrals["bkg"+str(icat)]
            
            lineorder = stacks["siglines"+str(icat)].keys()
            lineorder.sort()
            lineorder.reverse()
            for index in lineorder:
                sampleIndex=index[1]
                stacks["siglines"+str(icat)][index].Draw("histsame")
                if dodivide:
                    if (index == lineorder[0] and mcTot[icat] == ""):
                        mcTot[icat] = stacks["siglines"+str(icat)][index].Clone("mcTot")
            
            if dointegrals:
                oflowbin = int(stacks["siglines"+str(icat)][lineorder[0]].GetNbinsX()+1)
                stackintegrals["sig"+str(icat)]=stacks["siglines"+str(icat)][lineorder[0]].Integral(0,oflowbin)
 
            if dodata:
                lineorder = stacks["datalines"+str(icat)].keys()
                lineorder.sort()
                lineorder.reverse()
                for index in lineorder:
                    sampleIndex=index[1]
                    stacks["datalines"+str(icat)][index].Draw("epsame")
                
                if dointegrals:
                    oflowbin = int(stacks["datalines"+str(icat)][lineorder[0]].GetNbinsX()+1)
                    stackintegrals["data"+str(icat)]=stacks["datalines"+str(icat)][lineorder[0]].Integral(0,oflowbin)
                    if DebugNew:
                        print "stack integrals"
                        print "last lineorder",lineorder[-1]
                        print "data"+str(icat),stackintegrals["data"+str(icat)]
       
                if dochi2:
                    chi2Score = dataTot[icat].Chi2Test( mcTot[icat] , "UWCHI2/NDF")
                    print chi2Score
                    chi2text[icat] = ROOT.TPaveText(textx1,0.8,textx1+0.2,0.9,"brNDC");
                    chi2text[icat].SetTextFont(62);
                    chi2text[icat].SetTextSize(0.12)
                    chi2text[icat].SetBorderSize(0)
                    chi2text[icat].SetLineColor(0)
                    chi2text[icat].SetLineStyle(0)
                    chi2text[icat].SetLineWidth(0)
                    chi2text[icat].SetFillColor(0)
                    chi2text[icat].SetFillStyle(0)
                    chi2text[icat].AddText("#chi^{2}_{ }#lower[0.1]{/^{}#it{dof} = %.2f}"%(chi2Score))
                    if dodivide or dosoverb or dosoversqrtb:
                        if docats:
                            #    position in row  +  positions in prev rows
                            print "position in row, prev prositions",(icat%Ncol)+1,(icat/Ncol)*Ncol * 2 + Ncol
                            iPad = (icat%Ncol)+1  +  (icat/Ncol)*Ncol * 2 + Ncol
                            can.cd(iPad)
                        else:
                            can.cd(2)
                    chi2text[icat].Draw()
 
            if dodivide:
                if docats:
                    pad_id = icat%Ncol + 1 + (((icat)/Ncol) * 2 + 1) * Ncol # (icat%Ncol)+1 + ((icat/Ncol)+1)*Ncol
                else:
                    pad_id = 2
                can.cd(pad_id)
                can.GetPad(pad_id).SetGrid(True)
                can.GetPad(pad_id).SetGrid(True)
                dataTot[icat].Sumw2()
                mcTot[icat].Sumw2()
                dataTot[icat].Divide(mcTot[icat])
                dataTot[icat].SetMarkerColor(1)
                if docats:
                    dataTot[icat].SetMarkerSize(0.5)
                else:
                    dataTot[icat].SetMarkerSize(0.8)
                dataTot[icat].SetMarkerStyle(20)
                dataTot[icat].SetLineColor(1)
                dataTot[icat].SetLineWidth(2)
                dataTot[icat].SetMaximum(1.57)
                dataTot[icat].SetMinimum(0.5)
                dataTot[icat].GetYaxis().SetNdivisions(505)
                dataTot[icat].SetTitle('')
                dataTot[icat].GetXaxis().SetTitleOffset(1.)
                dataTot[icat].GetXaxis().SetTitleSize(0.13)
                dataTot[icat].GetXaxis().SetLabelOffset(0.014)
                dataTot[icat].GetXaxis().SetLabelSize(0.13)
                dataTot[icat].GetYaxis().SetTitleOffset(0.3)
                dataTot[icat].GetYaxis().SetTitleSize(0.13)
                dataTot[icat].GetYaxis().SetLabelOffset(0.005)
                dataTot[icat].GetYaxis().SetLabelSize(0.13)
                #dataTot[icat].GetYaxis().SetTitle('Data/MC')
                # FIXME need to make the font larger... how?
                dataTot[icat].GetYaxis().SetTitle('Data/MC')
                dataTot[icat].Draw('PE')

                # draw MC statistical uncertainty band
                mcErrBand[icat] = ROOT.TGraphErrors()
                for ibin in range(1,mcTot[icat].GetNbinsX()+1):
                    mcErrBand[icat].SetPoint(ibin,mcTot[icat].GetBinLowEdge(ibin)+(mcTot[icat].GetBinWidth(ibin)/2),1.)
                    #mcErrBand[icat].SetPoint(ibin,mcTot[icat].GetBinLowEdge(ibin),1.)
                    e = 0.0
                    if mcTot[icat].GetBinContent(ibin)>0:
                        e = mcTot[icat].GetBinError(ibin) / mcTot[icat].GetBinContent(ibin)
                    mcErrBand[icat].SetPointError(ibin,mcTot[icat].GetBinWidth(ibin)/2, e)
                mcErrBand[icat].SetFillColor(ROOT.kBlack)
                mcErrBand[icat].SetFillStyle(3013)
                mcErrBand[icat].Draw("SAME2") 
                if dochi2:    
                    chi2text[icat].Draw()

            if dosoverb or dosoversqrtb:
                pad_id = icat%Ncol + 1 + (((icat)/Ncol) * 2 + 1) * Ncol # (icat%Ncol)+1 + ((icat/Ncol)+1)*Ncol
                can.cd(pad_id)
                can.GetPad(pad_id).SetGrid(True)
                sigTot[icat].Sumw2()
                mcTot[icat].Sumw2()
                sigTot[icat].Divide(mcTot[icat])
                sigTot[icat].SetMarkerColor(4)
                if docats:
                    sigTot[icat].SetMarkerSize(0.5)
                else:
                    sigTot[icat].SetMarkerSize(0.8)
                sigTot[icat].SetMarkerStyle(20)
                sigTot[icat].SetLineColor(4)
                sigTot[icat].SetLineWidth(2)
                #sigTot[icat].SetMaximum(1.4)
                #sigTot[icat].SetMinimum(0.6)
                sigTot[icat].GetYaxis().SetNdivisions(505)
                sigTot[icat].SetTitle('')
                sigTot[icat].GetYaxis().SetLabelSize(0.15)
                sigTot[icat].GetXaxis().SetTitle('')
                # FIXME need to make the font larger... how?
                #if dosoverb:
                #    sigTot[icat].GetYaxis().SetTitle('S/B')
                #if dosoversqrtb:
                #    sigTot[icat].GetYaxis().SetTitle('S/#sqrt{B}')
                sigTot[icat].Draw('PE')

            if dolegend:
                if dodivide or dosoverb or dosoversqrtb:
                    if docats:
                        iPad = (icat%Ncol)+1+(icat/Ncol)*Ncol * 2
                        can.cd(iPad)
                    else:
                        can.cd(1)
                    #pad_id = icat%Ncol + 1 + (((icat)/Ncol) * 2 ) * Ncol # (icat%Ncol)+1 + ((icat/Ncol)+1)*Ncol
                    #can.cd(pad_id)
                legend.Draw()
                if nLeg==2:
                    legend2.Draw()
 
            if dotext:
                if dodivide or dosoverb or dosoversqrtb:
                    if docats:
                        iPad = icat%Ncol + 1 + ((icat)/Ncol) * 2 * Ncol
                        can.cd(iPad) 
                    else:
                        can.cd(1)
                plottext.Draw()
                if dotitles:
                    plottexttitle[icat] = SetTextTitle() # have to refresh to get rid of old category title
                    plottexttitle[icat].AddText(cattitles[icat])
                    plottexttitle[icat].Draw()
        
            if doline:
                if dodivide or dosoverb or dosoversqrtb:
                    if docats:
                        iPad = (icat%Ncol)+1+(icat/Ncol)*Ncol
                        can.cd(iPad)
                    else:
                        can.cd(1)
                cutline.SetLineWidth(4*plotscale)
                cutline.SetLineColor(633)
                cutline.DrawLine(linex,0.0,linex,stackmax/1.1);
 
            if docats: 
                if dodivide or dosoverb or dosoversqrtb:
                    pad_id = icat + 1 + (icat/Ncol)*Nrow  # (icat%Ncol)+1 + ((icat/Ncol)+1)*Ncol
                else:
                    pad_id = icat + 1
                can.cd(pad_id).Modified()
                can.cd(pad_id).SetTicks(1,1)
                can.cd(pad_id).RedrawAxis()
                can.cd(pad_id).Update()
            else:
                can.cd().Modified()
                can.cd().SetTicks(1,1)
                can.cd().RedrawAxis()
                can.cd().Update()
            
    if dointegrals:
        stackkeys = stackintegrals.keys()
        stackkeys.sort()
        for stackkey in stackkeys:
            print stackkey,'%.3e' %stackintegrals[stackkey]

    can.cd()    
    can.SetTicks(1,1)
    can.RedrawAxis()
    can.Update()
    if printsuffix is not "":
        if printcat is not -1:
            can.cd(icat).Print(str(cur_plot.plotvarname)+"_cat"+str(printcat)+"."+printsuffix)
        else:
            can.Print(str(cur_plot.plotvarname)+"."+printsuffix)

def AutoGroupNames():
    names=[]
    for sample in samples:
        name="group"+str(samples[sample].group)
        if name not in names:
            names.append(name)
    names=names.sort()
    print names
    return names

def SampleMap(sampletype):
    sampleIndexs={}  # key is order, maps to sampleIndex
    for sample in samples:
        if Debug:
            print "SampleMap, sampleIndex:",samples[sample].sampleIndex
        if samples[sample].plotsample==0:
            continue
        if sampletype == "sig" and sample < datasampleIndex:
            sampleIndexs[samples[sample].order]=samples[sample].sampleIndex
        elif sampletype == "data" and sample == datasampleIndex:
            sampleIndexs[samples[sample].order]=samples[sample].sampleIndex
        elif sampletype == "bkg" and sample >datasampleIndex:
            sampleIndexs[samples[sample].order]=samples[sample].sampleIndex
        elif sampletype == "all":
            sampleIndexs[samples[sample].order]=samples[sample].sampleIndex

    return sampleIndexs


def MakeMergeStack(sampletype):
        
    theStack = ROOT.THStack(sampletype,sampletype)
    if DebugNew:
        print "MakeMergeStack(sampletype)",sampletype
    
    sampleIndexs=SampleMap(sampletype)
    order = sampleIndexs.keys()
    order.sort()

    for index in order:
        sampleIndex = sampleIndexs[index]
        stacktop=False
        if DebugNew:
            print "SampleMap  index", index
            print "SampleMap samples[sampleIndex].plotsample",samples[sampleIndex].plotsample
        if samples[sampleIndex].plotsample == 1:
            for icat in xrange(cur_plot.ncat):
                if icat == xrange(cur_plot.ncat)[-1] and index == order[-1]:
                    if DebugNew:
                        print "icat,cur_plot.ncat,index,order[-1]",icat,cur_plot.ncat,index,order[-1]
                        print "icat == xrange(cur_plot.ncat)[-1]",icat == xrange(cur_plot.ncat)[-1]
                    stacktop=True
                histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[sampleIndex].inshortnames)
                if DebugNew:
                    print "SampleMap histname",histname
                rootfile.cd()
                hist=ROOT.gROOT.FindObject(histname).Clone()
                hist=FormatHistFull(hist,sampleIndex,sampletype,stacktop)
                theStack.Add(hist)
        
    return theStack
    


def MakeStack(sampletype, icat):

    sampleIndexs=SampleMap(sampletype)
    order = sampleIndexs.keys()
    order.sort()

    theStack = ROOT.THStack(sampletype+str(icat),sampletype+str(icat))
    for index in order:
        sampleIndex = sampleIndexs[index]
        if Debug:
            print "SampleMap  index", index
            print "SampleMap samples[sampleIndex].plotsample",samples[sampleIndex].plotsample
        if samples[sampleIndex].plotsample == 1:
            histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[sampleIndex].inshortnames)
            if Debug:
                print "SampleMap histname",histname
            rootfile.cd()
            hist=ROOT.gROOT.FindObject(histname).Clone()
            hist=FormatHist(hist,sampleIndex,sampletype)
            theStack.Add(hist)
        
    return theStack


def GetSampleIntegrals():
   
    sampleIndexs=SampleMap("all")
    order = sampleIndexs.keys()
    order.sort()
    
    integrals={}    

    if Debug:
        print "order",order
    for index in order:
        sampleIndex = sampleIndexs[index]
        integrals[sampleIndex]={}
        if Debug:
            print "samples[sampleIndex].plotsample",samples[sampleIndex].plotsample
            print "samples[sampleIndex].addtoleg",samples[sampleIndex].addtoleg
        if samples[sampleIndex].plotsample == 1:
            for icat in xrange(cur_plot.ncat):
                if Debug:
                    print "icat",icat
                histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[sampleIndex].inshortnames)
                if Debug:
                    print "histname",histname
                thishist=ROOT.gROOT.FindObject(histname)
                oflowbin = int(thishist.GetNbinsX()+1)
                if Debug:
                    print "oflowbin",oflowbin
                integrals[sampleIndex][icat]=thishist.Integral(0,oflowbin)

    for index in order:
        sampleIndex = sampleIndexs[index]
        summed=0
        for icat in xrange(cur_plot.ncat):
            if domergecats:
                summed=summed+integrals[sampleIndex][icat]
            else:
                print samples[sampleIndex].inshortnames,icat,integrals[sampleIndex][icat]
        if domergecats:
            print samples[sampleIndex].inshortnames,summed
                
        

def MakeMergeOverlayLines(sampletype):

    sampleIndexs=SampleMap(sampletype)
    if DebugNew:
        print "MakeOverlayLines", sampleIndexs
    order = sampleIndexs.keys()
    order.sort()
    if Debug:
        print "order",order
   
    # the key is a dictionary from order to sampleIndex
    # need sampleIndex info for plot properties 
    histmap={}

    first=1
    for index in order:
        if DebugNew:
            print "index",index
        sampleIndex = sampleIndexs[index]
        if DebugNew:
            print "samples[sampleIndex].plotsample",samples[sampleIndex].plotsample
            print "samples[sampleIndex].addtoleg",samples[sampleIndex].addtoleg
        if samples[sampleIndex].plotsample == 1:
            for icat in xrange(cur_plot.ncat):
                histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[sampleIndex].inshortnames)
                if DebugNew:
                    print "histname",histname
                rootfile.cd()
                if first == 1:
                    first=0
                    thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_index"+str(index))
                    if allintegrals:
                        print histname,thishist.Integral()
                    thishist=FormatHist(thishist,sampleIndex,sampletype)
                    fullhist=thishist.Clone("fullhist")
                else:
                    thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_index"+str(index))
                    thishist=FormatHist(thishist,sampleIndex,sampletype)
                    if allintegrals:
                        print histname,thishist.Integral()
                    thishist.Add(fullhist)
                    fullhist.Add(FormatHist(ROOT.gROOT.FindObject(histname).Clone(),sampleIndex,sampletype))
                        
                if  samples[sampleIndex].addtoleg is 1:
                    if DebugNew:
                        print "thishist.GetMaximum()",thishist.GetMaximum()
                        print "fullhist.GetMaximum()",fullhist.GetMaximum()
                        print "index,sampleIndex",index,sampleIndex
                    # individual histogram properties
                    #thishist.SetFillColor(samples[sampleIndex].color)
                    key=index,sampleIndex
                    histmap[key]=thishist

    if DebugNew:
        print "len(histmap)",len(histmap)
    return histmap


def MakeOverlayLines(sampletype, icat):

    sampleIndexs=SampleMap(sampletype)
    if DebugNew:
        print "MakeOverlayLines", sampleIndexs
    order = sampleIndexs.keys()
    order.sort()
    if Debug:
        print "order",order
   
    # the key is a dictionary from order to sampleIndex
    # need sampleIndex info for plot properties 
    histmap={}

    first=1
    for index in order:
        if DebugNew:
            print "index",index
        sampleIndex = sampleIndexs[index]
        if DebugNew:
            print "samples[sampleIndex].plotsample",samples[sampleIndex].plotsample
            print "samples[sampleIndex].addtoleg",samples[sampleIndex].addtoleg
        if samples[sampleIndex].plotsample == 1:
            histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[sampleIndex].inshortnames)
            if DebugNew:
                print "histname",histname
            rootfile.cd()
            if first == 1:
                first=0
                thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_cat"+str(icat)+"_index"+str(index))
                thishist=FormatHist(thishist,sampleIndex,sampletype)
                fullhist=thishist.Clone("fullhist")
            else:
                thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_cat"+str(icat)+"_index"+str(index))
                thishist=FormatHist(thishist,sampleIndex,sampletype)
                thishist.Add(fullhist)
                fullhist.Add(FormatHist(ROOT.gROOT.FindObject(histname).Clone(),sampleIndex,sampletype))
                    
            if  samples[sampleIndex].addtoleg is 1:
                if DebugNew:
                    print "thishist.GetMaximum()",thishist.GetMaximum()
                    print "fullhist.GetMaximum()",fullhist.GetMaximum()
                    print "index,sampleIndex",index,sampleIndex
                # individual histogram properties
                #thishist.SetFillColor(samples[sampleIndex].color)
                key=index,sampleIndex
                histmap[key]=thishist

    if DebugNew:
        print "len(histmap)",len(histmap)
    return histmap
            

def FormatHist(hist,sampleIndex,sampletype):
    return FormatHistFull(hist,sampleIndex,sampletype,True)


def FormatHistFull(hist,sampleIndex,sampletype,stacktop):
  
    if DebugNew:
        print "stacktop",stacktop
 
    if doscale: 
        hist.Scale(samples[sampleIndex].scale)
    
    if douflow:
        uflow = hist.GetBinContent(0)
        bin1content = hist.GetBinContent(1)
        hist.SetBinContent(1,bin1content+uflow)
    
    if dooflow:
        binN = hist.GetNbinsX()
        oflow = hist.GetBinContent(binN+1)
        binNcontent = hist.GetBinContent(binN)
        hist.SetBinContent(binN,binNcontent+oflow)
    
    if stacktop == True:
        if sampletype=="sig":
            hist.SetLineColor(int(samples[sampleIndex].color))
            hist.SetLineWidth(linewidth*plotscale)

        if sampletype=="data":
            hist.SetFillColor(ROOT.kBlack)
            hist.SetFillStyle(0)
            hist.SetMarkerStyle(20)
            if docats:
                hist.SetMarkerSize(0.6*plotscale)
            else:
                hist.SetMarkerSize(1.1*plotscale)
    else:
        hist.SetLineWidth(0)   
        hist.SetFillStyle(0)
        hist.SetMarkerSize(0)
 
    if dorebin:
        hist.Rebin(NReBin)

    return hist


def SetLegend():
    global legend

    if nLeg==1:
        legend = ROOT.TLegend(legx1,legy1,legx2,legy2)
    elif nLeg==2:
        legend = ROOT.TLegend(legx1,legy1,(legx2+legx1)/2.,legy2)
        global legend2
        legend2 = ROOT.TLegend((legx1+legx2)/2,legy1,legx2,legy2)
   
    legend.SetBorderSize(0)
    legend.SetTextFont(62)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(linewidth*plotscale)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
    if nLeg==2:
        legend2.SetBorderSize(0)
        legend2.SetTextFont(62)
        legend2.SetLineColor(0)
        legend2.SetLineStyle(1)
        legend2.SetLineWidth(linewidth*plotscale)
        legend2.SetFillColor(ROOT.kWhite)
        legend2.SetFillStyle(0)


def SetLine(newx):
    global linex,doline
    if newx==-1000:
        doline=False
    else:
        doline=True
    linex=newx
    ReDraw()


def MoveLegX(dx):
    global legx1,legx2
    legx1=legx1+dx
    legx2=legx2+dx
    ReDraw()


def MoveLegY(dy):
    global legy1,legy2
    legy1=legy1+dy
    legy2=legy2+dy
    ReDraw()


def SetText():
    global plottext
   
    plottext = ROOT.TPaveText(textx1,texty1,textx2,texty2,"brNDC");
    plottext.SetTextFont(62);
    plottext.SetTextSize(0.0425)
    plottext.SetBorderSize(0)
    plottext.SetLineColor(0)
    plottext.SetLineStyle(0)
    plottext.SetLineWidth(0)
    plottext.SetFillColor(0)
    plottext.SetFillStyle(0)

    if Debug:
        print text
        print text.find("\\n")
    if text.find("\\n") != -1:
        if Debug:
            print "text",text
        for line in text.split("\\n"):
            if Debug:
                print "line",line
            plottext.AddText(str(line))
    else:
        plottext.AddText(text)
    
def SetTextTitle():
    plottexttitle = ROOT.TPaveText(textx1,texty1-0.1,textx2,texty1,"brNDC");
    plottexttitle.SetTextFont(62);
    plottexttitle.SetTextSize(0.0425)
    plottexttitle.SetBorderSize(0)
    plottexttitle.SetLineColor(0)
    plottexttitle.SetLineStyle(0)
    plottexttitle.SetLineWidth(0)
    plottexttitle.SetFillColor(0)
    plottexttitle.SetFillStyle(0)
    return plottexttitle

def MoveTextX(dx):
    global textx1,textx2
    textx1=textx1+dx
    textx2=textx2+dx
    ReDraw()


def MoveTextY(dy):
    global texty1,texty2
    texty1=texty1+dy
    texty2=texty2+dy
    ReDraw()

def SetCat(cat):
    global singlecat, docats
    if int(cur_plot.ncat) < int(cat)+1:
        print "Cat",cat,"is larger than cur_plot.ncat",cur_plot.ncat
        print "Setting cat to ncat-1"
        cat=cur_plot.ncat-1
    if cat == -1:
        docats=True
    else:
        docats=False
    singlecat=cat
    ReDraw()


# Plot navigation
    
def PrevPlot():
    if globals().has_key("cur_plot"):
        if DebugNew:
            print "cur_plot.index",cur_plot.index
            print "len(plotinfos.keys())",len(plotinfos.keys())
        if cur_plot.index==0:
            prevnum=len(plotinfos.keys())-1
        else:
            prevnum=cur_plot.index-1
    else:
        prevnum=len(plotinfos.keys())-1
    Plot(prevnum)

def NextPlot():
    if globals().has_key("cur_plot"):
        if DebugNew:
            print "cur_plot.index",cur_plot.index
            print "len(plotinfos.keys())",len(plotinfos.keys())
        if cur_plot.index==len(plotinfos.keys())-1:
            nextnum=0
        else:
            nextnum=cur_plot.index+1
    else:
        nextnum=0
    Plot(nextnum)

def ReDraw():
    global plotscale
    plotscale=max(1, ysize/700)
    if doreplot==True:
        num=cur_plot.index
        Plot(num)

def Write(suffix):
    num=cur_plot.index
    Plot(num,suffix)

def WriteCat(suffix, icat):
    num=cur_plot.index
    Plot(num,suffix,icat)


option="START"
ENDLIST=[".q","q","quit","end","exit"]

while option not in ENDLIST:
    FindFunction(option)
    #try:
    #    FindFunction(option)
    #except:
    #    print "Option",option,"failed to run."
    option=str(raw_input("Enter option:  "))


