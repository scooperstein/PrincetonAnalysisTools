import ROOT
import sys
import numpy
import array
import math

f_LHC = 11245.6
t_LS=math.pow(2,18)/f_LHC

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


weightThreshold=1e-5

def GetWeightedValues(list):
    count=0
    sumOfWeights=0
    sumOfWeights2=0
    weightedSum=0

    for value,weight in list:
        #print value,weight
        if weight<weightThreshold:
            continue
        count=count+1
        sumOfWeights=sumOfWeights+weight
        sumOfWeights2=sumOfWeights2+math.pow(weight,2)
        weightedSum=weightedSum+weight*value

    return count,sumOfWeights,sumOfWeights2,weightedSum


def GetMean(list):
    #print "list length",len(list)
    count,sumOfWeights,sumOfWeights2,weightedSum=GetWeightedValues(list)
    mean=GetMeanFromWeightedValues(sumOfWeights,weightedSum)
    return mean


def GetMeanFromWeightedValues(sumOfWeights,weightedSum):
    mean=0
    if sumOfWeights>0:
        mean=weightedSum/sumOfWeights
    return mean


def GetMeanAndMeanError(list):
    count,sumOfWeights,sumOfWeights2,weightedSum=GetWeightedValues(list)
    if sumOfWeights2==0:
        return -99,-99
    neff=math.pow(sumOfWeights,2)/sumOfWeights2
    mean=GetMeanFromWeightedValues(sumOfWeights,weightedSum)

    #print neff,count,sumOfWeights
    
    weightedSumDiffFromAve2=0
    for value,weight in list:
        if weight<weightThreshold:
            continue
        weightedSumDiffFromAve2=weightedSumDiffFromAve2+weight*math.pow(value-mean,2) 

    stddev=0
    meanError=0
    if count>2:
        stddev=math.sqrt( weightedSumDiffFromAve2 / (sumOfWeights))
        meanError=stddev/math.sqrt(neff)

    #print "stddev",stddev

    return mean,meanError


csvfilenames=sys.argv[2:-1]
print csvfilenames

onlineLumi={} #(run,LS,LN)

fields = ['Fill','Run','LS','NB4','Mode','secs','msecs','deadfrac','PrimarySource','PrimaryLumi','HF','HFRaw','PLT','PLTRaw','PLTZero','PLTZeroRaw','BCMF','BCMFRaw','Ncol','Text','I','L[I]']
for csvfilename in csvfilenames:
    print csvfilename
    lines = open(csvfilename,'r').readlines()
    for line in lines : 
        vals = line.split(',')
        try:
            key=(int(vals[fields.index('Run')]),int(vals[fields.index('LS')]),int(vals[fields.index('NB4')]))
            if onlineLumi.has_key(key):
                print "onlineLumi ALREADY HAS KEY",key
                print onlineLumi[key]['line']
                print line
            
            onlineLumi[key]={}
            onlineLumi[key]['line']=line
            for field in fields:
                val=vals[fields.index(field)]
                if is_number(val):
                    val=float(val)
                onlineLumi[key][field]=val
        except:
            print "didn't work for line"
            print line

print "len(onlineLumi),",len(onlineLumi)


onlineLumiPerLSList={} #(run,LS)
onlineLumiPerLSMerged={} #(run,LS)

for LN in onlineLumi.keys():
    if not onlineLumiPerLSList.has_key((LN[0],LN[1])):
        onlineLumiPerLSList[(LN[0],LN[1])]=[]
        onlineLumiPerLSMerged[(LN[0],LN[1])]={}
    onlineLumiPerLSList[(LN[0],LN[1])].append(onlineLumi[LN])

for LS_key in onlineLumiPerLSList:
    unMergedLists={}
    mergedLists={}
    for field in fields:
        unMergedLists[field]=[]
        mergedLists[field]=0
    for NB_dict in onlineLumiPerLSList[LS_key]:
        for field in fields:
            try:
                unMergedLists[field].append(NB_dict[field])
            except:
                print "unMergedList has no",field,"for",LS_key


    for field in unMergedLists:
        item_list=unMergedLists[field]
        mergedLists[field]=unMergedLists[field]
        if len(item_list)>0:
            if is_number(item_list[0]):
                mergedLists[field]=sum(unMergedLists[field])/len(unMergedLists[field])

    onlineLumiPerLSMerged[LS_key]=mergedLists
    

filename=sys.argv[1]
tfile=ROOT.TFile(filename)
tree=tfile.Get("lumi/tree")

tree.SetBranchStatus("*",0)
#tree.SetBranchStatus("nVtx*",1)
#tree.SetBranchStatus("orbit*",1)
tree.SetBranchStatus("nPU",1)
tree.SetBranchStatus("run*",1)
tree.SetBranchStatus("LS*",1)
tree.SetBranchStatus("event*",1)
tree.SetBranchStatus("nPixelClusters*",1)
tree.SetBranchStatus("layer*",1)
tree.SetBranchStatus("BXNo",1)

newfile=ROOT.TFile("output_tree.root","recreate")
newtree=ROOT.TTree("newtree","validationtree")

run = array.array( 'l', [ 0 ] )
LS  = array.array( 'l', [ 0 ] )
nBX = array.array( 'l', [ 0 ] )
nCluster    = array.array( 'd', [ 0 ] )
nPCPerLayer = array.array( 'd', 5*[ 0 ] )

HFLumi    = array.array( 'd', [ 0 ] )
BCMFLumi  = array.array( 'd', [ 0 ] )
PLTLumi   = array.array( 'd', [ 0 ] )
BestLumi  = array.array( 'd', [ 0 ] )

HFLumi_integrated    = array.array( 'd', [ 0 ] )
BCMFLumi_integrated  = array.array( 'd', [ 0 ] )
PLTLumi_integrated   = array.array( 'd', [ 0 ] )
BestLumi_integrated  = array.array( 'd', [ 0 ] )

hasBrilData = array.array('b', [0])
hasCMSData  = array.array('b', [0])

pixel_xsec         = array.array( 'd', [ 0 ] )
pixel_xsec_layers  = array.array( 'd', 5*[ 0 ] )

nBX[0]=2
nPCPerBXid  = array.array( 'd', 3600*[ 0 ] )
BXid        = array.array( 'd', 3600*[ 0 ] )

newtree.Branch("run",run,"run/I")
newtree.Branch("LS",LS,"LS/I")
newtree.Branch("nBX",nBX,"nBX/I")
newtree.Branch("nCluster",nCluster,"nCluster/D")
newtree.Branch("nPCPerLayer",nPCPerLayer,"nPCPerLayer[5]/D")

newtree.Branch("pixel_xsec",pixel_xsec,"pixel_xsec/D")
newtree.Branch("pixel_xsec_layers",pixel_xsec_layers,"pixel_xsec_layers[5]/D")

newtree.Branch("BXid",BXid,"BXid[nBX]/D")
newtree.Branch("nPCPerBXid",nPCPerBXid,"nPCPerBXid[nBX]/D")

newtree.Branch("BestLumi",BestLumi,"BestLumi/D")
newtree.Branch("HFLumi",HFLumi,"HFLumi/D")
newtree.Branch("BCMFLumi",BCMFLumi,"BCMFLumi/D")
newtree.Branch("PLTLumi",PLTLumi,"PLTLumi/D")

newtree.Branch("BestLumi_integrated",BestLumi_integrated,"BestLumi_integrated/D")
newtree.Branch("HFLumi_integrated",HFLumi_integrated,"HFLumi_integrated/D")
newtree.Branch("BCMFLumi_integrated",BCMFLumi_integrated,"BCMFLumi_integrated/D")
newtree.Branch("PLTLumi_integrated",PLTLumi_integrated,"PLTLumi_integrated/D")

newtree.Branch("hasBrilData",hasBrilData,"hasBrilData/O")
newtree.Branch("hasCMSData",hasCMSData,"hasCMSData/O")

maxBin=100
maxLS={}

#can=ROOT.TCanvas("can","",600,600)

    
nentries=tree.GetEntries()

pixelCounts={}
lumiEstimate={}
# key is bx,LS and LS

print nentries
maxNBX=0
for iev in range(nentries):
#for iev in range(60):
    tree.GetEntry(iev)
    #if tree.nPU < 70:
    #    continue
    #print "PU",tree.nPU
    if iev%1000==0:
        print "iev,nVtx",iev#,tree.nVtx
        print "(tree.run,tree.LS)",tree.LS
        print "len(tree.nPixelClusters)",len(tree.nPixelClusters)
        print "len(tree.layers)",len(tree.layers)
    if len(tree.nPixelClusters)==0:
        continue
    pixelCount={}
    bxids={}
    if pixelCount.has_key((tree.run,tree.LS)) == 0:
        pixelCount[(tree.run,tree.LS)]=0
        pixelCount[(tree.run,tree.LS)]=[0]*6
        pixelCount[(tree.run,tree.LS)].append({}) # for bx->counts
   
    if pixelCounts.has_key((tree.run,tree.LS)) == 0:
        #nLists=1+5+2 #total+nLayers+nBX
        pixelCounts[(tree.run,tree.LS)]=[[] for x in xrange(6)]
        pixelCounts[(tree.run,tree.LS)].append({})

    if not maxLS.has_key(tree.run):
        maxLS[tree.run]=0
    
    if tree.LS>maxLS[tree.run]:
        maxLS[tree.run]=tree.LS+5

    layerNumbers=[]
    for item in tree.layers:
        layerNumbers.append(item[1])

    #for layer in range(1,6):
    #    print "layer,number of modules,",layer,layerNumbers.count(layer)

    counter=0
    for item in tree.nPixelClusters:
        #if counter>20:
        #    break
        counter=counter+1
        bxid=item[0][0]
        module=item[0][1]
        layer=tree.layers[module]
        clusters=item[1]

        if layer==6:
            layer=1

        pixelCount[(tree.run,tree.LS)][layer]=pixelCount[(tree.run,tree.LS)][layer]+clusters
        if not pixelCount[(tree.run,tree.LS)][6].has_key(bxid):
            pixelCount[(tree.run,tree.LS)][6][bxid]=0

        if layer!=1:
            pixelCount[(tree.run,tree.LS)][6][bxid]=pixelCount[(tree.run,tree.LS)][6][bxid]+clusters
            pixelCount[(tree.run,tree.LS)][0]=pixelCount[(tree.run,tree.LS)][0]+clusters

    
        if bxids.has_key(bxid)==0:
           bxids[bxid]=1
        else:
           bxids[bxid]=bxids[bxid]+1

    pixelCounts[(tree.run,tree.LS)][0].append([pixelCount[(tree.run,tree.LS)][0]/float(tree.eventCounter),1])
    for layer in range(1,6):
        #print pixelCount[(tree.run,tree.LS,layer)], tree.eventCounter
        pixelCounts[(tree.run,tree.LS)][layer].append([pixelCount[(tree.run,tree.LS)][layer]/float(tree.eventCounter),1])
    for bxid in bxids:
        if not pixelCounts[(tree.run,tree.LS)][6].has_key(bxid):
            pixelCounts[(tree.run,tree.LS)][6][bxid]=[]
        pixelCounts[(tree.run,tree.LS)][6][bxid].append([pixelCount[(tree.run,tree.LS)][6][bxid]/float(tree.BXNo[bxid]),1])



cmskeys=pixelCounts.keys()
brilkeys=onlineLumiPerLSMerged.keys()
keys=list(set(cmskeys+brilkeys))

keys.sort()

hists={}
PCCPerLayer=[118.,44.3,39.2,34.9,22.3,23.9]
for key in keys:
    run[0]=key[0]
    LS[0]=key[1]
    #print key
    hasBrilData[0]=False
    hasCMSData[0]=False
    
    HFLumi[0]=-1
    BestLumi[0]=-1
    PLTLumi[0] =-1
    BCMFLumi[0]=-1
    
    HFLumi_integrated[0]=-1
    BestLumi_integrated[0]=-1
    PLTLumi_integrated[0] =-1
    BCMFLumi_integrated[0]=-1
        
    pixel_xsec[0]=-1

    try:
        for layer in range(0,5):
            pixel_xsec_layers[layer]=-1
            nPCPerLayer[layer]=-1

        if key in brilkeys:
            hasBrilData[0]=True
            HFLumi[0]=onlineLumiPerLSMerged[key]['HF']
            BestLumi[0]=onlineLumiPerLSMerged[key]['PrimaryLumi']
            PLTLumi[0] =onlineLumiPerLSMerged[key]['PLT']
            BCMFLumi[0]=onlineLumiPerLSMerged[key]['BCMF']
            
            HFLumi_integrated[0]=onlineLumiPerLSMerged[key]['HF']*t_LS
            BestLumi_integrated[0]=onlineLumiPerLSMerged[key]['PrimaryLumi']*t_LS
            PLTLumi_integrated[0] =onlineLumiPerLSMerged[key]['PLT']*t_LS
            BCMFLumi_integrated[0]=onlineLumiPerLSMerged[key]['BCMF']*t_LS
            

        if key in cmskeys:
            hasCMSData[0]=True
            count=0
            for PCCs in pixelCounts[key]:
                if count==0:
                    mean,error=GetMeanAndMeanError(PCCs)
                    nCluster[0]=mean
                elif count<6:
                    mean,error=GetMeanAndMeanError(PCCs)
                    nPCPerLayer[count-1]=mean
                else:
                    ibx=0
                    nBX[0]=len(PCCs)
                    for bxid in PCCs:
                        mean,error=GetMeanAndMeanError(PCCs[bxid])
                        BXid[ibx]=bxid
                        print ibx,bxid,BXid[ibx]
                        nPCPerBXid[ibx]=mean
                        ibx=ibx+1
                        if ibx>nBX[0]:
                            print "ibx,nBX[0],",ibx,nBX[0],", but WHY?!!!"

                count=count+1 
  
        if hasCMSData[0] and hasBrilData[0]: 
            pixel_xsec[0]=nCluster[0]/BestLumi_integrated[0]*math.pow(2,18)
            for layer in range(0,5):
                pixel_xsec_layers[layer]=nPCPerLayer[layer]/BestLumi_integrated[0]*math.pow(2,18)


        newtree.Fill()

    except:
        print "I've failed me for the last time",key

newfile.Write()
newfile.Close()

sys.exit(0)




#
#if len(key) ==4:
#        NBX=1
#    elif len(key) == 3:
#        if key[0]>230000 and key[0]<400000:
#            layer=key[2]
#            NBX=2
#            run=key[0]
#            ls=key[1]
#            histkey=str(key[0])+"_layer"+str(key[2])
#        else:
#            NBX=1
#            run=key[1]
#            ls=key[2]
#            histkey=str(key[1])+"_bx"+str(key[0])
#    elif len(key) == 2:
#        histkey=str(key[0])
#        run=key[0]
#        ls=key[1]
#        NBX=2
#    else:
#        print "len(key)",len(key),"?"
#
#    if hists.has_key(histkey)==0:
#        title="PCC/LS in Run  "+str(run)
#        if histkey.find("bx") !=-1:
#            title=title+" (BX="+histkey.split("bx")[1]+")"
#        if histkey.find("layer") !=-1:
#            title=title+" (layer="+histkey.split("layer")[1]+")"
#        title=title+";LumiSection;"
#        titlePU=title+"AveragePU/LS"
#        title=title+"Average PCC/LS"
#        hists[histkey]=ROOT.TH1F(histkey,title,maxLS[run],0,maxLS[run])
#        hists[histkey+"PU"]=ROOT.TH1F(histkey+"PU",titlePU,maxLS[run],0,maxLS[run])
#
#    mean,error=GetMeanAndMeanError(pixelCounts[key])
#    hists[histkey].Fill(ls,mean)
#
#    meanPU=mean/float(PCCPerLayer[layer])
#    print run,ls,layer,mean,meanPU
#    hists[histkey+"PU"].Fill(ls,meanPU)
#    
#    lumi=meanPU * math.pow(2,18) * NBX / 78400


sys.exit(0)
tcan=ROOT.TCanvas("tcan","",900,700)

for histkey in hists:
    hists[histkey].Draw("hist")
    tcan.Update()
    tcan.SaveAs("hists/"+histkey+".png")
    #raw_input()
    


