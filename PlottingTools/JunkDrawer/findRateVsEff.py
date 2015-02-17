#! /usr/bin/env python

#   
#   This program takes a flat ROOT tree, 
#   applies selection on a few variables,
#   and writes a new tree that contains 
#   the efficiencies and rates at various
#   thresholds, which are also saved.
#   
#   Development is encouraged for more 
#   general usage.
#   
#   Written by Chris Palmer (Princeton)
#


import ROOT
import getopt, sys, array
import numpy
  

if __name__ == "__main__":  

#    if (rootInputFile == ""):
#        print "You have to select one input file."
#        usage()
#        sys.exit(3)

    file = ROOT.TFile("../VHbbTrigger/output_tightersig_phicut_v2.root")
    tree = file.Get("tree")

    optset=["elPt","jetPt1","jetPt2"]
    
    cuts={}
    #cuts["elPt"]=numpy.arange(10,46,2.5)
    #cuts["jetPt1"]=numpy.arange(10,121,10)
    #cuts["jetPt2"]=numpy.arange(10,101,10)
    cuts["elPt"]=numpy.arange(10,46,2.5)
    cuts["jetPt1"]=numpy.arange(10,121,10)
    cuts["jetPt2"]=numpy.arange(10,101,10)
    cuts["jetPt2"][0]=-10


    basesel="((notPassBaseline==1&&sampleIndex<0)||(passCurrentL1==0&&sampleIndex>0))"
    
    fullcuts={}

    finalHistos=[]
    for elpt in cuts["elPt"]:
        for jetpt1 in cuts["jetPt1"]:
            for jetpt2 in cuts["jetPt2"]:
                #fullcuts.["elIsoPt","jetPt1","jetPt2"]
                fullcuts[elpt,jetpt1,jetpt2]=[-1.,-1.,-1.,-1.]

    #types=["sampleIndex<0"]
    types=["sampleIndex<0","(sampleIndex<0&&genHPt>100)","(sampleIndex<0&&genHPt>200)","sampleIndex>0"]

    print "Looping over",str(len(fullcuts)*len(types))
    counter=0
    for itype in xrange(len(types)):
        for fullcut in fullcuts:
            if int(counter)%40 ==0:
                print "\t\t",str(counter),"of",str(len(fullcuts)*len(types))
            print fullcut
            fullstr="("+basesel+"&&elPt>"+str(fullcut[0])+"&&(hltElePt>"+str(fullcut[0])+"||sampleIndex>0)&&jetPt1>"+str(fullcut[1])+"&&jetPt2>"+str(fullcut[2])+"&&"+types[itype]+")*weight"
            temp_hist = "temp_hist(70,0,70)"
           
            print fullstr 
            tree.Draw("elPt>> "+temp_hist, fullstr, "goff")
            final_h = ROOT.gDirectory.Get("temp_hist")
            final_name = "elP_"+str(fullcut[0])+"_"+str(fullcut[1])+"_"+str(fullcut[2])
            final_h.SetName(final_name)
            final_h.SetTitle(final_name)
            finalHistos.append(final_h)

            fullcuts[fullcut][itype]=final_h.Integral()
            #print "integral, fullcut, itype,",str(fullcuts[fullcut][itype]),str(fullcut),str(itype)
            counter=counter+1 

    output = ROOT.TFile("testing_phicut_new.root", "recreate")
    newtree = ROOT.TTree("newtree","newtree")

    # FIXME make local variables with name of branch as key
    localVar={}
   
    for var in optset:
        localVar[var]=numpy.zeros(1, dtype=float) 

    localVar["eff"]=numpy.zeros(1, dtype=float)
    localVar["effhpt100"]=numpy.zeros(1, dtype=float)
    localVar["effhpt200"]=numpy.zeros(1, dtype=float)
    localVar["rate"]=numpy.zeros(1, dtype=float)

    for var in localVar:
        newtree.Branch(var,localVar[var], var+"/D")
    
    for cut in fullcuts:
        # FIXME access by key
        localVar["elPt"][0]=cut[0]   
        localVar["jetPt1"][0]=cut[1]   
        localVar["jetPt2"][0]=cut[2]   

        localVar["eff"][0]=fullcuts[cut][0]
        localVar["effhpt100"][0]=fullcuts[cut][1]
        localVar["effhpt200"][0]=fullcuts[cut][2]
        localVar["rate"][0]=fullcuts[cut][3]
        newtree.Fill()

    newtree.Write()
    output.Close()
