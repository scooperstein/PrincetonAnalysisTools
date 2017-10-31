import ROOT
import sys

ifile = ROOT.TFile(sys.argv[1],"r")
ofile = ROOT.TFile(sys.argv[2],"RECREATE")
#ofile = ROOT.TFile(sys.argv[1],"UPDATE")

cats = []
variables = []

for key in ifile.GetListOfKeys():
     name = key.GetName()
     if (name == "inputfiles" or name == "plotvariables"):
         ofile.cd()
         obj = key.ReadObj()
         #obj.Write(name)
         continue
     print name
     sample = name[name.find("cat")+5:]
     cat = int(name[name.find("cat")+3:name.find("cat")+4])
     var = name[:name.find("cat")-1]
     print var,cat,sample
     obj = key.ReadObj()
     ofile.cd()
     #obj.Write(key.GetName())
     if cat not in cats:
         cats.append(cat)
     if var not in variables:
         variables.append(var)

print cats
print variables

#WJetsScale = 1.08
#WJetsScale = [0.9958*1.021, 0.9668*1.021]
#WJetsScale = [0.9807,0.9775,]
WJetsScale = [1.0,1.0]

hQCDs = []
for cat in cats:
    for var in variables:
         print cat,var
         hData = ifile.Get("%s_cat%i_Data" % (var,cat))
         hQCD = hData.Clone("%s_cat%i_QCD_HT100To200" % (var,cat))
         print hQCD.Integral()
         hQCD.Sumw2(True)
         for key in ifile.GetListOfKeys():
             name = key.GetName()
             if (name.find("QCD")!=-1 or name.find("Data") != -1) or name.find("DYToLL") != -1: continue
             curr_var = name[:name.find("cat")-1]
             if (curr_var != var): continue
             if (name.find("_cat%i" % cat)==-1): continue
             obj = key.ReadObj()
             if (name.find("WJets")==-1):
                 hQCD.Add(obj,-1)
             elif (name.find("WJets_madgraph")!=-1 or name.find("EWKWJets_herwig")!=-1): continue
             elif (name.find("WJets")!=-1 and name.find("WJets_0J")==-1 and name.find("WJets_1J")==-1 and name.find("WJets_2J")==-1): continue
             else:
                 hQCD.Add(obj,-1.0*WJetsScale[cat])
             #elif (name.find("WJets_0J")!=-1):
             #    #hQCD.Add(obj,-0.8*WJetsScale[cat])
             #    #hQCD.Add(obj,-0.86)
             #elif (name.find("WJets_1J")!=-1):
             #    #hQCD.Add(obj,-1.0*WJetsScale[cat])
             #    #hQCD.Add(obj,-1.08)
             #elif (name.find("WJets_2J")!=-1):
             #    #hQCD.Add(obj,-1.18)
             #    #hQCD.Add(obj,-1.09*WJetsScale[cat])
             #else:
             #    hQCD.Add(obj,-1)
               
                 
             print "after sub. ",name,obj.Integral()
             print hQCD.Integral()
         if (cat==0): hQCD.Scale(1.31683)
         if (cat==1): hQCD.Scale(3.28502)
         #if (cat==0): hQCD.Scale(1.67074)
         #if (cat==1): hQCD.Scale(3.98831)

         # don't let the QCD prediction be negative in any bins
         nbin = hQCD.GetNbinsX()
         for i in range(1,nbin+1):
             val = hQCD.GetBinContent(i)
             if val < 0.:
                 hQCD.SetBinContent(i,0.0)
         hQCDs.append(hQCD)
ofile.cd()
for hQCD in hQCDs:
    hQCD.Sumw2(True)
    hQCD.Write()

