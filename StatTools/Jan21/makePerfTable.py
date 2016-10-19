import sys

ifile = open(sys.argv[1])
sigs = ""
mus = ""
rs = ""
for line in ifile:
    line = line.strip() 
    if (line.find("150_4")!=-1 or line.find("150_5")!=-1):
        sigs += "\n"
        mus += "\n"
        rs += "\n"
    if (line.find("Significance")!=-1):
      sigs += line.split()[1] + " "
    if (line.find("Best fit")!=-1):
      mus += line.split()[4].split('/')[1][1:] + " "
    if (line.find("Expected")!=-1):
      rs += line.split()[4] + " "
print "Significance"
print sigs
print "Mu Unc."
print mus
print "Limit"
print rs
