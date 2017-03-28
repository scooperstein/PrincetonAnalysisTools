# Copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TopQuarkAnalysis/SingleTop/src/TopProducer.cc?revision=1.9&view=markup
# and converted from C++ to Python
from ROOT import *
from math import pow
from math import sqrt
from math import acos
from math import cos
from math import sin

##
# \Function EquationSolver:
#
# Solves 3rd degree equations
#
# \Author A. Orso M. Iorio
# 
#
# \version  $Id: EquationSolver.h,v 1.1 2013/02/27 12:18:42 degrutto Exp $
#
def EquationSolve(a, b, c, d):


  result = []



  if (a != 0):
    
    q = (3*a*c-b*b)/(9*a*a)
    r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a)
    Delta = q*q*q + r*r

    rho=0.
    theta=0.
    
    if( Delta<=0):
      rho = sqrt(-(q*q*q))

      theta = acos(r/rho)
  
      s = complex(sqrt(-q)*cos(theta/3.0),sqrt(-q)*sin(theta/3.0))
      t = complex(sqrt(-q)*cos(-theta/3.0),sqrt(-q)*sin(-theta/3.0))
    
    if(Delta>0):
      #print r, sqrt(Delta)
      if (r+sqrt(Delta) > 0): 
          s = complex(pow((r+sqrt(Delta)),(1./3)),0)
      else:
          s = complex(-pow(abs((r+sqrt(Delta))),(1./3)),0)
      if (r-sqrt(Delta) > 0):
          t = complex(pow((r-sqrt(Delta)),(1./3)),0)
      else:
          t = complex(-pow(abs((r-sqrt(Delta))),(1./3)),0)
    
    i = complex(0.,1.0) 
    
    x1 = s+t+complex(-b/(3.0*a),0)
    x2 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)+(s-t)*i*complex(sqrt(3)/2.0,0)
    x3 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)-(s-t)*i*complex(sqrt(3)/2.0,0)

    if(abs(x1.imag)<0.0001): result.append(x1.real)
    if(abs(x2.imag)<0.0001): result.append(x2.real)
    if(abs(x3.imag)<0.0001): result.append(x3.real)

    #print x1,x2,x3
    return result
  else:
      return result


  return result

def getNu4Momentum(TLepton, TMET):

  Lepton = TLorentzVector()
  Lepton.SetPxPyPzE(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
  MET = TLorentzVector()
  MET.SetPxPyPzE(TMET.Px(), TMET.Py(), 0., TMET.E());

  mW = 80.38;

  result = []

  MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py());
  mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py();
  a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
  a2 = pow(a,2);
  b  = -10*(pow(Lepton.Energy(),2.)*(MisET2) - pow(mu,2.))/(pow(Lepton.Energy(),2) - pow(Lepton.Pz(),2));
  pz1 = 0.
  pz2 = 0.
  pznu = 0.
  nNuSol = 0

  p4nu_rec = TLorentzVector()
  p4W_rec = TLorentzVector()
  p4b_rec = TLorentzVector()
  p4Top_rec = TLorentzVector()
  p4lep_rec = TLorentzVector()

  p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());
 
  #print a2,b 
  if(a2-b > 0 ):
    root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
    nNuSol = 2;

    pznu = pz1;

    Enu = sqrt(MisET2 + pznu*pznu);

    p4nu_rec.SetPxPyPzE(MET.Px(), MET.Py(), pznu, Enu);

    result.append(p4nu_rec);

  else:

    ptlep = Lepton.Pt()
    pxlep=Lepton.Px()
    pylep=Lepton.Py()
    metpx=MET.Px()
    metpy=MET.Py()

    EquationA = 1.
    EquationB = -3.*pylep*mW/(ptlep)
    EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep)
    EquationD = 4.*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep

    solutions = EquationSolve(EquationA,EquationB,EquationC,EquationD)

    solutions2 = EquationSolve(EquationA,EquationB,EquationC,EquationD);

    deltaMin = 14000*14000
    zeroValue = -mW*mW/(4*pxlep)
    minPx=0
    minPy=0

    for i in range(len(solutions)):
        if(solutions[i]<0 ): continue
        p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep)
        p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep)
        Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)

        if(Delta2< deltaMin and Delta2 > 0):
            deltaMin = Delta2
            minPx=p_x
            minPy=p_y

    for i in range(len(solutions2)):
        if(solutions2[i]<0 ): continue
        p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep)
        p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep)
        Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)
        if(Delta2< deltaMin and Delta2 > 0):
            deltaMin = Delta2
            minPx=p_x
            minPy=p_y

    pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue)
    delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy)

    if(deltaMin==14000*14000): return TLorentzVector(0,0,0,0)

    if(delta2ZeroValue < deltaMin):
        deltaMin = delta2ZeroValue
        minPx=zeroValue
        minPy=pyZeroValue

    mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep
    a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz())
    pznu = a_Minimum

    Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu)
    p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu)
    result.append(p4nu_rec)
  return result[0]

def computeTopMass(lep, met, jets):
    neutrino = getNu4Momentum(lep, met)
    bjet = TLorentzVector()
    minDR = 99
    for jet in jets:
        if (jet.Pt() > 30): # and jet.bTagCSV > CSVL and jet.puID > 0 and jet.Id > 0? how?
            dR = jet.DeltaR(lep)
            if (dR < minDR):
                minDR = dR
                bjet = jet
    if (bjet.Pt() <= 0):
        return -99    
    top = lep + met + bjet
    return top.M()


# just some basic tests

import random

nIter = 1000000
for i in range(nIter):
    if (i % 10000 == 0):
        print "processing entry: %i" % i

    random.seed(i)
    rand1 = random.uniform(-2,2)
    random.seed(i+nIter)
    rand2 = random.uniform(-2,2)
    random.seed(i+2*nIter)
    rand3 = random.uniform(-2,2)
    random.seed(i+3*nIter)
    rand4 = random.uniform(-2,2)
    
    lep = TLorentzVector()
    met = TLorentzVector()
    lep.SetPtEtaPhiM(rand1*86.575485, rand1*-0.370986, rand1*-1.694283, rand1*0.1057000)
    met.SetPtEtaPhiM(rand2*39.530349, 0., rand2*-2.810159, 0.0)
    #neutrino = getNu4Momentum(lep,met)
    #print neutrino.Px(),neutrino.Py(),neutrino.Pz()

    bjet1 = TLorentzVector()
    bjet2 = TLorentzVector()
    bjet1.SetPtEtaPhiM(rand3*156.64633, rand3*0.6337512, rand3*-1.150040, rand3*19.657625)
    bjet2.SetPtEtaPhiM(rand4*257.87811, rand4*-1.495820, rand4*.5050827,  rand4*39.884072)
    #hbb = bjet1 + bjet2
    #print hbb.M()
    #W = lep + neutrino
    #deta = abs(hbb.Eta() - W.Eta())
    #print deta
    #print W.M()
    jets = [bjet1, bjet2]
    tmp = computeTopMass(lep,met,jets)
    #print computeTopMass(lep,met,jets)

