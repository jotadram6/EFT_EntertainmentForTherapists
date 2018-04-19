import commands as cmd
import argparse
import ROOT
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--DIR', help='Directory to write output')
parser.add_argument('--XS', help='Sample cross section')
parser.add_argument('--MT', help='Minimum MT')
parser.add_argument('--OUT', help='ROOT output file')
args = parser.parse_args()

def DeltaPhi(phi1,phi2):
    PHI=phi1-phi2
    if PHI >= pi:
        PHI -= 2*pi
    elif PHI < -1*pi:
        PHI += 2*pi
    return PHI

def MT(Lpt,MET,Lphi,METphi):
    return sqrt(2*Lpt*MET*(1-cos(DeltaPhi(METphi,Lphi))))

#Delphes_Path="/home/joser/Pheno_Studies/MG5_aMC_v2_6_1/Delphes/"
Delphes_Path="/home/joser/Dropbox/Vandy/EFT/ElectronNeutrino/Delphes-3.4.1/"
ROOT.gSystem.AddDynamicPath(Delphes_Path)
#ROOT.gROOT.ProcessLine("#include <math.h>")
ROOT.gSystem.Load("libDelphes");

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
  print "Delphes classes imported"
except:
  pass

TreeName="Delphes"
DataChain=ROOT.TChain(TreeName)
DataChain.Add(args.DIR+"*.root")

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(DataChain)
numberOfEntries = treeReader.GetEntries()

branchElectron = treeReader.UseBranch("Electron")
branchMET = treeReader.UseBranch("MissingET")

#ListOfEventNumber=[]

NumberOfEventsToCheck=DataChain.GetEntries()

print "--------------------------------------> Going to analyze", NumberOfEventsToCheck, "events!!!!!!!!!!!"

CountingEvents=0

if args.OUT is not None:
    RootFile = ROOT.TFile(args.OUT,"recreate")
    MThisto = ROOT.TH1F("MT","MT",100,500.,3000.)

for entry in xrange(NumberOfEventsToCheck):
    treeReader.ReadEntry(entry)
    #Selection from EXO-15-006
    if branchElectron.GetEntries() != 1: continue
    if branchElectron.At(0).PT<=130: continue
    if branchElectron.At(0).SumPt>=5: continue
    if abs(DeltaPhi(branchElectron.At(0).Phi,branchMET.At(0).Phi))<=2.5: continue
    if (branchElectron.At(0).PT/branchMET.At(0).MET)>0.4 and (branchElectron.At(0).PT/branchMET.At(0).MET)<1.5:
        EventMT=MT(branchElectron.At(0).PT,branchMET.At(0).MET,branchElectron.At(0).Phi,branchMET.At(0).Phi)
        if EventMT>float(args.MT):
            CountingEvents=CountingEvents+1
            MThisto.Fill(EventMT)

"""for event in xrange(NumberOfEventsToCheck):
    DataChain.GetEntry(event)
    #print DataChain.eventNumber
    #Selection from EXO-15-006
    if DataChain.Electron_size!=1: continue
    if DataChain.Electron.PT<=130: continue
    if DataChain.Electron.SumPt>=5: continue
    if abs(DeltaPhi(DataChain.Electron.Phi,DataChain.MissingET.Phi))<=2.5: continue
    if (DataChain.Electron.PT/DataChain.MissingET.MET)>0.4 and (DataChain.Electron.PT/DataChain.MissingET.MET)<1.5:
        CountingEvents=CountingEvents+1
"""



print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=10150 #pb
TotalXS=float(args.XS)
Lumi=2200 #pb-1
Weight=(Lumi*TotalXS/NumberOfEventsToCheck)
print "Weighted number of events that passed full selection:", CountingEvents*Weight, "+-", sqrt(CountingEvents)*Weight
if args.OUT is not None:
    MThisto.Sumw2()
    MThisto.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(MThisto.GetNbinsX()):
        print str(MThisto.GetBinLowEdge(i+1))+"-"+str(MThisto.GetBinLowEdge(i+1)+MThisto.GetBinWidth(i+1))+",", MThisto.GetBinContent(i+1)
    MThisto.Write()
    RootFile.Close()
