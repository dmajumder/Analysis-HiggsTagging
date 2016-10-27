#!/usr/bin/env python

import sys, os, shutil, re, subprocess

subprocess.call("eval scramv1 runtime -csh", shell=True)
subprocess.call("source /afs/cern.ch/user/d/devdatta/.tcshrc", shell=True)

import ROOT
#from rootpy.vector import LorentzVector
#from rootpy.tree import Tree, TreeModel, FloatArrayCol, IntCol
#from rootpy.io import root_open

def getfiles(input_dir):
 proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'ls', input_dir ], 
   stdout = subprocess.PIPE, 
   stderr = subprocess.STDOUT )
 output = proc.communicate()[0]
 if proc.returncode != 0:
   print output
   sys.exit(1)
 return output.splitlines()

def main():

  fout = ROOT.TFile("AK8Regressed.root", "RECREATE")
  factory = ROOT.TMVA.Factory("TMVARegression", fout,
      "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Regression")
  factory.AddVariable("pt_AK8MatchedToHbb", "pt_AK8MatchedToHbb", "GeV", "F")
  factory.AddVariable("eta_AK8MatchedToHbb", "eta_AK8MatchedToHbb", "", "F")
  factory.AddVariable("nsv_AK8MatchedToHbb", "nsv_AK8MatchedToHbb", "", "F")
  factory.AddVariable("sv0mass_AK8MatchedToHbb", "sv0mass_AK8MatchedToHbb", "", "F")
  factory.AddVariable("sv1mass_AK8MatchedToHbb", "sv1mass_AK8MatchedToHbb", "", "F")
  factory.AddVariable("nch_AK8MatchedToHbb", "nch_AK8MatchedToHbb", "", "F")
  factory.AddVariable("nmu_AK8MatchedToHbb", "nmu_AK8MatchedToHbb", "", "F")
  factory.AddVariable("nel_AK8MatchedToHbb", "nel_AK8MatchedToHbb", "", "F")
  factory.AddVariable("muenfr_AK8MatchedToHbb", "muenfr_AK8MatchedToHbb", "", "F")
  factory.AddVariable("emenfr_AK8MatchedToHbb", "emenfr_AK8MatchedToHbb", "", "F")

  factory.AddTarget("pt_MatchedHbb/pt_AK8MatchedToHbb")
  #factory.SetWeightExpression("n_pv")
  factory.AddSpectator("n_pv")
  factory.AddSpectator("msoftdrop_AK8MatchedToHbb")
  
  fin = ROOT.TFile.Open("TrainingTree_BulkGravTohhTohbbhbb_narrow.root")  
  traintree = fin.Get("TrainingTree")
  testtree = fin.Get("TrainingTree")
  factory.AddRegressionTree(traintree, 1.0, ROOT.TMVA.Types.kTraining)
  factory.AddRegressionTree(testtree, 1.0, ROOT.TMVA.Types.kTesting)

  cut = ROOT.TCut()

  factory.PrepareTrainingAndTestTree( cut, "V:nTrain_Regression=100000:nTest_Regression=100000:SplitMode=Random:SplitSeed=0:NormMode=NumEvents" );
  factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDTG",
                                   "!H:V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.7:nCuts=20:MaxDepth=3:NNodesMax=15" );
  factory.TrainAllMethods()
  factory.TestAllMethods()
  factory.EvaluateAllMethods()

  fout.Close()

  #for event in tree:
  #  print event.pt_Hbb
  
if __name__ == "__main__":
  main()
