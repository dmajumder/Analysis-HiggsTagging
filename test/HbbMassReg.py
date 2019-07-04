#!/usr/bin/env python

import sys, os, shutil, re, subprocess
import argparse
import numpy as np
from array import array

subprocess.call("eval scramv1 runtime -csh", shell=True)
subprocess.call("source /afs/cern.ch/user/d/devdatta/.tcshrc", shell=True)

import ROOT

def train():

  fout = ROOT.TFile("HbbRegressed_DiPho_GenJetsWithNu_No_reco_Mjj_and_pT_Cut_dR_p3_matching.root", "RECREATE")
  factory = ROOT.TMVA.Factory("TMVARegression", fout,
      "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Regression")

  dataloader = ROOT.TMVA.DataLoader("dataset")
  dataloader.AddVariable ("reg_recoJet_1_pt"    , "reg_recoJet_1_pt"    , "GeV", "F")
  dataloader.AddVariable ("reg_recoJet_1_eta"   , "reg_recoJet_1_eta"   , ""   , "F")
  dataloader.AddVariable ("reg_recoJet_1_phi"   , "reg_recoJet_1_phi"   , ""   , "F")
  dataloader.AddVariable ("reg_recoJet_1_mass"  , "reg_recoJet_1_mass"  , "GeV", "F")
  dataloader.AddVariable ("reg_recoJet_1_e"     , "reg_recoJet_1_e"     , "GeV", "F")
  dataloader.AddVariable ("reg_recoJet_2_pt"    , "reg_recoJet_2_pt"    , "GeV", "F")
  dataloader.AddVariable ("reg_recoJet_2_eta"   , "reg_recoJet_2_eta"   , ""   , "F")
  dataloader.AddVariable ("reg_recoJet_2_phi"   , "reg_recoJet_2_phi"   , ""   , "F")
  dataloader.AddVariable ("reg_recoJet_2_mass"  , "reg_recoJet_2_mass"  , "GeV", "F")
  dataloader.AddVariable ("reg_recoJet_2_e"     , "reg_recoJet_2_e"     , "GeV", "F")
  dataloader.AddVariable ("Met_CorPt"           , "Met_CorPt"           , "GeV", "F")
  dataloader.AddVariable ("Met_CorPhi"          , "Met_CorPhi"          , ""   , "F")

  dataloader.AddSpectator("reg_reco_mjj"        , "reg_reco_mjj"        , "GeV"     )
  dataloader.AddSpectator("mbb"                 , "mbb"                 , "GeV"     )
  dataloader.AddSpectator("mbbNu"               , "mbbNu"               , "GeV"     )

  dataloader.AddRegressionTarget("mbbNu/reg_reco_mjj")

  #fin = ROOT.TFile.Open("/afs/cern.ch/work/l/lata/public/ForDevdatta/Regression/2017/root_files_dR_p3_matching/node_SM_2017.root", "READ")
  fin = ROOT.TFile.Open("/afs/cern.ch/work/l/lata/public/ForDevdatta/Regression/2017/root_files_dR_p3_matching/DiPho_2017.root", "READ")
  traintree = fin.Get("myRegTree")
  testtree = fin.Get("myRegTree")
  dataloader.AddRegressionTree(traintree, 1.0, ROOT.TMVA.Types.kTraining)
  dataloader.AddRegressionTree(testtree, 1.0, ROOT.TMVA.Types.kTesting)

  cut = ROOT.TCut("dRmin_Jet1<0.2 && dRmin_Jet2 < 0.2")

  dataloader.PrepareTrainingAndTestTree( cut, "V:nTrain_Regression=10000:nTest_Regression=9000:SplitMode=Random:SplitSeed=0:NormMode=NumEvents" );
  factory.BookMethod( dataloader,  ROOT.TMVA.Types.kBDT, "BDTG",
      "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
  factory.TrainAllMethods()
  factory.TestAllMethods()
  factory.EvaluateAllMethods()

  fout.Write()
  fout.Close()




def test(verbose=False):

  reg_recoJet_1_pt   = np.empty(1, dtype='float32')  
  reg_recoJet_1_eta  = np.empty(1, dtype='float32')  
  reg_recoJet_1_phi  = np.empty(1, dtype='float32')  
  reg_recoJet_1_mass = np.empty(1, dtype='float32')  
  reg_recoJet_1_e    = np.empty(1, dtype='float32')  
  reg_recoJet_2_pt   = np.empty(1, dtype='float32')  
  reg_recoJet_2_eta  = np.empty(1, dtype='float32')  
  reg_recoJet_2_phi  = np.empty(1, dtype='float32')  
  reg_recoJet_2_mass = np.empty(1, dtype='float32')  
  reg_recoJet_2_e    = np.empty(1, dtype='float32')  
  Met_CorPt          = np.empty(1, dtype='float32')  
  Met_CorPhi         = np.empty(1, dtype='float32')  
  reg_reco_mjj       = np.empty(1, dtype='float32')  
  mbb                = np.empty(1, dtype='float32')  
  mbbNu              = np.empty(1, dtype='float32')  

  #fin = ROOT.TFile.Open("/afs/cern.ch/work/l/lata/public/ForDevdatta/Regression/2017/root_files_dR_p3_matching/DiPho_2017.root", "READ")
  fin = ROOT.TFile.Open("/afs/cern.ch/work/l/lata/public/ForDevdatta/Regression/2017/root_files_dR_p3_matching/node_SM_2017.root", "READ")
  tree = fin.Get('myRegTree')

  tree.SetBranchAddress("reg_recoJet_1_pt"    , reg_recoJet_1_pt   )
  tree.SetBranchAddress("reg_recoJet_1_eta"   , reg_recoJet_1_eta  )
  tree.SetBranchAddress("reg_recoJet_1_phi"   , reg_recoJet_1_phi  )
  tree.SetBranchAddress("reg_recoJet_1_mass"  , reg_recoJet_1_mass )
  tree.SetBranchAddress("reg_recoJet_1_e"     , reg_recoJet_1_e    )
  tree.SetBranchAddress("reg_recoJet_2_pt"    , reg_recoJet_2_pt   )
  tree.SetBranchAddress("reg_recoJet_2_eta"   , reg_recoJet_2_eta  )
  tree.SetBranchAddress("reg_recoJet_2_phi"   , reg_recoJet_2_phi  )
  tree.SetBranchAddress("reg_recoJet_2_mass"  , reg_recoJet_2_mass )
  tree.SetBranchAddress("reg_recoJet_2_e"     , reg_recoJet_2_e    )
  tree.SetBranchAddress("Met_CorPt"           , Met_CorPt          )
  tree.SetBranchAddress("Met_CorPhi"          , Met_CorPhi         )
  tree.SetBranchAddress("reg_reco_mjj"        , reg_reco_mjj       )
  tree.SetBranchAddress("mbb"                 , mbb                )
  tree.SetBranchAddress("mbbNu"               , mbbNu              )

  fout = ROOT.TFile('node_SM_reg_on_DiPho_No_reco_Mjj_and_pT_Cut_dR_p3_matching_jetEnergy.root', 'RECREATE')
  fout.cd()

  h_reco_mjj       = ROOT.TH1D('h_reco_mjj'    , ';reco m_{jj} [GeV];Events;;'          ,100, 0, 1000)
  h_reco_mjj_reg   = ROOT.TH1D('h_reco_mjj_reg', ';reco+regressed m_{jj} [GeV];Events;;',100, 0, 1000)
  h_regwt          = ROOT.TH1D('h_regwt'       , ';Regression weight;;'                 , 50, 0, 5)
  h_mjjres_reg     = ROOT.TH1D('h_mjjres_reg'  , ';(m_{jj}^{Reco+Reg.} - m_{jj})^{Gen}/m_{jj}^{Gen}', 50, -2, 2)
  h_mjjres         = ROOT.TH1D('h_mjjres'      , ';(m_{jj}^{Reco} - m_{jj})^{Gen}/m_{jj}^{Gen}', 50, -2, 2)

  methodname = "BDTG"
  #weightfile = "/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_9_4_13_patch2/src/Analysis/HiggsTagging/test/dataset/weights/TMVARegression_On_node_SM_2017_No_reco_Mjj_and_pT_Cut_dR_p3_matching_BDTG.weights.xml"
  weightfile = "/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_9_4_13_patch2/src/Analysis/HiggsTagging/test/dataset/weights/TMVARegression_On_DiPho_2017_No_reco_Mjj_and_pT_Cut_dR_p3_matching_jetEnergy_BDTG.weights.xml"

  ROOT.TMVA.Tools.Instance()
  reader = ROOT.TMVA.Reader( "!Color:Silent" );

  reader.AddVariable ("reg_recoJet_1_pt"    , reg_recoJet_1_pt   )
  reader.AddVariable ("reg_recoJet_1_eta"   , reg_recoJet_1_eta  )
  reader.AddVariable ("reg_recoJet_1_phi"   , reg_recoJet_1_phi  )
  reader.AddVariable ("reg_recoJet_1_mass"  , reg_recoJet_1_mass )
  reader.AddVariable ("reg_recoJet_1_e"     , reg_recoJet_1_e    )
  reader.AddVariable ("reg_recoJet_2_pt"    , reg_recoJet_2_pt   )
  reader.AddVariable ("reg_recoJet_2_eta"   , reg_recoJet_2_eta  )
  reader.AddVariable ("reg_recoJet_2_phi"   , reg_recoJet_2_phi  )
  reader.AddVariable ("reg_recoJet_2_mass"  , reg_recoJet_2_mass )
  reader.AddVariable ("reg_recoJet_2_e"     , reg_recoJet_2_e    )
  reader.AddVariable ("Met_CorPt"           , Met_CorPt          )
  reader.AddVariable ("Met_CorPhi"          , Met_CorPhi         )
  reader.AddSpectator("reg_reco_mjj"        , reg_reco_mjj       )
  reader.AddSpectator("mbb"                 , mbb                )
  reader.AddSpectator("mbbNu"               , mbbNu              )

  reader.BookMVA(methodname, weightfile)

  for ievt in range(0, tree.GetEntries()):
    if ievt%100 == 0: print('Processed {} events'.format(ievt))
    if ievt > 1000000: break
    tree.GetEntry(ievt)

    reg = reader.EvaluateRegression(methodname)[0]

    if verbose: 
      print("ievt = {0} reco_mjj = {1} reg = {2} reco_mjj_reg = {3}".format(ievt, reg_reco_mjj[0], reg, reg_reco_mjj[0]*reg))

    #if reg_reco_mjj[0] < 98 or reg_reco_mjj[0] > 143:
    h_reco_mjj.Fill(reg_reco_mjj[0])
    h_reco_mjj_reg.Fill(reg_reco_mjj[0]*reg)
    h_regwt.Fill(reg)
    h_mjjres_reg.Fill( (reg_reco_mjj[0]*reg - mbbNu)/mbbNu )
    h_mjjres.Fill( (reg_reco_mjj[0] - mbbNu)/mbbNu )

    if verbose:
      print('reg wt = {}'.format(reg))

  fout.Write()
  fout.Close()


def main():

  parser = argparse.ArgumentParser(description='Train or test')
  subparsers = parser.add_subparsers(dest='command')
  parser_train = subparsers.add_parser('train', help='train regression')
  parser_test  = subparsers.add_parser('test' , help='test  regression')

  args = parser.parse_args()
  if args.command == 'train':
    train()
  elif args.command == 'test':
    test()



if __name__ == "__main__":
  main()
