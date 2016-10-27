#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys, os, shutil, re, subprocess
import ROOT
import glob
from DataFormats.FWLite import Events, Handle
from rootpy.vector import LorentzVector
from rootpy.tree import Tree, TreeModel, FloatArrayCol, IntCol, FloatCol
from rootpy.io import root_open

subprocess.call("source /afs/cern.ch/user/d/devdatta/.tcshrc", shell=True)

class Event(TreeModel):
  nevents = IntCol()
  n_pv  = IntCol()
  n_Hbb = IntCol()
  m_HH  = FloatCol()
  n_AK8MatchedToHbb= IntCol()
  pt_Hbb                    = FloatArrayCol(10, length_name='n_Hbb')
  eta_Hbb                   = FloatArrayCol(10, length_name='n_Hbb')
  phi_Hbb                   = FloatArrayCol(10, length_name='n_Hbb')
  en_Hbb                    = FloatArrayCol(10, length_name='n_Hbb')
  pt_MatchedHbb             = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  eta_MatchedHbb            = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  phi_MatchedHbb            = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  en_MatchedHbb             = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  jec0_AK8MatchedToHbb      = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  pt_AK8MatchedToHbb        = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  eta_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  phi_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  en_AK8MatchedToHbb        = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  msoftdrop_AK8MatchedToHbb = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  nsv_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  sv0mass_AK8MatchedToHbb   = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  sv1mass_AK8MatchedToHbb   = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  nch_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  nmu_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  nel_AK8MatchedToHbb       = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  muenfr_AK8MatchedToHbb    = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')
  emenfr_AK8MatchedToHbb    = FloatArrayCol(10, length_name='n_AK8MatchedToHbb')  

inputs = [
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v3_B2GAnaFW_80X_V2p1/161013_161650/0000/B2GEDMNtuple_1.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v3_B2GAnaFW_80X_V2p1/161013_161650/0000/B2GEDMNtuple_2.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v3_B2GAnaFW_80X_V2p1/161013_161650/0000/B2GEDMNtuple_3.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161705/0000/B2GEDMNtuple_1.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161705/0000/B2GEDMNtuple_2.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161720/0000/B2GEDMNtuple_1.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161735/0000/B2GEDMNtuple_1.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161750/0000/B2GEDMNtuple_1.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161750/0000/B2GEDMNtuple_2.root', 
  '/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161013_161750/0000/B2GEDMNtuple_3.root', 
  ]

def getfiles(inputs):
 proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'ls', inputs ], 
   stdout = subprocess.PIPE, 
   stderr = subprocess.STDOUT )
 output = proc.communicate()[0]
 if proc.returncode != 0:
   print '>>>>Output:', output
   sys.exit(1)
 return output.splitlines()

def main():

  masses = [1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 4000, 4500]
  models = ['BulkGravTohhTohbbhbb', 'RadionTohhTohbbhbb']
  
  #inputs= glob.glob('/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/BulkGravTohhTohbbhbb_narrow_M-*/*/*/*/*/*/*')
  #files = getfiles(inputs)
  #print ">>>> >>>> Files", files

  h_npv         = Handle("std::int")

  h_genparticles= Handle("vector<reco::GenParticle>")

  h_ak8_jec0      = Handle("std::vector<float>")
  h_ak8_pt        = Handle("std::vector<float>")
  h_ak8_eta       = Handle("std::vector<float>")
  h_ak8_phi       = Handle("std::vector<float>")
  h_ak8_en        = Handle("std::vector<float>")
  h_ak8_msoftdrop = Handle("std::vector<float>")
  h_ak8_nsv       = Handle("std::vector<float>")
  h_ak8_sv0mass   = Handle("std::vector<float>")
  h_ak8_sv1mass   = Handle("std::vector<float>")
  h_ak8_nch       = Handle("std::vector<float>")
  h_ak8_nmu       = Handle("std::vector<float>")
  h_ak8_nel       = Handle("std::vector<float>")
  h_ak8_muen      = Handle("std::vector<float>")
  h_ak8_emenfr    = Handle("std::vector<float>")

  l_npv           = ("vertexInfo", "npv")

  l_genparticles  = ("filteredPrunedGenParticles")

  l_ak8_jec0      = ("jetsAK8Puppi", "jetAK8PuppijecFactor0")
  l_ak8_pt        = ("jetsAK8Puppi", "jetAK8PuppiPt")
  l_ak8_eta       = ("jetsAK8Puppi", "jetAK8PuppiEta")
  l_ak8_phi       = ("jetsAK8Puppi", "jetAK8PuppiPhi")
  l_ak8_en        = ("jetsAK8Puppi", "jetAK8PuppiE")
  l_ak8_msoftdrop = ("jetsAK8Puppi", "jetAK8PuppisoftDropMass")
  l_ak8_nsv       = ("jetsAK8Puppi", "jetAK8PuppinSV")
  l_ak8_sv0mass   = ("jetsAK8Puppi", "jetAK8PuppiSV0mass")
  l_ak8_sv1mass   = ("jetsAK8Puppi", "jetAK8PuppiSV1mass")
  l_ak8_nch       = ("jetsAK8Puppi", "jetAK8PuppichargedMultiplicity")
  l_ak8_nmu       = ("jetsAK8Puppi", "jetAK8PuppimuonMultiplicity")
  l_ak8_nel       = ("jetsAK8Puppi", "jetAK8PuppielectronMultiplicity")
  l_ak8_muen      = ("jetsAK8Puppi", "jetAK8PuppiMuonEnergy")
  l_ak8_emenfr    = ("jetsAK8Puppi", "jetAK8PuppichargedEmEnergyFrac")

  fout = root_open("TrainingTree_BulkGravTohhTohbbhbb_narrow.root", "RECREATE")
  tree = Tree("TrainingTree",
      "Training trees for AK8 jet mass regression",
      model=Event)

  for f in inputs:

    if f[0] == '#' or (not f.endswith('.root')): continue
    #if not 'BulkGravTohhTohbbhbb_narrow_M' in f: continue
    print ">>>> Opening file ", f

    events = Events('root://cms-xrd-global//'+f)
    nevents = 0
    for event in events:

      event.getByLabel(l_npv         , h_npv         )

      event.getByLabel(l_genparticles, h_genparticles)

      event.getByLabel(l_ak8_jec0      , h_ak8_jec0      )
      event.getByLabel(l_ak8_pt        , h_ak8_pt        )
      event.getByLabel(l_ak8_eta       , h_ak8_eta       )
      event.getByLabel(l_ak8_phi       , h_ak8_phi       )
      event.getByLabel(l_ak8_en        , h_ak8_en        )
      event.getByLabel(l_ak8_msoftdrop , h_ak8_msoftdrop )
      event.getByLabel(l_ak8_nsv       , h_ak8_nsv       )
      event.getByLabel(l_ak8_sv0mass   , h_ak8_sv0mass   )
      event.getByLabel(l_ak8_sv1mass   , h_ak8_sv1mass   )
      event.getByLabel(l_ak8_nch       , h_ak8_nch       )
      event.getByLabel(l_ak8_nmu       , h_ak8_nmu       )
      event.getByLabel(l_ak8_nel       , h_ak8_nel       )
      event.getByLabel(l_ak8_muen      , h_ak8_muen      )
      event.getByLabel(l_ak8_emenfr    , h_ak8_emenfr    )

      if len(h_npv.product()) < 1: continue

      if len(h_ak8_pt.product()) < 1: continue
     
      tree.n_pv = h_npv.product()[0]

      m_HH = 0
      n_Hbb = 0
      n_AK8MatchedToHbb = 0

      p4_Higgses = []

      if h_genparticles.product().size() > 0:
        for igen in range(0, h_genparticles.product().size()):

          gen_pdgid = h_genparticles.product()[igen].pdgId()  
          if gen_pdgid != 25: continue

          gen_ndau = h_genparticles.product()[igen].numberOfDaughters()     
          if gen_ndau != 2: continue

          dau0 = h_genparticles.product()[igen].daughter(0).pdgId()
          dau1 = h_genparticles.product()[igen].daughter(1).pdgId()
          if abs(dau0) != 5 or abs(dau1) != 5: continue

          gen_pt   = h_genparticles.product()[igen].pt()  
          gen_eta  = h_genparticles.product()[igen].eta()  
          gen_phi  = h_genparticles.product()[igen].phi()  
          gen_en   = h_genparticles.product()[igen].energy() 

          p4_hbb = ROOT.TLorentzVector()
          p4_hbb.SetPtEtaPhiE(gen_pt, gen_eta, gen_phi, gen_en)

          p4_Higgses.append(p4_hbb)

          tree.pt_Hbb [n_Hbb] = p4_hbb.Pt()
          tree.eta_Hbb[n_Hbb] = p4_hbb.Eta()
          tree.phi_Hbb[n_Hbb] = p4_hbb.Phi()
          tree.en_Hbb [n_Hbb] = p4_hbb.Energy()

          n_Hbb = n_Hbb + 1

          for ijet in range(0, len(h_ak8_pt.product())):
            jet_pt  = h_ak8_pt.product()[ijet]
            jet_eta = h_ak8_eta.product()[ijet]
            jet_phi = h_ak8_phi.product()[ijet]
            jet_en  = h_ak8_en.product()[ijet]
            
            p4_ak8 = ROOT.TLorentzVector()
            p4_ak8.SetPtEtaPhiE(jet_pt, jet_eta, jet_phi, jet_en)

            if jet_pt < 300: continue
            if abs(jet_eta) > 2.4: continue

            if p4_ak8.DeltaR(p4_hbb) > 0.3:
              continue

            ### Get variables for regression
            tree.pt_MatchedHbb [n_AK8MatchedToHbb] = p4_hbb.Pt()
            tree.eta_MatchedHbb[n_AK8MatchedToHbb] = p4_hbb.Eta()
            tree.phi_MatchedHbb[n_AK8MatchedToHbb] = p4_hbb.Phi()
            tree.en_MatchedHbb [n_AK8MatchedToHbb] = p4_hbb.Energy()

            tree.jec0_AK8MatchedToHbb      [n_AK8MatchedToHbb] = h_ak8_jec0.product()[ijet]
            tree.pt_AK8MatchedToHbb        [n_AK8MatchedToHbb] = p4_ak8.Pt()
            tree.eta_AK8MatchedToHbb       [n_AK8MatchedToHbb] = p4_ak8.Eta()
            tree.phi_AK8MatchedToHbb       [n_AK8MatchedToHbb] = p4_ak8.Phi()
            tree.en_AK8MatchedToHbb        [n_AK8MatchedToHbb] = p4_ak8.Energy()
            tree.msoftdrop_AK8MatchedToHbb [n_AK8MatchedToHbb] = h_ak8_msoftdrop.product()[ijet] 
            tree.nsv_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nsv    .product()[ijet] 
            tree.sv0mass_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_sv0mass.product()[ijet] 
            tree.sv1mass_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_sv1mass.product()[ijet] 
            tree.nch_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nch    .product()[ijet] 
            tree.nmu_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nmu    .product()[ijet] 
            tree.nel_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nel    .product()[ijet] 
            tree.muenfr_AK8MatchedToHbb    [n_AK8MatchedToHbb] = h_ak8_muen   .product()[ijet]/p4_ak8.Energy() 
            tree.emenfr_AK8MatchedToHbb    [n_AK8MatchedToHbb] = h_ak8_emenfr .product()[ijet] 

            n_AK8MatchedToHbb = n_AK8MatchedToHbb + 1

      if len(p4_Higgses) > 1:
        m_HH = (p4_Higgses[0] + p4_Higgses[1]).Mag()

      tree.n_Hbb = n_Hbb
      tree.m_HH = m_HH
      tree.n_AK8MatchedToHbb = n_AK8MatchedToHbb
      tree.nevents = nevents
      tree.fill()

      nevents += 1

  tree.Write()
  fout.Close()

if __name__ == "__main__":
  main()
