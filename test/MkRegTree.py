#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys, os, shutil, re, subprocess
import ROOT
import glob
import numpy as np
from array import array
from inputFiles import inputs
from DataFormats.FWLite import Events, Handle

MAX=1000

class RegisterTree:
  nevents                   = np.zeros(1, dtype=int)
  n_pv                      = np.zeros(1, dtype=int)
  n_Hbb                     = np.zeros(1, dtype=int)
  m_HH                      = np.zeros(1, dtype=int)
  n_AK8MatchedToHbb         = np.zeros(1, dtype=int)
  pt_Hbb                    = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  eta_Hbb                   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  phi_Hbb                   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  en_Hbb                    = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  pt_MatchedHbb             = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  eta_MatchedHbb            = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  phi_MatchedHbb            = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  en_MatchedHbb             = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  jec0_AK8MatchedToHbb      = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  pt_AK8MatchedToHbb        = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  eta_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  phi_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  en_AK8MatchedToHbb        = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  mass_AK8MatchedToHbb      = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  msdCHS_AK8MatchedToHbb    = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  msdPuppi_AK8MatchedToHbb  = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  nsv_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  sv0mass_AK8MatchedToHbb   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  sv1mass_AK8MatchedToHbb   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  nch_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  nmu_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  nel_AK8MatchedToHbb       = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  muenfr_AK8MatchedToHbb    = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float) 
  emenfr_AK8MatchedToHbb    = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)   
  DoubleB_AK8MatchedToHbb   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)
  tau1CHS_AK8MatchedToHbb   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)
  tau2CHS_AK8MatchedToHbb   = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)
  tau1Puppi_AK8MatchedToHbb = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)
  tau2Puppi_AK8MatchedToHbb = array('f', MAX*[0.]) #np.zeros(MAX, dtype=float)

  def __init__(self, tree): 
    tree.Branch('n_pv'               ,self.n_pv               ,'n_pv/I'               )
    tree.Branch('n_Hbb'              ,self.n_Hbb              ,'n_Hbb/I'              )
    tree.Branch('m_HH'               ,self.m_HH               ,'m_HH/I'               )
    tree.Branch('n_AK8MatchedToHbb'  ,self.n_AK8MatchedToHbb  ,'n_AK8MatchedToHbb/I'  )
    tree.Branch('pt_Hbb'             ,self.pt_Hbb             ,'pt_Hbb[n_Hbb]/F'      )
    tree.Branch('eta_Hbb'            ,self.eta_Hbb            ,'eta_Hbb[n_Hbb]/F'     )
    tree.Branch('phi_Hbb'            ,self.phi_Hbb            ,'phi_Hbb[n_Hbb]/F'     )
    tree.Branch('en_Hbb'             ,self.en_Hbb             ,'en_Hbb[n_Hbb]/F'      )

    tree.Branch('pt_MatchedHbb'             ,self.pt_MatchedHbb             ,'pt_MatchedHbb[n_AK8MatchedToHbb]/F'             )     
    tree.Branch('eta_MatchedHbb'            ,self.eta_MatchedHbb            ,'eta_MatchedHbb[n_AK8MatchedToHbb]/F'            )     
    tree.Branch('phi_MatchedHbb'            ,self.phi_MatchedHbb            ,'phi_MatchedHbb[n_AK8MatchedToHbb]/F'            )     
    tree.Branch('en_MatchedHbb'             ,self.en_MatchedHbb             ,'en_MatchedHbb[n_AK8MatchedToHbb]/F'             )     
    tree.Branch('jec0_AK8MatchedToHbb'      ,self.jec0_AK8MatchedToHbb      ,'jec0_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'      )     

    tree.Branch('pt_AK8MatchedToHbb'        ,self.pt_AK8MatchedToHbb        ,'pt_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'        )     
    tree.Branch('eta_AK8MatchedToHbb'       ,self.eta_AK8MatchedToHbb       ,'eta_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       )     
    tree.Branch('phi_AK8MatchedToHbb'       ,self.phi_AK8MatchedToHbb       ,'phi_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       )     
    tree.Branch('en_AK8MatchedToHbb'        ,self.en_AK8MatchedToHbb        ,'en_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'        )     
    tree.Branch('mass_AK8MatchedToHbb'      ,self.mass_AK8MatchedToHbb      ,'mass_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'      )     
    tree.Branch('msdCHS_AK8MatchedToHbb'    ,self.msdCHS_AK8MatchedToHbb    ,'msdCHS_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'    ) 
    tree.Branch('msdPuppi_AK8MatchedToHbb'  ,self.msdPuppi_AK8MatchedToHbb  ,'msdPuppi_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'  ) 

    tree.Branch('nsv_AK8MatchedToHbb'       ,self.nsv_AK8MatchedToHbb       ,'nsv_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       ) 
    tree.Branch('sv0mass_AK8MatchedToHbb'   ,self.sv0mass_AK8MatchedToHbb   ,'sv0mass_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'   ) 
    tree.Branch('sv1mass_AK8MatchedToHbb'   ,self.sv1mass_AK8MatchedToHbb   ,'sv1mass_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'   ) 

    tree.Branch('nch_AK8MatchedToHbb'       ,self.nch_AK8MatchedToHbb       ,'nch_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       ) 
    tree.Branch('nmu_AK8MatchedToHbb'       ,self.nmu_AK8MatchedToHbb       ,'nmu_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       ) 
    tree.Branch('nel_AK8MatchedToHbb'       ,self.nel_AK8MatchedToHbb       ,'nel_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'       ) 
    tree.Branch('muenfr_AK8MatchedToHbb'    ,self.muenfr_AK8MatchedToHbb    ,'muenfr_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'    ) 
    tree.Branch('emenfr_AK8MatchedToHbb'    ,self.emenfr_AK8MatchedToHbb    ,'emenfr_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'    ) 

    tree.Branch('DoubleB_AK8MatchedToHbb'   ,self.DoubleB_AK8MatchedToHbb   ,'DoubleB_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'   ) 
    tree.Branch('tau1CHS_AK8MatchedToHbb'   ,self.tau1CHS_AK8MatchedToHbb   ,'tau1CHS_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'   ) 
    tree.Branch('tau2CHS_AK8MatchedToHbb'   ,self.tau2CHS_AK8MatchedToHbb   ,'tau2CHS_AK8MatchedToHbb[n_AK8MatchedToHbb]/F'   ) 
    tree.Branch('tau1Puppi_AK8MatchedToHbb' ,self.tau1Puppi_AK8MatchedToHbb ,'tau1Puppi_AK8MatchedToHbb[n_AK8MatchedToHbb]/F' ) 
    tree.Branch('tau2Puppi_AK8MatchedToHbb' ,self.tau2Puppi_AK8MatchedToHbb ,'tau2Puppi_AK8MatchedToHbb[n_AK8MatchedToHbb]/F' ) 

def main():

  maxEvents = sys.argv[1] if len(sys.argv) > 1 else -1

  masses = [1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 4000, 4500]
  models = ['BulkGravTohhTohbbhbb', 'RadionTohhTohbbhbb']

  h_npv         = Handle("std::int")

  h_genparticles= Handle("vector<reco::GenParticle>")

  h_ak8_jec0      = Handle("std::vector<float>")
  h_ak8_pt        = Handle("std::vector<float>")
  h_ak8_eta       = Handle("std::vector<float>")
  h_ak8_phi       = Handle("std::vector<float>")
  h_ak8_en        = Handle("std::vector<float>")
  h_ak8_msdCHS    = Handle("std::vector<float>")
  h_ak8_msdPuppi  = Handle("std::vector<float>")
  h_ak8_nsv       = Handle("std::vector<float>")
  h_ak8_sv0mass   = Handle("std::vector<float>")
  h_ak8_sv1mass   = Handle("std::vector<float>")
  h_ak8_nch       = Handle("std::vector<float>")
  h_ak8_nmu       = Handle("std::vector<float>")
  h_ak8_nel       = Handle("std::vector<float>")
  h_ak8_muen      = Handle("std::vector<float>")
  h_ak8_emenfr    = Handle("std::vector<float>")
  h_ak8_DoubleB   = Handle("std::vector<float>")  
  h_ak8_tau1CHS   = Handle("std::vector<float>")  
  h_ak8_tau2CHS   = Handle("std::vector<float>")  
  h_ak8_tau1Puppi = Handle("std::vector<float>")  
  h_ak8_tau2Puppi = Handle("std::vector<float>") 

  l_npv           = ("vertexInfo", "npv")

  l_genparticles  = ("filteredPrunedGenParticles")

  l_ak8_jec0      = ("jetsAK8CHS", "jetAK8CHSjecFactor0")
  l_ak8_pt        = ("jetsAK8CHS", "jetAK8CHSPt")
  l_ak8_eta       = ("jetsAK8CHS", "jetAK8CHSEta")
  l_ak8_phi       = ("jetsAK8CHS", "jetAK8CHSPhi")
  l_ak8_en        = ("jetsAK8CHS", "jetAK8CHSE")
  l_ak8_msdCHS    = ("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS")
  l_ak8_msdPuppi  = ("jetsAK8CHS", "jetAK8CHSsoftDropMassPuppi")
  l_ak8_nsv       = ("jetsAK8CHS", "jetAK8CHSnSV")
  l_ak8_sv0mass   = ("jetsAK8CHS", "jetAK8CHSSV0mass")
  l_ak8_sv1mass   = ("jetsAK8CHS", "jetAK8CHSSV1mass")
  l_ak8_nch       = ("jetsAK8CHS", "jetAK8CHSchargedMultiplicity")
  l_ak8_nmu       = ("jetsAK8CHS", "jetAK8CHSmuonMultiplicity")
  l_ak8_nel       = ("jetsAK8CHS", "jetAK8CHSelectronMultiplicity")
  l_ak8_muen      = ("jetsAK8CHS", "jetAK8CHSMuonEnergy")
  l_ak8_emenfr    = ("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac")

  l_ak8_DoubleB   = ("jetsAK8CHS", "jetAK8CHSDoubleBAK8")
  l_ak8_tau1CHS   = ("jetsAK8CHS", "jetAK8CHStau1CHS"   )
  l_ak8_tau2CHS   = ("jetsAK8CHS", "jetAK8CHStau2CHS"   )
  l_ak8_tau1Puppi = ("jetsAK8CHS", "jetAK8CHStau1Puppi" )
  l_ak8_tau2Puppi = ("jetsAK8CHS", "jetAK8CHStau2Puppi" )

  fout = ROOT.TFile("TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root", "RECREATE")
  fout.cd()
  tree = ROOT.TTree("TrainingTree","TrainingTree")
  treereg = RegisterTree(tree)

  nevents = 0
  for f in inputs:

    if f[0] == '#' or (not f.endswith('.root')): continue
    f = 'root://eoscms.cern.ch/'+f
    #if not 'BulkGravTohhTohbbhbb_narrow_M' in f: continue
    print ">>>> Opening file ", f

    events = Events(f)
    for event in events:

      if maxEvents > 0 and nevents > maxEvents: break

      event.getByLabel(l_npv         , h_npv         )

      event.getByLabel(l_genparticles, h_genparticles)

      event.getByLabel(l_ak8_jec0      , h_ak8_jec0      )
      event.getByLabel(l_ak8_pt        , h_ak8_pt        )
      event.getByLabel(l_ak8_eta       , h_ak8_eta       )
      event.getByLabel(l_ak8_phi       , h_ak8_phi       )
      event.getByLabel(l_ak8_en        , h_ak8_en        )
      event.getByLabel(l_ak8_msdCHS    , h_ak8_msdCHS    )
      event.getByLabel(l_ak8_msdPuppi  , h_ak8_msdPuppi  )
      event.getByLabel(l_ak8_nsv       , h_ak8_nsv       )
      event.getByLabel(l_ak8_sv0mass   , h_ak8_sv0mass   )
      event.getByLabel(l_ak8_sv1mass   , h_ak8_sv1mass   )
      event.getByLabel(l_ak8_nch       , h_ak8_nch       )
      event.getByLabel(l_ak8_nmu       , h_ak8_nmu       )
      event.getByLabel(l_ak8_nel       , h_ak8_nel       )
      event.getByLabel(l_ak8_muen      , h_ak8_muen      )
      event.getByLabel(l_ak8_emenfr    , h_ak8_emenfr    )
      event.getByLabel(l_ak8_DoubleB   , h_ak8_DoubleB   )
      event.getByLabel(l_ak8_tau1CHS   , h_ak8_tau1CHS   )
      event.getByLabel(l_ak8_tau2CHS   , h_ak8_tau2CHS   )
      event.getByLabel(l_ak8_tau1Puppi , h_ak8_tau1Puppi )
      event.getByLabel(l_ak8_tau2Puppi , h_ak8_tau2Puppi )

      if len(h_npv.product()) < 1: continue

      if len(h_ak8_pt.product()) < 1: continue

      treereg.n_pv[0] = h_npv.product()[0]

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

          dau0_pt   = h_genparticles.product()[igen].daughter(0).pt()  
          dau0_eta  = h_genparticles.product()[igen].daughter(0).eta()  
          dau0_phi  = h_genparticles.product()[igen].daughter(0).phi()  
          dau0_en   = h_genparticles.product()[igen].daughter(0).energy() 

          dau1_pt   = h_genparticles.product()[igen].daughter(1).pt()  
          dau1_eta  = h_genparticles.product()[igen].daughter(1).eta()  
          dau1_phi  = h_genparticles.product()[igen].daughter(1).phi()  
          dau1_en   = h_genparticles.product()[igen].daughter(1).energy() 

          p4_hbb = ROOT.TLorentzVector()
          p4_hbb.SetPtEtaPhiE(gen_pt, gen_eta, gen_phi, gen_en)

          p4_dau0 = ROOT.TLorentzVector()
          p4_dau0.SetPtEtaPhiE(dau0_pt, dau0_eta, dau0_phi, dau0_en)

          p4_dau1 = ROOT.TLorentzVector()
          p4_dau1.SetPtEtaPhiE(dau1_pt, dau1_eta, dau1_phi, dau1_en)

          p4_Higgses.append(p4_hbb)

          treereg.pt_Hbb [n_Hbb] = p4_hbb.Pt()
          treereg.eta_Hbb[n_Hbb] = p4_hbb.Eta()
          treereg.phi_Hbb[n_Hbb] = p4_hbb.Phi()
          treereg.en_Hbb [n_Hbb] = p4_hbb.Energy()

          n_Hbb = n_Hbb + 1

          ### For each H->bb find a matched jet
          for ijet in range(0, len(h_ak8_pt.product())):
            jet_pt  = h_ak8_pt.product()[ijet]
            jet_eta = h_ak8_eta.product()[ijet]
            jet_phi = h_ak8_phi.product()[ijet]
            jet_en  = h_ak8_en.product()[ijet]

            p4_ak8 = ROOT.TLorentzVector()
            p4_ak8.SetPtEtaPhiE(jet_pt, jet_eta, jet_phi, jet_en)

            if jet_pt < 300: continue
            if abs(jet_eta) > 2.4: continue
            if p4_ak8.DeltaR(p4_hbb)  > 0.3: continue
            if p4_ak8.DeltaR(p4_dau0) > 0.8: continue
            if p4_ak8.DeltaR(p4_dau1) > 0.8: continue

            ### Get variables for regression
            treereg.pt_MatchedHbb [n_AK8MatchedToHbb] = p4_hbb.Pt()
            treereg.eta_MatchedHbb[n_AK8MatchedToHbb] = p4_hbb.Eta()
            treereg.phi_MatchedHbb[n_AK8MatchedToHbb] = p4_hbb.Phi()
            treereg.en_MatchedHbb [n_AK8MatchedToHbb] = p4_hbb.Energy()

            treereg.jec0_AK8MatchedToHbb      [n_AK8MatchedToHbb] = h_ak8_jec0.product()[ijet]
            treereg.pt_AK8MatchedToHbb        [n_AK8MatchedToHbb] = p4_ak8.Pt()
            treereg.eta_AK8MatchedToHbb       [n_AK8MatchedToHbb] = p4_ak8.Eta()
            treereg.phi_AK8MatchedToHbb       [n_AK8MatchedToHbb] = p4_ak8.Phi()
            treereg.en_AK8MatchedToHbb        [n_AK8MatchedToHbb] = p4_ak8.Energy()
            treereg.mass_AK8MatchedToHbb      [n_AK8MatchedToHbb] = p4_ak8.Mag()
            treereg.msdCHS_AK8MatchedToHbb    [n_AK8MatchedToHbb] = h_ak8_msdCHS.product()[ijet] 
            treereg.msdPuppi_AK8MatchedToHbb  [n_AK8MatchedToHbb] = h_ak8_msdPuppi.product()[ijet] 
            treereg.nsv_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nsv    .product()[ijet] 
            treereg.sv0mass_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_sv0mass.product()[ijet] 
            treereg.sv1mass_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_sv1mass.product()[ijet] 
            treereg.nch_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nch    .product()[ijet] 
            treereg.nmu_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nmu    .product()[ijet] 
            treereg.nel_AK8MatchedToHbb       [n_AK8MatchedToHbb] = h_ak8_nel    .product()[ijet] 
            treereg.muenfr_AK8MatchedToHbb    [n_AK8MatchedToHbb] = h_ak8_muen   .product()[ijet]/p4_ak8.Energy() 
            treereg.emenfr_AK8MatchedToHbb    [n_AK8MatchedToHbb] = h_ak8_emenfr .product()[ijet] 
            treereg.DoubleB_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_DoubleB.product()[ijet] 

            treereg.tau1CHS_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_tau1CHS  .product()[ijet] 
            treereg.tau2CHS_AK8MatchedToHbb   [n_AK8MatchedToHbb] = h_ak8_tau2CHS  .product()[ijet] 
            treereg.tau1Puppi_AK8MatchedToHbb [n_AK8MatchedToHbb] = h_ak8_tau1Puppi.product()[ijet] 
            treereg.tau2Puppi_AK8MatchedToHbb [n_AK8MatchedToHbb] = h_ak8_tau2Puppi.product()[ijet] 

            n_AK8MatchedToHbb = n_AK8MatchedToHbb + 1

      if len(p4_Higgses) > 1:
        treereg.m_HH[0] = (p4_Higgses[0] + p4_Higgses[1]).Mag()

      treereg.n_Hbb[0] = n_Hbb
      treereg.n_AK8MatchedToHbb[0] = n_AK8MatchedToHbb
      treereg.nevents[0] = nevents
      tree.Fill()

      nevents += 1

  fout.Write()
  fout.Close()

if __name__ == "__main__":
    main()
