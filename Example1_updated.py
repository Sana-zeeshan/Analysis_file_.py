#!/usr/bin/env python

import math
import sys
import time
import ROOT
import numpy as np
from array import array

from definitions_final import delta_phi, delta_r, calculate_event_shape, define_histograms, define_tree, calculate_fox_wolfram, transverse_mass
# Initialize histograms and tree
histograms = define_histograms()
my_tree, branches = define_tree()

try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Example1.py input_file")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass
# Create a ROOT file to save histograms
file = ROOT.TFile("output_file.root", "RECREATE")

inputFile = sys.argv[1]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()
#numberOfEntries = 4
# Get pointers to branches used in this analysis
branchWeight   = treeReader.UseBranch("Weight")
branchEvent    = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")
branchScalarHT = treeReader.UseBranch("ScalarHT")
branchMET = treeReader.UseBranch("MissingET")


#n_leps = 4
#lep_pt_br = [array('f', [0.]) for _ in range(n_leps)]  # Array to hold PT values of 4 leptons for lepton branch

# Create branches dynamically using a loop for leptons
#for i in range(n_leps):
#    my_tree.Branch(f"lep{i+1}_pt", lep_pt_br[i], f"lep{i+1}_pt/F")

# Book histograms
#########histJetsSize = ROOT.TH1F("jet_Size", "jet Size", 18, 0.0, 17.0)

#histbJetsSize = ROOT.TH1F("bjets_Size", "bjets Size", 17, 0.0, 17.0)
#histScalarHT = ROOT.TH1F("Scalar_HT", "Scalar HT", 100, 100.0, 1800.0)
#histMET = ROOT.TH1F("MET", "MET", 100, 0.0, 600.0)
bjets_hPt_l = []
bjets_hEta_l = []
bjets_hPhi_l = []
bjets_hMass_l = []
bjets_hdPhi_l = []
bjets_hdR_l = []

lep_hPt_l = []
lep_hEta_l = []
lep_hPhi_l = []
for i in range (4):
    bjets_Pt_hName = f"bjet{i+1}_hPt"
    hist_bjet_Pt = ROOT.TH1F(bjets_Pt_hName, bjets_Pt_hName, 100, 0, 400)
    bjets_hPt_l.append(hist_bjet_Pt)

    bjets_Eta_hName = f"bjet{i+1}_hEta"
    hist_bjet_Eta = ROOT.TH1F(bjets_Eta_hName, bjets_Eta_hName, 100, -4.0, 4.0)
    bjets_hEta_l.append(hist_bjet_Eta)

    bjets_Phi_hName = f"bjet{i+1}_hPhi"
    hist_bjet_Phi = ROOT.TH1F(bjets_Phi_hName, bjets_Phi_hName, 100, 0.0, 3.20)
    bjets_hPhi_l.append(hist_bjet_Phi)

    bjets_Mass_hName = f"bjet{i+1}_hMass"
    hist_bjet_Mass = ROOT.TH1F(bjets_Mass_hName, bjets_Mass_hName, 100, 0.0, 200.0)
    bjets_hMass_l.append(hist_bjet_Mass)
    for j in range (i+1, 4):
        bjets_dPhi_hName = f"bjet{i}_{j}_dPhi"
        hist_bjet_dPhi = ROOT.TH1F(bjets_dPhi_hName, bjets_dPhi_hName, 100, -3.15, 3.15)
        bjets_hdPhi_l.append(hist_bjet_dPhi)


    lep_Pt_hName = f"lep{i+1}_hPt"
    hist_lep_Pt = ROOT.TH1F(lep_Pt_hName, lep_Pt_hName, 100, 0, 400)
    lep_hPt_l.append(hist_lep_Pt)

    lep_Eta_hName = f"lep{i+1}_hEta"
    hist_lep_Eta = ROOT.TH1F(lep_Eta_hName, lep_Eta_hName, 100, -4.0, 4.0)
    lep_hEta_l.append(hist_lep_Eta)

    lep_Phi_hName = f"lep{i+1}_hPhi"
    hist_lep_Phi = ROOT.TH1F(lep_Phi_hName, lep_Phi_hName, 100, 0.0, 3.2)
    lep_hPhi_l.append(hist_lep_Phi)


#histMass = ROOT.TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0)

# Loop over all events
for entry in range(0, numberOfEntries):
#for entry in range(0, 2):
  if (entry+1)%100 == 0:
    print (' ... processed {} events ...'.format(entry+100))

  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)
#  print("No of entries: ",entry)

  ## main MC event weight
  w =  branchEvent[0].Weight
  for met in range(0,branchMET.GetEntries()):
    histograms["histMET"].Fill(branchMET.At(met).MET,w)
#    print('met information: ',branchMET.At(met).Phi)
    branches["MET_"][0] = branchMET.At(met).MET
  for sHT in range(0,branchScalarHT.GetEntries()):
    histograms["histScalarHT"].Fill(branchScalarHT.At(sHT).HT,w)
    branches["SHT"][0] = branchScalarHT.At(sHT).HT

  # If event contains at least 4 jet
  histograms["histJetsSize"].Fill(branchJet.GetEntries(), w)
  branches["jets_size"][0] = branchJet.GetEntries()
  my_tree.Fill()
  if branchJet.GetEntries() > 3:
    momentum_vectors_jets = [[branchJet.At(jet).P4().Px(), branchJet.At(jet).P4().Py(), branchJet.At(jet).P4().Pz()] for jet in range (0, branchJet.GetEntries())]
    sphericity, aplanarity, circularity = calculate_event_shape(momentum_vectors_jets)
    branches["sphericity_jets"][0] = sphericity
    branches["aplanarity_jets"][0] = aplanarity
    branches["circularity_jets"][0] = circularity

#    branches["sphericity_jets"][0], branches["aplanarity_jets"][0], branches["circularity_jets"][0] = calculate_event_shape(momentum_vectors_jets)

    ## 0 - Loose , 1 - Medium, 2 - Tight
    wp = 1
    bjets_list = []
    dphi_values = []
    dR_values = []
    for bj in range(0,branchJet.GetEntries()):
        #    BtagOk = ( jet1.BTag & (1 << wp) )
        BtagOk = ( branchJet.At(bj).BTag )
        pt = branchJet.At(bj).PT
        eta = abs(branchJet.At(bj).Eta)
        phi = branchJet.At(bj).Phi
#    print("btag>>>>>>>>>>>>>>>>:",BtagOk)
    # Plot bjets transverse momentum
        if (BtagOk and pt > 30. and eta < 5.):
            bjets_list.append(branchJet.At(bj))
    histograms["histbJetsSize"].Fill(len(bjets_list), w)
    branches["bjets_size"][0] = len(bjets_list)
#    my_tree.Fill()
    print("Jets no: ",branchJet.GetEntries())
    if len(bjets_list) > 3:
        momentum_vectors_bjets = [[bjet.P4().Px(), bjet.P4().Py(), bjet.P4().Pz()] for bjet in bjets_list]
        sphericity, aplanarity, circularity = calculate_event_shape(momentum_vectors_bjets)
        branches["sphericity_bjets"][0] = sphericity
        branches["aplanarity_bjets"][0] = aplanarity
        # Compute Fox-Wolfram Moments
        max_order = 4
        bjets_fox_wolfram_moments = calculate_fox_wolfram(momentum_vectors_bjets, max_order=max_order)

    # Update Fox-Wolfram branches
        for l, moment in enumerate(bjets_fox_wolfram_moments):
            branches[f"bjets_fox_wolfram_H{l}"][0] = moment
        for hist in bjets_hPt_l:
            idx = bjets_hPt_l.index(hist)
            hist.Fill(branchJet.At(idx).PT, w)
# Fill Branch of tree
            #bjet_pt_br[idx][0] = branchJet.At(idx).PT
            branches[f"bjet{idx+1}_pt_br"][0] = branchJet.At(idx).PT
#        my_tree.Fill()   

        for hist in bjets_hEta_l:
            idx = bjets_hEta_l.index(hist)
            hist.Fill(branchJet.At(idx).Eta, w)
        bj_pair_count = 0
        for idx, hist in enumerate(bjets_hPhi_l):
#            idx = bjets_hPhi_l.index(hist)
            print("bjets index in phi loop: ", idx);
            hist.Fill(branchJet.At(idx).Phi, w)
            for idx2 in range(idx+1, len(bjets_hPhi_l)):
                dphi = delta_phi(branchJet.At(idx).Phi, branchJet.At(idx2).Phi)
                dphi_values.append(dphi)
                branches[f"bjet_dphi_br{idx}_{idx2}"][0] = dphi
                print("idx1: ", idx,",  idx2: ", idx2, ",  dPhi:  ",dphi)
                dR = delta_r(branchJet.At(idx).Eta, branchJet.At(idx2).Eta, branchJet.At(idx).Phi, branchJet.At(idx2).Phi)
                dR_values.append(dR)
                branches[f"bjet_dr_br{idx}_{idx2}"][0] = dR
                print("deta R:  ", dR)
                bj_pair_count += 1
        for idx_dphi, (hist_dphi, dphi_val) in enumerate(zip(bjets_hdPhi_l, dphi_values)):
                hist_dphi.Fill(dphi_val, w)



        for hist in bjets_hMass_l:
            idx = bjets_hMass_l.index(hist)
            hist.Fill(branchJet.At(idx).Mass, w)




  # If event contains at least 4 leptons
    leptons_list = []
    print(f"Event {entry}: Electrons = {branchElectron.GetEntries()}, Muons = {branchMuon.GetEntries()}")

    if (branchElectron.GetEntries() + branchMuon.GetEntries() ) > 0:
    # Take first two electrons
        for e in range(0,branchElectron.GetEntries()):
            leptons_list.append(branchElectron.At(e))
        for m in range(0,branchMuon.GetEntries()):
            leptons_list.append(branchMuon.At(m))
            
            
        leptons_list.sort(key=lambda lep: lep.PT, reverse=True)
        print(f"Event {entry}: Found {len(leptons_list)} leptons after sorting")

    # Compute transverse mass (mT) for the leading lepton
    
    
    if len(leptons_list) > 0 and branchMET.GetEntries() > 0:
        leading_lepton = leptons_list[0]  # Highest pT lepton
        for met in range(0,branchMET.GetEntries()):
        	mt_value = transverse_mass(leading_lepton.PT, leading_lepton.Phi, branchMET.At(0).MET, branchMET.At(0).Phi)
        	branches["mT"][0] = float(mt_value)
#        	print("mt_value; ", mt_value, "lep pt: ", leading_lepton.PT, ", ", leading_lepton.Phi, ", ", branchMET.At(0).MET, ", ", branchMET.At(0).Phi)
        
        
       #leptons_list.sort()
        momentum_vectors_lep = [[lep.P4().Px(), lep.P4().Py(), lep.P4().Pz()] for lep in leptons_list]
        sphericity, aplanarity, circularity = calculate_event_shape(momentum_vectors_lep)
        branches["sphericity_leps"][0] = sphericity
        branches["aplanarity_leps"][0] = aplanarity
    # Compute Fox-Wolfram Moments
        max_order = 4
        leps_fox_wolfram_moments = calculate_fox_wolfram(momentum_vectors_lep, max_order=max_order)

    # Update Fox-Wolfram branches
        for l, moment in enumerate(leps_fox_wolfram_moments):
            branches[f"leps_fox_wolfram_H{l}"][0] = moment
    if len(leptons_list) == 0:
        print(f"Skipping lepton filling: No leptons in this event.")
    else:
        for idx, lep in enumerate(lep_hPt_l):
          if idx < len(leptons_list):
            lep.Fill(leptons_list[idx].PT, w)
          
# Fill leptons brach
            branches[f"lep{idx}_pt_br"][0] = leptons_list[idx].PT
#        my_tree.Fill()
    if len(leptons_list) == 0:
        print(f"Skipping lepton filling: No leptons in this event.")
    else:   
        for idx, lep in enumerate(lep_hEta_l):
            if idx < len(leptons_list):
             lep.Fill(leptons_list[idx].Eta, w)
            
            
    if len(leptons_list) == 0:
        print(f"Skipping lepton filling: No leptons in this event.")
    else:
        for idx, lep in enumerate(lep_hPhi_l):
            if idx < len(leptons_list):
             lep.Fill(leptons_list[idx].Phi, w)
#      print ("leptons list", leptons_list)
#      print ("lepton index: ", idx)
#    print("Electron no: ",branchElectron.GetEntries())
#    print("Muon no: ",branchMuon.GetEntries())
    if len(bjets_list) >= 2 and len(leptons_list) >= 2:
        for i in range(2):
            branches["br_dphi_bjet_lep_leading"][i][0] = delta_phi(bjets_list[i].Phi, leptons_list[i].Phi)
            branches["br_dr_bjet_lep_leading"][i][0] = delta_r(bjets_list[i].Eta, leptons_list[i].Eta, bjets_list[i].Phi, leptons_list[i].Phi)

    # Plot their invariant mass
 #   histMass.Fill(((elec1.P4()) + (elec2.P4())).M())
#    my_tree.Fill()
#Write tree
#my_tree.Fill()
my_tree.Write()
# Write the histograms to the ROOT file
for hist in histograms.values():
    hist.Write()
file.Write()

# Close the file
file.Close()


