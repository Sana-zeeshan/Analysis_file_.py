import ROOT
import numpy as np
from array import array
import math

# Function definitions
def delta_phi(phi1, phi2):
    dphi = phi2 - phi1
    if abs(dphi) > math.pi:
        dphi = 2 * math.pi - abs(dphi)
    return dphi

def delta_r(eta1, eta2, phi1, phi2):
    deta = eta2 - eta1
    dphi = delta_phi(phi1, phi2)  # Reuse delta_phi function
    return math.sqrt(dphi**2 + deta**2)

def transverse_mass(pt, phi, met, met_phi):
    """Compute transverse mass (m_T)"""
    dp = delta_phi(phi, met_phi)
    return math.sqrt(2 * pt * met * (1 - math.cos(dp)))

def calculate_event_shape(momentum_vectors):
    """
    Calculate sphericity, aplanarity, and circularity for a given list of momentum vectors.
    
    Args:
        momentum_vectors (list of list): List of [px, py, pz] for particles.
    
    Returns:
        tuple: (sphericity, aplanarity, circularity)
    """
    if len(momentum_vectors) == 0:  # Need at least 2 particles for meaningful calculations
        return -1, -1, -1

    # Calculate the full 3D sphericity tensor
    S = np.zeros((3, 3))
    norm_factor_3D = sum(np.linalg.norm(p) ** 2 for p in momentum_vectors)
    for p in momentum_vectors:
        for i in range(3):
            for j in range(3):
                S[i, j] += p[i] * p[j]
    S /= norm_factor_3D

    # Diagonalize to find eigenvalues for 3D sphericity tensor
    eigenvalues_3D = np.linalg.eigvalsh(S)
    eigenvalues_3D = sorted(eigenvalues_3D, reverse=True)  # Descending order

    sphericity = 1.5 * (eigenvalues_3D[1] + eigenvalues_3D[2])
    aplanarity = 1.5 * eigenvalues_3D[2]

    # Calculate the se 2D momentum tensor for circularity
    T = np.zeros((2, 2))
    norm_factor_2D = sum(p[0]**2 + p[1]**2 for p in momentum_vectors)
    for px, py, _ in momentum_vectors:  # Only consider transverse components
        T[0, 0] += px * px
        T[0, 1] += px * py
        T[1, 0] += py * px
        T[1, 1] += py * py
    T /= norm_factor_2D

    # Diagonalize to find eigenvalues for transverse tensor
    eigenvalues_2D = np.linalg.eigvalsh(T)
    circularity = 2 * min(eigenvalues_2D)

    return sphericity, aplanarity, circularity

from scipy.special import legendre
max_order=4
def calculate_fox_wolfram(momentum_vectors, max_order=4):
    """
    Calculate Fox-Wolfram Moments for a given set of momentum vectors.
    
    Args:
        momentum_vectors (list of list): List of [px, py, pz] for particles.
        max_order (int): Maximum order of Legendre polynomials to calculate (default: 4).
    
    Returns:
        list: Fox-Wolfram moments [H_0, H_1, ..., H_max_order].
    """
    if len(momentum_vectors) < 2:  # Need at least 2 particles for meaningful calculation
        return [0.0] * (max_order + 1)

    # Compute the total momentum magnitude
    total_momentum = sum(np.linalg.norm(p) for p in momentum_vectors)

    # Initialize moments
    moments = [0.0] * (max_order + 1)

    # Calculate pairwise contributions
    for i, p_i in enumerate(momentum_vectors):
        norm_i = np.linalg.norm(p_i)
        for j, p_j in enumerate(momentum_vectors):
            if i >= j:  # Avoid double counting
                continue

            norm_j = np.linalg.norm(p_j)
            cos_theta_ij = np.dot(p_i, p_j) / (norm_i * norm_j)  # Cosine of the angle

            # Add contributions to each moment
            for l in range(max_order + 1):
                P_l = legendre(l)
                moments[l] += norm_i * norm_j * P_l(cos_theta_ij)

    # Normalize moments
    normalization = total_momentum**2
    moments = [moment / normalization for moment in moments]

    return moments

# Define histograms
def define_histograms():
    histograms = {
        "histJetsSize": ROOT.TH1F("jet_Size", "Jet Size", 18, 0.0, 17.0),
        "histbJetsSize": ROOT.TH1F("bjets_Size", "B-Jets Size", 17, 0.0, 17.0),
        "histScalarHT": ROOT.TH1F("Scalar_HT", "Scalar HT", 100, 100.0, 1800.0),
        "histMET": ROOT.TH1F("MET", "MET", 100, 0.0, 600.0),
    }

    for i in range(4):
        histograms[f"bjet{i+1}_hPt"] = ROOT.TH1F(f"bjet{i+1}_hPt", f"bjet{i+1}_hPt", 100, 0, 400)
        histograms[f"bjet{i+1}_hEta"] = ROOT.TH1F(f"bjet{i+1}_hEta", f"bjet{i+1}_hEta", 100, -4.0, 4.0)
        histograms[f"bjet{i+1}_hPhi"] = ROOT.TH1F(f"bjet{i+1}_hPhi", f"bjet{i+1}_hPhi", 100, 0.0, 3.2)
        histograms[f"bjet{i+1}_hMass"] = ROOT.TH1F(f"bjet{i+1}_hMass", f"bjet{i+1}_hMass", 100, 0.0, 200.0)

    return histograms

# Define tree and branches
def define_tree():
    my_tree = ROOT.TTree("my_Tree", "Tree with multiple branches")
    

    # Define sphericity, aplanarity and circularity
    branches = {
        "sphericity_jets": array('f', [0.]),
        "aplanarity_jets": array('f', [0.]),
        "circularity_jets": array('f', [0.]),
        "sphericity_bjets": array('f', [0.]),
        "aplanarity_bjets": array('f', [0.]),
        "sphericity_leps": array('f', [0.]),
        "aplanarity_leps": array('f', [0.]),
        "jets_size":  array('i', [0]),
        "bjets_size": array('i', [0]),
        "MET_": array('f', [0.]),
        "SHT": array('f', [0.]),
        "mT": array('f', [0.])

    }

    for name, branch in branches.items():
        my_tree.Branch(name, branch, f"{name}/F")
    #branches["mT"] = array('f', [0.0])
    #my_tree.Branch("mT", branches["mT"], "mT/F")

# Define dynamic branches for b-jets and leptons
    n_bjets, n_leps = 4, 4
    branches.update({
        f"bjet{i+1}_pt_br": array('f', [0.]) for i in range(n_bjets)})
    branches.update({
        f"mT_bjet{i+1}_met_br": array('f', [0.]) for i in range(n_bjets)
    })
    for i in range(n_bjets):
        my_tree.Branch(f"bjet{i+1}_pt", branches[f"bjet{i+1}_pt_br"], f"bjet{i+1}_pt/F")
        my_tree.Branch(f"mT_bjet{i+1}_met", branches[f"mT_bjet{i+1}_met_br"], f"mT_bjet{i+1}_met/F")

# Define pairwise dPhi and dR branches for b-jets
    for i in range(n_bjets):
        for j in range(i + 1, n_bjets):
            branches[f"bjet_dphi_br{i}_{j}"] = array('f', [0.])
            my_tree.Branch(f"dPhi_bjet{i}_{j}", branches[f"bjet_dphi_br{i}_{j}"], f"dPhi_bjet{i}_{j}/F")
    for i in range(n_bjets):
        for j in range(i + 1, n_bjets):
            branches[f"bjet_dr_br{i}_{j}"] = array('f', [0.])
            my_tree.Branch(f"dR_bjet{i}_{j}", branches[f"bjet_dr_br{i}_{j}"], f"dR_bjet{i}_{j}/F")

# for leptons
    branches.update({
        f"lep{i}_pt_br": array('f', [0.]) for i in range(n_leps)
    })
    for i in range(n_leps):
        my_tree.Branch(f"lep{i}_pt", branches[f"lep{i}_pt_br"], f"lep{i}_pt/F")

# leading and sub-leading bjets and leptons branches    
    n_lead_bjets = 2
    branches[f"br_dphi_bjet_lep_leading"] = [array('f', [0.]) for _ in range(n_lead_bjets)] # dphi and dr between leading bjet and leading lepton
    branches[f"br_dr_bjet_lep_leading"] = [array('f', [0.]) for _ in range(n_lead_bjets)] 
# Create branches dynamically using a loop for bjets

    for i, label in enumerate(["leading", "subleading"]):
        my_tree.Branch(f"dphi_{label}_bjet_lep", branches[f"br_dphi_bjet_lep_leading"][i], f"dphi_{label}_bjet_lep/F")
        my_tree.Branch(f"dr_{label}_bjet_lep", branches[f"br_dr_bjet_lep_leading"][i], f"dr_{label}_bjet_lep/F")

    # Fox-Wolfram Moments
    for l in range(max_order + 1):
        branches[f"bjets_fox_wolfram_H{l}"] = array('f', [0.])
        my_tree.Branch(f"bjets_fox_wolfram_H{l}", branches[f"bjets_fox_wolfram_H{l}"], f"bjets_fox_wolfram_H{l}/F")
    # Fox-Wolfram Moments
    for l in range(max_order + 1):
        branches[f"leps_fox_wolfram_H{l}"] = array('f', [0.])
        my_tree.Branch(f"leps_fox_wolfram_H{l}", branches[f"leps_fox_wolfram_H{l}"], f"leps_fox_wolfram_H{l}/F")

    return my_tree, branches

