#!/usr/bin/env python3
"""
Analyze CMS data: cluster jets, assign labels, calculate impact parameters, and plot
"""

import numpy as np
import awkward as ak
import uproot
import matplotlib.pyplot as plt
from jet_helper import Jet, get_cluster_sequence
import fastjet as fj
from tqdm import tqdm
from multiprocessing import Pool

import argparse


def to_ak(pt, eta, phi, mass=None):
    """Convert pt, eta, phi to awkward Momentum4D array for a single event"""
    import numpy as np

    if mass is None:
        mass = np.zeros_like(pt)

    # Convert to numpy arrays if needed
    pt = np.asarray(pt)
    eta = np.asarray(eta)
    phi = np.asarray(phi)
    mass = np.asarray(mass)

    return ak.zip(
        {
            "px": pt * np.cos(phi),
            "py": pt * np.sin(phi),
            "pz": pt * np.sinh(eta),
            "E": np.sqrt(mass**2 + (pt * np.cosh(eta)) ** 2),
        },
        with_name="Momentum4D",
    )


def cluster_jets_event(pt, eta, phi, jetdef, ptmin=20):
    """Cluster jets for a single event using FastJet"""
    if len(pt) == 0:
        return [], np.array([], dtype=int)

    particles = to_ak(pt, eta, phi)
    cs = get_cluster_sequence(
        jetdef, particles, user_indices=list(range(len(particles)))
    )
    jets = cs.inclusive_jets(ptmin)
    jets = fj.sorted_by_pt(jets)
    jets = [Jet(j, 0.5, calc_substructure=True) for j in jets]
    jets = [j for j in jets if (j.nconstituents >= 2 and abs(j.eta()) < 2.5)]
    # jets = jets[:2]  # keep only leading 2 jets

    used_indices = set()

    jet_idxs = np.zeros(len(pt), dtype=int)
    for jet_idx, jet in enumerate(jets):
        particle_idx = jet.constituents_idx
        jet_idxs[particle_idx] = jet_idx
        used_indices.update(particle_idx)
    particle_idx = np.arange(len(pt))
    particle_idx = particle_idx[~np.isin(particle_idx, list(used_indices))]
    jet_idxs[particle_idx] = -1

    return jets, jet_idxs


def match_jets_to_hadrons(gen_jets, hadron_pt, hadron_eta, hadron_phi, hadron_labels, dr_threshold=0.3):
    """Match gen jets to heavy hadrons using unique dR matching
    
    Each hadron is matched to at most one jet (the closest one).
    Uses greedy matching with B-hadron precedence: B-hadrons are matched first,
    then C-hadrons, ensuring B-jets take priority.
    
    Args:
        gen_jets: List of fastjet Jet objects
        hadron_pt, hadron_eta, hadron_phi: Arrays of hadron kinematics
        hadron_labels: Array of hadron labels (5 for B-hadron, 4 for C-hadron)
        dr_threshold: Maximum dR for matching (default 0.3)
    
    Returns:
        Array of jet labels (5 for b-jet, 4 for c-jet, 0 for light-jet)
    """
    if len(gen_jets) == 0:
        return np.array([], dtype=int)
    
    # If no hadrons in event, all jets are light jets
    if len(hadron_pt) == 0:
        return np.zeros(len(gen_jets), dtype=int)
    
    jet_labels = np.zeros(len(gen_jets), dtype=int)  # Default to light-jet
    
    # Get jet kinematics
    jet_eta = np.array([j.eta() for j in gen_jets])
    jet_phi = np.array([j.phi() for j in gen_jets])
    
    # Calculate dR between each jet and each hadron
    # deta: shape (n_jets, n_hadrons)
    deta = jet_eta[:, np.newaxis] - hadron_eta[np.newaxis, :]
    dphi = jet_phi[:, np.newaxis] - hadron_phi[np.newaxis, :]
    # Wrap phi difference
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))
    dr_matrix = np.sqrt(deta**2 + dphi**2)
    
    # Create priority key: B-hadrons (label=5) first, then C-hadrons (label=4)
    # For each (jet, hadron) pair, create tuple (priority, dR, jet_idx, hadron_idx)
    # where priority = 0 for B-hadrons, 1 for C-hadrons
    matches = []
    for jet_idx in range(len(gen_jets)):
        for hadron_idx in range(len(hadron_pt)):
            dr_val = dr_matrix[jet_idx, hadron_idx]
            if dr_val <= dr_threshold:
                priority = 0 if hadron_labels[hadron_idx] == 5 else 1  # B-hadrons first
                matches.append((priority, dr_val, jet_idx, hadron_idx))
    
    # Sort by priority (B first), then by dR (closest first)
    matches.sort(key=lambda x: (x[0], x[1]))
    
    used_jets = set()
    used_hadrons = set()
    
    for priority, dr_val, jet_idx, hadron_idx in matches:
        # Check if already matched
        if jet_idx in used_jets or hadron_idx in used_hadrons:
            continue
        
        # Make the match
        jet_labels[jet_idx] = hadron_labels[hadron_idx]
        used_jets.add(jet_idx)
        used_hadrons.add(hadron_idx)
        
        # Optimization: stop if all jets or hadrons are matched
        if len(used_jets) == len(gen_jets) or len(used_hadrons) == len(hadron_pt):
            break
    
    return jet_labels


def match_jets(gen_jets, pflow_jets, dr_threshold=0.3):
    """
    Match pflow jets to gen jets uniquely (1-to-1) using greedy DeltaR minimization.
    """
    if len(gen_jets) == 0 or len(pflow_jets) == 0:
        return np.full(len(pflow_jets), -1)

    # 1. Calculate deltaR matrix for ALL pairs
    gen_eta = np.array([j.eta() for j in gen_jets])
    gen_phi = np.array([j.phi() for j in gen_jets])
    pflow_eta = np.array([j.eta() for j in pflow_jets])
    pflow_phi = np.array([j.phi() for j in pflow_jets])

    # deta: shape (n_pflow, n_gen)
    deta = pflow_eta[:, np.newaxis] - gen_eta[np.newaxis, :]
    dphi = pflow_phi[:, np.newaxis] - gen_phi[np.newaxis, :]
    # Wrap phi difference
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))
    dr_matrix = np.sqrt(deta**2 + dphi**2)

    # 2. Flatten and sort by deltaR
    # We want indices that sort the flattened array
    flat_indices = np.argsort(dr_matrix, axis=None)

    # 3. Iterate and assign unique matches
    matched_gen_indices = np.full(len(pflow_jets), -1)
    used_gen = set()
    used_pflow = set()

    for idx in flat_indices:
        # Convert flat index back to (pflow_idx, gen_idx)
        pflow_idx, gen_idx = np.unravel_index(idx, dr_matrix.shape)

        dr_val = dr_matrix[pflow_idx, gen_idx]

        # Stop if we exceeded threshold
        if dr_val > dr_threshold:
            break

        # Check if already matched
        if pflow_idx in used_pflow or gen_idx in used_gen:
            continue

        # Make the match
        matched_gen_indices[pflow_idx] = gen_idx
        used_pflow.add(pflow_idx)
        used_gen.add(gen_idx)

        # Optimization: If we matched everything possible, stop early
        if len(used_pflow) == len(pflow_jets) or len(used_gen) == len(gen_jets):
            break

    return matched_gen_indices


def plot_impact_parameters(
    b_jet_data, c_jet_data, light_jet_data, output_prefix="impact_params"
):
    """Plot impact parameter distributions for different jet flavors"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot settings
    # For dxy and ip3d: bin width = 0.04 mm = 0.004 cm
    # For significance: bin width = 0.7
    bin_width_mm = 0.02  # mm
    bin_width_cm = bin_width_mm / 10.0  # convert to cm (0.004 cm)
    bin_width_sig = 0.7  # significance bin width

    # Calculate number of bins for dxy and ip3d with 0.04 mm bin width
    dxy_range = (-0.1, 0.1)  # cm
    ip3d_range = (-0.1, 0.1)  # cm
    dxy_sig_range = (-20, 40)
    ip3d_sig_range = (-30, 30)

    n_bins_dxy = int((dxy_range[1] - dxy_range[0]) / bin_width_cm)
    n_bins_ip3d = int((ip3d_range[1] - ip3d_range[0]) / bin_width_cm)
    n_bins_dxy_sig = int((dxy_sig_range[1] - dxy_sig_range[0]) / bin_width_sig)
    n_bins_ip3d_sig = int((ip3d_sig_range[1] - ip3d_sig_range[0]) / bin_width_sig)

    bins_dxy = np.linspace(dxy_range[0], dxy_range[1], n_bins_dxy + 1)
    bins_dxy_sig = np.linspace(dxy_sig_range[0], dxy_sig_range[1], n_bins_dxy_sig + 1)
    bins_ip3d = np.linspace(ip3d_range[0], ip3d_range[1], n_bins_ip3d + 1)
    bins_ip3d_sig = np.linspace(
        ip3d_sig_range[0], ip3d_sig_range[1], n_bins_ip3d_sig + 1
    )

    variables = ["dxy", "dxy_sig", "ip3d", "ip3d_sig"]
    bins_dict = {
        "dxy": bins_dxy,
        "dxy_sig": bins_dxy_sig,
        "ip3d": bins_ip3d,
        "ip3d_sig": bins_ip3d_sig,
    }
    axes = axes.flatten()
    x_label_dict = {
        "dxy": "dxy (cm)",
        "dxy_sig": "dxy Significance",
        "ip3d": "3D Impact Parameter (cm)",
        "ip3d_sig": "3D Impact Parameter Significance",
    }
    y_label_dict = {
        "dxy": f"Fraction of tracks / {bin_width_mm:.2f} mm",
        "dxy_sig": f"Tracks / {bin_width_sig:.1f}",
        "ip3d": f"Fraction of tracks / {bin_width_mm:.2f} mm",
        "ip3d_sig": f"Tracks / {bin_width_sig:.1f}",
    }

    for i in range(4):
        var = variables[i]
        # Calculate weights to get fraction of tracks per bin
        b_weights = (
            np.ones_like(b_jet_data[var]) / len(b_jet_data[var])
            if len(b_jet_data[var]) > 0
            else None
        )
        c_weights = (
            np.ones_like(c_jet_data[var]) / len(c_jet_data[var])
            if len(c_jet_data[var]) > 0
            else None
        )
        light_weights = (
            np.ones_like(light_jet_data[var]) / len(light_jet_data[var])
            if len(light_jet_data[var]) > 0
            else None
        )

        axes[i].hist(
            b_jet_data[var],
            bins=bins_dict[var],
            weights=b_weights,
            alpha=0.5,
            label="b-jets",
            color="red",
            histtype="step",
        )
        axes[i].hist(
            c_jet_data[var],
            bins=bins_dict[var],
            weights=c_weights,
            alpha=0.5,
            label="c-jets",
            color="green",
            histtype="step",
        )
        axes[i].hist(
            light_jet_data[var],
            bins=bins_dict[var],
            weights=light_weights,
            alpha=0.5,
            label="light-jets",
            color="blue",
            histtype="step",
        )

        axes[i].set_xlabel(x_label_dict[var])
        axes[i].set_ylabel(y_label_dict[var])
        axes[i].set_title(f"{x_label_dict[var]} Distribution")
        axes[i].legend()
        axes[i].set_yscale("log")
        axes[i].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches="tight")
    print(f"Saved plot to {output_prefix}.png")
    plt.close()


def get_3d_pca_vector(d0, z0, track_px, track_py, track_pz, jet_px, jet_py):
    """
    Approximates the 3D vector from the Primary Vertex to the track's Point of
    Closest Approach (PCA).

    NOTE: This is a simplified approximation for pseudocode demonstration.
    A real implementation requires a full helical track model. In frameworks
    like CMSSW, this is handled by dedicated tools like IPTools.

    The main simplification here is assuming the z-coordinate of the 3D PCA
    is simply z0 and that the xy-projection of the PCA vector can be
    determined from the 2D geometry alone.

    Args:
        d0 (float): Magnitude of the transverse impact parameter.
        z0 (float): Magnitude of the longitudinal impact parameter.
        track_px (float): Particle's momentum x-component.
        track_py (float): Particle's momentum y-component.
        track_pz (float): Particle's momentum z-component.
        jet_px (float): Jet's momentum x-component.
        jet_py (float): Jet's momentum y-component.

    Returns:
        np.ndarray: The approximated 3D vector to the PCA.
    """
    # The vector to the PCA in the transverse plane is perpendicular to the
    # track's transverse momentum.
    track_pt = np.sqrt(track_px**2 + track_py**2)
    if track_pt == 0:
        return np.array([0, 0, z0])  # Cannot determine direction

    # Unit vector perpendicular to transverse track momentum
    # This corresponds to the standard definition where d0 = (r_pca x p_track)_z / p_t
    perp_vec_xy_unit = np.array([-track_py / track_pt, track_px / track_pt])

    # The PCA vector in XY is simply d0 * perp_direction
    # We do NOT use jet direction here; the vector is a property of the track/PV only.
    vec_pca_xy = d0 * perp_vec_xy_unit

    # Combine with z0 to form a simplified 3D PCA vector.
    # This is the primary approximation.
    vec_pca_3d = np.array([vec_pca_xy[0], vec_pca_xy[1], z0])

    return vec_pca_3d


def calculate_signed_impact_parameter(
    d0,
    z0,
    track_px,
    track_py,
    track_pz,
    jet_px,
    jet_py,
    jet_pz,
    d0_error=None,
    z0_error=None,  # Optional for significance
    track_vx=None,
    track_vy=None,

):
    """
    Calculates the signed 2D (dxy) and 3D impact parameters and their significance.

    The sign is determined by the dot product of the jet direction and the
    vector from the primary vertex (PV) to the track's point of closest
    approach (PCA). A positive sign indicates the track originated from a
    decay displaced in the direction of the jet.

    Args:
        d0 (float): Transverse impact parameter magnitude (distance to PV in xy-plane).
        z0 (float): Longitudinal impact parameter (z-distance to PV at point of closest xy-approach).
        track_px, track_py, track_pz (float): Particle's momentum components.
        jet_px, jet_py, jet_pz (float): Jet's momentum components.
        d0_error (float, optional): Uncertainty on d0.
        z0_error (float, optional): Uncertainty on z0.

    Returns:
        dict: A dictionary containing the signed parameters and their significance.
    """
    # --- 1. Signed Transverse Impact Parameter (d_xy) ---
    # This calculation is geometrically straightforward in the 2D plane.

    track_phi = np.arctan2(track_py, track_px)
    track_theta = np.arctan2(np.sqrt(track_px**2 + track_py**2), track_pz)
    jet_phi = np.arctan2(jet_py, jet_px)

    # # The signed d_xy is given by d0 * sin(delta_phi), where delta_phi is the
    # # angle between the jet and the track. This formula elegantly combines
    # # the magnitude (d0) with the sign from the relative orientation.
    # dxy_signed = d0 * np.sin(jet_phi - track_phi)

    # # # --- 2. Signed 3D Impact Parameter (IP3D) ---
    # # # The principle is the same, but in 3D. We need the dot product of the
    # # # 3D jet vector and the 3D PCA vector.

    # # # Construct 3D vectors
    # # vec_jet = np.array([jet_px, jet_py, jet_pz]).reshape(-1)

    # # # Get the (approximated) 3D vector to the PCA
    # # vec_pca_3d = get_3d_pca_vector(
    # #     d0, z0, track_px, track_py, track_pz, jet_px, jet_py
    # # ).reshape(-1)
    # # # Calculate the dot product to determine the sign
    # # dot_product_3d = np.dot(vec_pca_3d, vec_jet)

    # sign_2d = np.sign(dxy_signed)

    # # Combine the magnitude and the sign
    # ip3d_magnitude = np.sqrt(d0**2 + z0**2)
    # # ip3d_signed = np.sign(dot_product_3d) * ip3d_magnitude
    # ip3d_signed = sign_2d * ip3d_magnitude


    # Calculate the "Physics Sign"
    # This check determines if the track crosses the jet axis downstream (+) or upstream (-)
    # We use d0 * sin(...) to get the correct geometric crossing side.
    crossing_geometry = d0 * np.sin(jet_phi - track_phi)
    phys_sign = np.sign(crossing_geometry)

    # Delphes formula
    # phys_sign = np.sign(jet_px * track_vx + jet_py * track_vy) 

    # The Signed IP is the MAGNITUDE (abs) times the SIGN
    dxy_signed = np.abs(d0) * phys_sign

    # --- 2. Signed 3D Impact Parameter (IP3D) ---
    # Use the same 2D sign for the 3D parameter
    # Version with sin(theta) factor:
    # ip3d_magnitude = np.sqrt(d0**2 + (z0 * np.sin(track_theta))**2)
    # Original version without sin(theta):
    ip3d_magnitude = np.sqrt(d0**2 + z0**2)
    ip3d_signed = ip3d_magnitude * phys_sign

    # --- 3. Significance (optional) ---
    # Significance is a powerful metric, defined as value / error.
    dxy_significance = None
    ip3d_significance = None

    if d0_error is not None:
        # Assuming dxy error is approximately d0_error
        dxy_significance = dxy_signed / d0_error if d0_error > 0 else np.array([0])

    if d0_error is not None and z0_error is not None:
        # Propagate errors for 3D IP. Assuming d0 and z0 are uncorrelated.
        # ip3d_error = sqrt( (d0/ip3d * d0_err)^2 + (z0/ip3d * z0_err)^2 )
        if ip3d_magnitude > 0:
            ip3d_error = (
                np.sqrt((d0 * d0_error) ** 2 + (z0 * z0_error) ** 2) / ip3d_magnitude
            )
            ip3d_significance = (
                ip3d_signed / ip3d_error if ip3d_error > 0 else np.array([0])
            )
        else:
            ip3d_significance = np.array([0])

    return {
        "dxy_signed": dxy_signed,
        "ip3d_signed": ip3d_signed,
        "dxy_significance": dxy_significance,
        "ip3d_significance": ip3d_significance,
    }


def process_event(event_data):
    """Process a single event: cluster jets, match, and calculate impact parameters

    Args:
        event_data: tuple containing (event_idx, gen_event_data, pfc_event_data, hadron_event_data)
    Returns:
        dict with results for this event
    """
    event_idx, gen_event, pfc_event, hadron_event = event_data

    # Define anti-kT jet algorithm with R=0.5
    jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.5)

    # Extract gen particle data for this event
    gen_pt = np.asarray(gen_event["pt"])
    gen_eta = np.asarray(gen_event["eta"])
    gen_phi = np.asarray(gen_event["phi"])
    gen_labels = np.asarray(gen_event["labels"])

    # Extract hadron data for this event
    had_pt = np.asarray(hadron_event["pt"])
    had_eta = np.asarray(hadron_event["eta"])
    had_phi = np.asarray(hadron_event["phi"])
    had_labels = np.asarray(hadron_event["labels"])

    # Extract PF candidate data for this event
    pfc_pt = np.asarray(pfc_event["pt"])
    pfc_eta = np.asarray(pfc_event["eta"])
    pfc_phi = np.asarray(pfc_event["phi"])
    pfc_d0 = np.asarray(pfc_event["d0"])
    pfc_d0_err = np.asarray(pfc_event["d0_err"])
    pfc_z0 = np.asarray(pfc_event["z0"])
    pfc_z0_err = np.asarray(pfc_event["z0_err"])

    pfc_vx = np.asarray(pfc_event["vx"])
    pfc_vy = np.asarray(pfc_event["vy"])
    pfc_vz = np.asarray(pfc_event["vz"])

    # Cluster gen jets for this event
    gen_jets, gen_jet_idxs = cluster_jets_event(
        gen_pt, gen_eta, gen_phi, jetdef, ptmin=20
    )

    # Assign labels to gen jets using dR matching to hadrons
    gen_jet_labels = match_jets_to_hadrons(
        gen_jets, had_pt, had_eta, had_phi, had_labels, dr_threshold=0.3
    )

    # Cluster pflow jets for this event
    pflow_jets, pflow_jet_idxs = cluster_jets_event(
        pfc_pt, pfc_eta, pfc_phi, jetdef, ptmin=20
    )

    # Match pflow jets to gen jets
    matched_gen_idx = match_jets(gen_jets, pflow_jets, dr_threshold=0.3)

    # Assign labels to pflow jets
    pflow_jet_labels = np.full(len(pflow_jets), -1)
    for i, gen_idx in enumerate(matched_gen_idx):
        if gen_idx >= 0:
            pflow_jet_labels[i] = gen_jet_labels[gen_idx]

    # Storage for this event's impact parameters
    event_results = {
        "b_jet": {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []},
        "c_jet": {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []},
        "light_jet": {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []},
        "n_gen_jets": len(gen_jets),
        "n_pflow_jets": len(pflow_jets),
        "gen_jet_labels": gen_jet_labels,
        "pflow_jet_labels": pflow_jet_labels,
    }

    # Calculate impact parameters for each pflow jet
    for jet_idx, (jet, label) in enumerate(zip(pflow_jets, pflow_jet_labels)):
        if label == -1:
            continue

        # Get constituents of this jet
        constituent_mask = pflow_jet_idxs == jet_idx

        # Extract impact parameter data and kinematics for constituents
        d0 = pfc_d0[constituent_mask]
        d0_err = pfc_d0_err[constituent_mask]
        z0 = pfc_z0[constituent_mask]
        z0_err = pfc_z0_err[constituent_mask]
        pt = pfc_pt[constituent_mask]
        eta = pfc_eta[constituent_mask]
        phi = pfc_phi[constituent_mask]

        vx = pfc_vx[constituent_mask]
        vy = pfc_vy[constituent_mask]
        vz = pfc_vz[constituent_mask]

        # Filter for charged particles only (neutral particles have d0 = -1000)
        charged_mask = d0 > -999
        if not np.any(charged_mask):
            continue  # Skip jets with no charged constituents

        d0 = d0[charged_mask]
        d0_err = d0_err[charged_mask]
        z0 = z0[charged_mask]
        z0_err = z0_err[charged_mask]
        pt = pt[charged_mask]
        eta = eta[charged_mask]
        phi = phi[charged_mask]

        vx = vx[charged_mask]
        vy = vy[charged_mask]
        vz = vz[charged_mask]

        # Calculate track momentum components
        track_px = pt * np.cos(phi)
        track_py = pt * np.sin(phi)
        track_pz = pt * np.sinh(eta)

        # Get jet direction (px, py, pz)
        jet_px = jet.px()
        jet_py = jet.py()
        jet_pz = jet.pz()

        # Calculate signed impact parameters for each constituent
        dxy = []
        dxy_sig = []
        ip3d_signed = []
        ip3d_sig_signed = []

        for i in range(len(d0)):
            result = calculate_signed_impact_parameter(
                d0[i],
                z0[i],
                track_px[i],
                track_py[i],
                track_pz[i],
                jet_px,
                jet_py,
                jet_pz,
                d0_error=d0_err[i],
                z0_error=z0_err[i],
                track_vx=vx[i],
                track_vy=vy[i],
            )
            dxy.append(result["dxy_signed"])
            dxy_sig.append(
                result["dxy_significance"]
                if result["dxy_significance"] is not None
                else 0
            )
            ip3d_signed.append(result["ip3d_signed"])
            ip3d_sig_signed.append(
                result["ip3d_significance"]
                if result["ip3d_significance"] is not None
                else 0
            )

        dxy = np.array(dxy)
        dxy_sig = np.array(dxy_sig)
        ip3d_signed = np.array(ip3d_signed)
        ip3d_sig_signed = np.array(ip3d_sig_signed)

        # Store in appropriate category
        if label == 5:  # b-jet
            event_results["b_jet"]["dxy"].extend(dxy)
            event_results["b_jet"]["dxy_sig"].extend(dxy_sig)
            event_results["b_jet"]["ip3d"].extend(ip3d_signed)
            event_results["b_jet"]["ip3d_sig"].extend(ip3d_sig_signed)
        elif label == 4:  # c-jet
            event_results["c_jet"]["dxy"].extend(dxy)
            event_results["c_jet"]["dxy_sig"].extend(dxy_sig)
            event_results["c_jet"]["ip3d"].extend(ip3d_signed)
            event_results["c_jet"]["ip3d_sig"].extend(ip3d_sig_signed)
        elif label == 0:  # light-jet
            event_results["light_jet"]["dxy"].extend(dxy)
            event_results["light_jet"]["dxy_sig"].extend(dxy_sig)
            event_results["light_jet"]["ip3d"].extend(ip3d_signed)
            event_results["light_jet"]["ip3d_sig"].extend(ip3d_sig_signed)

    return event_results


def main(input_file, num_events=1000, n_jobs=-1):
    """Main analysis function"""

    print(f"Loading data from {input_file}...")

    # Load data from ROOT file
    file = uproot.open(input_file)
    gen_tree = file["gens/Events"]
    hadron_tree = file["hadrons/Events"]
    pfc_tree = file["pfcs/Events"]
    if num_events > gen_tree.num_entries:
        num_events = gen_tree.num_entries
        print(
            f"Warning: num_events exceeds available entries. Using {num_events} events."
        )
    if num_events < 1:
        num_events = gen_tree.num_entries
        print("Using all available events -- num_events set to", num_events)
    # Load gen particle data as awkward arrays (jagged)
    gen_data = gen_tree.arrays(
        ["GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_label"],
        library="ak",
        entry_stop=num_events,
    )
    hadron_data = hadron_tree.arrays(
        [
            "HeavyHadron_pt",
            "HeavyHadron_eta",
            "HeavyHadron_phi",
            "HeavyHadron_mass",
            "HeavyHadron_label",
        ],
        library="ak",
        entry_stop=num_events,
    )

    # Load PF candidate data as awkward arrays (jagged)
    pfc_data = pfc_tree.arrays(
        [
            "PFCand_pt",
            "PFCand_eta",
            "PFCand_phi",
            "PFCand_mass",
            "PFCand_d0",
            "PFCand_d0Error",
            "PFCand_z0",
            "PFCand_z0Error",
            "PFCand_ip3d",
            "PFCand_ip3dError",
            "PFCand_vx",
            "PFCand_vy",
            "PFCand_vz",
        ],
        library="ak",
        entry_stop=num_events,
    )

    n_events = len(gen_data["GenPart_pt"])
    print(f"Loaded {n_events} events")

    # Determine number of parallel workers
    if n_jobs == -1:
        import multiprocessing

        n_jobs = multiprocessing.cpu_count()
    elif n_jobs <= 0:
        n_jobs = 1

    print(f"Using {n_jobs} parallel workers")

    # Impact parameter data by jet category
    b_jet_data = {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []}
    c_jet_data = {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []}
    light_jet_data = {"dxy": [], "dxy_sig": [], "ip3d": [], "ip3d_sig": []}

    # Storage for statistics
    all_gen_jet_labels = []
    all_pflow_jet_labels = []

    # Prepare data for multiprocessing
    print("\nPreparing data for parallel processing...")
    event_data_list = []
    for event_idx in range(n_events):
        gen_event = {
            "pt": gen_data["GenPart_pt"][event_idx],
            "eta": gen_data["GenPart_eta"][event_idx],
            "phi": gen_data["GenPart_phi"][event_idx],
            "labels": gen_data["GenPart_label"][event_idx],
        }
        pfc_event = {
            "pt": pfc_data["PFCand_pt"][event_idx],
            "eta": pfc_data["PFCand_eta"][event_idx],
            "phi": pfc_data["PFCand_phi"][event_idx],
            "d0": pfc_data["PFCand_d0"][event_idx],
            "d0_err": pfc_data["PFCand_d0Error"][event_idx],
            "z0": pfc_data["PFCand_z0"][event_idx],
            "z0_err": pfc_data["PFCand_z0Error"][event_idx],
            "ip3d": pfc_data["PFCand_ip3d"][event_idx],
            "ip3d_err": pfc_data["PFCand_ip3dError"][event_idx],
            "vx": pfc_data["PFCand_vx"][event_idx],
            "vy": pfc_data["PFCand_vy"][event_idx],
            "vz": pfc_data["PFCand_vz"][event_idx],
        }
        hadron_event = {
            "pt": hadron_data["HeavyHadron_pt"][event_idx],
            "eta": hadron_data["HeavyHadron_eta"][event_idx],
            "phi": hadron_data["HeavyHadron_phi"][event_idx],
            "mass": hadron_data["HeavyHadron_mass"][event_idx],
            "labels": hadron_data["HeavyHadron_label"][event_idx],
        }
        event_data_list.append((event_idx, gen_event, pfc_event, hadron_event))

    # Process events in parallel
    print(f"\nProcessing {n_events} events with {n_jobs} workers...")
    if n_jobs == 1:
        # Single-threaded for debugging
        results = []
        for event_data in tqdm(event_data_list, desc="Processing events"):
            results.append(process_event(event_data))
    else:
        # Multi-threaded processing
        with Pool(n_jobs) as pool:
            results = list(
                tqdm(
                    pool.imap(process_event, event_data_list),
                    total=n_events,
                    desc="Processing events",
                )
            )

    # Aggregate results from all events
    print("\nAggregating results...")
    for result in results:
        # Accumulate impact parameters by jet flavor
        for key in ["dxy", "dxy_sig", "ip3d", "ip3d_sig"]:
            b_jet_data[key].extend(result["b_jet"][key])
            c_jet_data[key].extend(result["c_jet"][key])
            light_jet_data[key].extend(result["light_jet"][key])

        # Collect labels for statistics
        all_gen_jet_labels.extend(result["gen_jet_labels"])
        all_pflow_jet_labels.extend(result["pflow_jet_labels"])

    # Convert to numpy arrays
    for data in [b_jet_data, c_jet_data, light_jet_data]:
        for key in data:
            data[key] = np.array(data[key])

    all_gen_jet_labels = np.array(all_gen_jet_labels)
    all_pflow_jet_labels = np.array(all_pflow_jet_labels)

    # Print statistics
    print("\nTotal jets found:")
    n_gen = len(all_gen_jet_labels)
    n_b = np.sum(all_gen_jet_labels == 5)
    n_c = np.sum(all_gen_jet_labels == 4)
    n_light = np.sum(all_gen_jet_labels == 0)
    print(f"  Gen jets: {n_gen}")
    print(f"    b-jets: {n_b}, c-jets: {n_c}, light-jets: {n_light}")

    n_pflow = len(all_pflow_jet_labels)
    n_matched = np.sum(all_pflow_jet_labels >= 0)
    n_b_pf = np.sum(all_pflow_jet_labels == 5)
    n_c_pf = np.sum(all_pflow_jet_labels == 4)
    n_light_pf = np.sum(all_pflow_jet_labels == 0)
    print(f"  Pflow jets: {n_pflow}")
    print(
        f"    Matched: {n_matched}, b-jets: {n_b_pf}, c-jets: {n_c_pf}, light-jets: {n_light_pf}"
    )

    print("\nImpact parameters calculated:")
    print(f"  b-jets: {len(b_jet_data['dxy'])} constituents")
    print(f"  c-jets: {len(c_jet_data['dxy'])} constituents")
    print(f"  light-jets: {len(light_jet_data['dxy'])} constituents")

    # Plot results
    print("\nGenerating plots...")
    plot_impact_parameters(b_jet_data, c_jet_data, light_jet_data)

    print("\nAnalysis complete!")


parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="Input ROOT file with CMS data")
parser.add_argument(
    "-n", "--num_events", type=int, default=1000, help="Number of events to process"
)
parser.add_argument(
    "-j",
    "--n_jobs",
    type=int,
    default=-1,
    help="Number of parallel workers (-1 for all CPUs, 1 for sequential)",
)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.input_file, args.num_events, args.n_jobs)
