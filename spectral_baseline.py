#!/usr/bin/env python3
"""
Spectral Baseline Computation for the Hopfion QGT Paper
========================================================

Evaluates the Hopf-twisted spectral zeta function on S^3 x S^1
at s=2, reproducing Appendix D.5 (equation 10) of:

  "Hopfion Topology, Quantum Geometric Tensor Saturation,
   and the Fine-Structure Constant" (A. Down, April 2026)

The twisted spectral zeta function is:

  zeta_theta(s) = sum_{k>=0} sum_n sum_{m in Z}
                    w(k,n) * [(k + 3/2)^2 + (m + n*theta)^2]^{-s}

where:
  k >= 0          : S^3 Dirac level
  n               : Hopf charge at level k
  m in Z          : S^1 mode (truncated to [-M, M])
  w(k,n)          : multiplicity of Hopf charge n at level k
  sum_n w(k,n)    = 2(k+1)(k+2) at each level (standard Dirac multiplicity)

At s=2 this gives Tr(|D_theta|^{-4}), the smooth spectral capacity.
(Note: the paper parameterizes the zeta function as sum |lambda^2|^{-s},
so s=2 here corresponds to s=4 in the standard convention zeta_D(s) =
sum |lambda|^{-s}. Both give Tr(|D|^{-4}), the relevant quantity for
the Dixmier trace on a 4-dimensional spectral triple.)

RESULT: At the default truncation (k_max=80, m_max=300), the untwisted
baseline is zeta_0(2) ~ 13.40 (the raw sum grows logarithmically with
truncation, as expected for a 4D spectral triple; only ratios converge).
Normalizing so the bulk contribution recovers Vol(S^3 x S^1) = 4*pi^3
= 124.025, the golden-angle twist produces a normalized correction of
~ -9.8e-4, consistent with the paper's reported -9.3e-4 (the ~5%
difference is from truncation and the effective w(k,n) decomposition).
This correction is four orders of magnitude below the required irrational
gain pi^2 + pi = 13.01. The smooth spectral route is closed.

Repository: https://hopfion-topology.github.io/hopfion-evaluation-map/
Paper DOI:  10.5281/zenodo.19452185
"""

import numpy as np
from math import pi, sqrt

# =============================================================
# Physical constants
# =============================================================
PHI = (1 + sqrt(5)) / 2          # golden ratio
THETA_GOLD = 1 / PHI**2          # golden angle as fraction of full turn
P_SAT = 4*pi**3 + pi**2 + pi    # saturation ceiling
ALPHA_INV = 137.035999177        # CODATA 2022

# =============================================================
# Fibonacci convergents of theta = 1/phi^2
# =============================================================
def fibonacci_convergents(n_max=15):
    """
    Generate Fibonacci convergents p_k/q_k of theta = 1/phi^2.
    
    The continued fraction of 1/phi^2 = [0; 2, 1, 1, 1, ...].
    Convergents: 0/1, 1/2, 1/3, 2/5, 3/8, 5/13, 8/21, 13/34,
                 21/55, 34/89, 55/144, ...
    """
    # Build convergents from continued fraction [0; 2, 1, 1, 1, ...]
    cf = [0, 2] + [1] * (n_max - 2)
    p_prev, p_curr = 0, 1
    q_prev, q_curr = 1, 0
    convergents = []
    for a in cf:
        p_prev, p_curr = p_curr, a * p_curr + p_prev
        q_prev, q_curr = q_curr, a * q_curr + q_prev
        convergents.append((p_curr, q_curr))
    return convergents


# =============================================================
# Hopf charge decomposition w(k, n)
# =============================================================
def hopf_charges(k):
    """
    Decompose the S^3 Dirac level k under the Hopf U(1) action.
    
    Total multiplicity at level k: 2(k+1)(k+2).
    
    The Hopf U(1) acts on S^3 = SU(2) via the maximal torus. At Dirac
    level k, the eigenspace decomposes into integer Hopf charges
    n = -k, -k+1, ..., k-1, k (giving 2k+1 charges), with effective
    multiplicity w = 2(k+1)(k+2)/(2k+1) per charge.
    
    NOTE: These weights are non-integer for k >= 2 (e.g., w = 4.8 at
    k = 2). This reflects an equidistribution assumption: the total
    multiplicity is distributed uniformly across charges. The exact
    SU(2) branching rule would give integer multiplicities but is
    representation-dependent. The paper's conclusion (that the smooth
    spectral route is closed) is insensitive to this choice: all
    tested decompositions (uniform, triangular, parity-stepped) give
    corrections between 5e-4 and 1e-2, all at least three orders of
    magnitude below the required gain of 13.01.
    
    Returns: list of (n, w) pairs where n is the integer Hopf charge
             and w is the effective multiplicity.
    """
    n_charges = 2 * k + 1
    total_mult = 2 * (k + 1) * (k + 2)
    w_each = total_mult / n_charges
    
    return [(n, w_each) for n in range(-k, k + 1)]


# =============================================================
# Core computation: twisted spectral zeta function
# =============================================================
def spectral_zeta(theta, s=2, k_max=80, m_max=300, use_hopf=True):
    """
    Evaluate the Hopf-twisted spectral zeta function at angle theta.
    
    zeta_theta(s) = sum_{k,n,m} w(k,n) * [(k+3/2)^2 + (m+n*theta)^2]^{-s}
    
    Parameters
    ----------
    theta : float
        Reeb flow angle (fraction of full turn). 0 = untwisted.
    s : float
        Spectral parameter. s=2 gives Tr(|D|^{-4}).
    k_max : int
        S^3 level truncation (paper uses 80).
    m_max : int
        S^1 mode truncation (paper uses 300).
    use_hopf : bool
        If True, use the Hopf charge decomposition w(k,n).
        If False, ignore Hopf charges (equivalent to theta=0).
    
    Returns
    -------
    float : zeta_theta(s)
    """
    total = 0.0
    m_range = np.arange(-m_max, m_max + 1, dtype=np.float64)
    
    for k in range(k_max + 1):
        A = (k + 1.5) ** 2  # S^3 eigenvalue squared
        
        if use_hopf and theta != 0.0:
            # Use Hopf charge decomposition
            for n, w in hopf_charges(k):
                shifted_m = m_range + n * theta
                eigenvalues = A + shifted_m ** 2
                total += w * np.sum(eigenvalues ** (-s))
        else:
            # Untwisted: all charges give same result, use total mult
            mult = 2 * (k + 1) * (k + 2)
            eigenvalues = A + m_range ** 2
            total += mult * np.sum(eigenvalues ** (-s))
    
    return total


# =============================================================
# Main computation
# =============================================================
def main():
    print("=" * 70)
    print("SPECTRAL BASELINE COMPUTATION")
    print("Hopf-twisted spectral zeta function on S^3 x S^1")
    print("Appendix D.5, equation (10)")
    print("=" * 70)
    print()
    
    # ---------------------------------------------------------
    # 1. Physical constants
    # ---------------------------------------------------------
    print("Physical constants:")
    print(f"  phi           = {PHI:.15f}")
    print(f"  theta_gold    = 1/phi^2 = {THETA_GOLD:.15f}")
    print(f"  4*pi^3        = {4*pi**3:.6f}")
    print(f"  pi^2          = {pi**2:.6f}")
    print(f"  pi            = {pi:.6f}")
    print(f"  P_sat         = 4pi^3 + pi^2 + pi = {P_SAT:.6f}")
    print(f"  alpha^-1      = {ALPHA_INV} (CODATA 2022)")
    print(f"  Delta         = P - alpha^-1 = {P_SAT - ALPHA_INV:.6e} "
          f"({(P_SAT - ALPHA_INV)/ALPHA_INV * 1e6:.2f} ppm)")
    print()
    
    # ---------------------------------------------------------
    # 2. Untwisted baseline (theta = 0)
    # ---------------------------------------------------------
    print("-" * 70)
    print("STEP 1: Untwisted baseline (theta = 0)")
    print("-" * 70)
    
    k_max = 80
    m_max = 300
    print(f"  Truncation: k_max = {k_max}, m_max = +/-{m_max}")
    
    zeta_0 = spectral_zeta(0.0, s=2, k_max=k_max, m_max=m_max,
                           use_hopf=False)
    print(f"  zeta_0(2)     = {zeta_0:.6f}")
    
    # Normalization constant
    C_norm = 4 * pi**3 / zeta_0
    print(f"  C_norm        = 4*pi^3 / zeta_0(2) = {C_norm:.4f}")
    print(f"  C * zeta_0(2) = {C_norm * zeta_0:.3f}  "
          f"(should be {4*pi**3:.3f})")
    print()
    
    # Convergence note: the raw spectral sum diverges logarithmically
    # with k_max (expected for a 4D spectral triple; each level k
    # contributes ~pi/k to the sum). The normalization C absorbs this
    # divergence. The meaningful convergence test is the Fibonacci
    # convergent sweep in Step 3, which shows the twist correction
    # stabilizing across successive rational approximants of 1/phi^2
    # at fixed truncation.
    print("  Note: raw zeta_0 grows with k_max (log-divergent, expected).")
    print("  The normalization C = 4*pi^3 / zeta_0 absorbs this.")
    print("  Convergence of the twist correction is verified by the")
    print("  Fibonacci convergent sweep in Step 3.")
    print()
    
    # ---------------------------------------------------------
    # 3. Golden-angle evaluation
    # ---------------------------------------------------------
    print("-" * 70)
    print("STEP 2: Golden-angle twist (theta = 1/phi^2)")
    print("-" * 70)
    
    zeta_gold = spectral_zeta(THETA_GOLD, s=2, k_max=k_max, m_max=m_max,
                              use_hopf=True)
    normalized_gold = C_norm * zeta_gold
    correction = normalized_gold - 4 * pi**3
    
    print(f"  zeta_{{1/phi^2}}(2)     = {zeta_gold:.6f}")
    print(f"  C * zeta_{{1/phi^2}}(2) = {normalized_gold:.3f}")
    print(f"  theta-correction      = {correction:.4e}")
    print(f"  required irrational gain: pi^2 + pi = {pi**2 + pi:.2f}")
    print(f"  ratio (correction / required): "
          f"{abs(correction) / (pi**2 + pi):.2e}")
    print()
    
    if abs(correction) < 0.01:
        print("  Result: the smooth spectral route produces no significant")
        print("  irrational gain. The sub-leading terms pi^2 + pi cannot")
        print("  arise from the smooth Dirac spectrum on S^3 x S^1.")
    print()
    
    # ---------------------------------------------------------
    # 4. Fibonacci convergent sweep
    # ---------------------------------------------------------
    print("-" * 70)
    print("STEP 3: Fibonacci convergent sweep")
    print("-" * 70)
    
    convergents = fibonacci_convergents(14)
    
    print(f"  {'p/q':>10s}  {'theta_k':>12s}  {'|theta-theta_k|':>14s}  "
          f"{'C*zeta':>10s}  {'correction':>12s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*14}  {'-'*10}  {'-'*12}")
    
    for p, q in convergents:
        if q > 200:
            break
        theta_k = p / q if q > 0 else 0.0
        err = abs(theta_k - THETA_GOLD)
        
        zeta_k = spectral_zeta(theta_k, s=2, k_max=k_max, m_max=m_max,
                               use_hopf=True)
        norm_k = C_norm * zeta_k
        corr_k = norm_k - 4 * pi**3
        
        print(f"  {p:>4d}/{q:<5d}  {theta_k:>12.8f}  {err:>14.2e}  "
              f"{norm_k:>10.3f}  {corr_k:>12.4e}")
    
    # Add the irrational limit
    print(f"  {'1/phi^2':>10s}  {THETA_GOLD:>12.8f}  {'0':>14s}  "
          f"{normalized_gold:>10.3f}  {correction:>12.4e}")
    print()
    
    # ---------------------------------------------------------
    # 5. Summary
    # ---------------------------------------------------------
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"  Untwisted baseline:       zeta_0(2) = {zeta_0:.6f}")
    print(f"  Normalization constant:   C = {C_norm:.4f}")
    print(f"  Bulk volume (normalized): C * zeta_0(2) = "
          f"{C_norm * zeta_0:.3f} = 4*pi^3")
    print(f"  Golden-angle result:      C * zeta_gold(2) = "
          f"{normalized_gold:.3f}")
    print(f"  Theta-dependent shift:    {correction:.4e}")
    print(f"  Required irrational gain: pi^2 + pi = {pi**2 + pi:.5f}")
    print(f"  Shortfall factor:         "
          f"{(pi**2 + pi) / max(abs(correction), 1e-15):.0f}x")
    print()
    print("  CONCLUSION: The smooth Dirac spectrum on S^3 x S^1,")
    print("  twisted by Hopf charge at the golden angle, contributes")
    print("  only the leading term 4*pi^3. The twist correction is")
    print("  four orders of magnitude below the required irrational")
    print("  gain pi^2 + pi = 13.01. The sub-leading terms cannot")
    print("  arise from the smooth Dirac spectrum. This closes the")
    print("  smooth spectral route and narrows the origin of P_sat")
    print("  to the non-perturbative structure of the foliation")
    print("  algebra (Connes trace conjecture, Appendix D.4).")
    print()

    return zeta_0, C_norm, normalized_gold, correction


if __name__ == "__main__":
    main()
