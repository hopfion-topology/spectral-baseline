# Spectral Baseline Computation

Companion code for
*Hopfion Topology, Quantum Geometric Tensor Saturation, and the Fine-Structure Constant*
(A. Down, April 2026; [DOI: 10.5281/zenodo.19452185](https://doi.org/10.5281/zenodo.19452185)).

## Overview

This script evaluates the Hopf-twisted spectral zeta function on S³ × S¹ at s = 2, reproducing the numerical spectral baseline of Appendix D.5 (equation 10).

The twisted spectral zeta function is:

```
ζ_θ(s) = Σ_{k≥0} Σ_n Σ_{m∈Z} w(k,n) · [(k + 3/2)² + (m + nθ)²]^{-s}
```

where `w(k,n)` is the multiplicity of Hopf charge `n` at S³ Dirac level `k`, satisfying `Σ_n w(k,n) = 2(k+1)(k+2)`.

## Result

The untwisted baseline ζ₀(2) is normalized so that the bulk contribution recovers Vol(S³ × S¹) = 4π³ ≈ 124.025. The golden-angle twist (θ = 1/φ²) produces a normalized correction of approximately -9.8 × 10⁻⁴, consistent with the paper's reported -9.3 × 10⁻⁴. This correction is four orders of magnitude below the required irrational gain π² + π ≈ 13.01.

The sub-leading terms π² + π cannot arise from the smooth Dirac spectrum on S³ × S¹. This closes the smooth spectral route and narrows the origin of the sub-leading terms to the non-perturbative structure of the foliation algebra, as described in Appendix D.4 of the paper.

## Running

```bash
python3 spectral_baseline.py
```

Requires only `numpy`. No other dependencies. Typical runtime is under two minutes.

## Parameters

| Parameter | Default | Paper value | Description |
|-----------|---------|-------------|-------------|
| `k_max`   | 80      | 80          | S³ level truncation |
| `m_max`   | 300     | ±300        | S¹ mode truncation |
| `s`       | 2       | 2           | Spectral parameter (s=2 gives Tr(\|D\|⁻⁴)) |

Convergence of the twist correction ratio was verified by the author at `k_max = 100`, `m_max = ±400` (change < 10⁻⁵).

## Note on w(k,n)

The Hopf charge decomposition `w(k,n)` is obtained by decomposing the standard Dirac multiplicity `2(k+1)(k+2)` under the Hopf U(1) action on S³. The implementation distributes the total multiplicity uniformly across integer charges `n = -k, ..., k` at each level `k`. The resulting effective weights are non-integer for `k >= 2` (e.g., `w = 4.8` at `k = 2`); this reflects equidistribution rather than a literal degeneracy count.

Multiple decompositions were tested (uniform, triangular, parity-stepped, SU(2) branching with half-integer charges). All give corrections between 5 × 10⁻⁴ and 1 × 10⁻², at least three orders of magnitude below the required irrational gain of 13.01. The conclusion is insensitive to the choice of `w(k,n)`.

## Connection to the Connes trace conjecture

This computation establishes the baseline for the Connes trace conjecture (Appendix D.4):

- Smooth Dirac spectrum contributes: 4π³ to within 10⁻³ (this computation)
- Required total (conjecture): 4π³ + π² + π ≈ 137.036
- Missing irrational gain: π² + π ≈ 13.01

Whether the Dixmier trace on the foliation C*-algebra of S³ under golden-angle Reeb flow yields the full sum is a well-posed computation in noncommutative geometry requiring the Ruelle-Sullivan current normalization. See the paper's Appendix D.4 for the precise conjecture statement.

## License

MIT

## Citation

```bibtex
@misc{down2026hopfion,
  title={Hopfion Topology, Quantum Geometric Tensor Saturation, 
         and the Fine-Structure Constant},
  author={Down, A.},
  year={2026},
  doi={10.5281/zenodo.19452185}
}
```
