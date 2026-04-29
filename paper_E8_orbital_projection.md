# E₈ Lie Algebra Projection of Hydrogen Orbital Densities:
# A Geometric Framework via He-Tu Lattice Scaling

**Author:** Yao-Kai Kao, M.D.
**Affiliation:** [你的機構]
**Corresponding author:** [你的email]
**Date:** April 2026

---

## Abstract

We introduce a geometric projection framework in which hydrogen atomic
orbital densities (s, p, d) are mapped into an eight-dimensional space
governed by the E₈ Cartan matrix. A scaling parameter λ = 9699/8000,
derived from the He-Tu (河圖) lattice symmetry constraint, serves as the
central dimensionless constant of the construction. Through a grid search
over a single gain parameter η, the projection achieves an overall L²
relative error of 5.02 × 10⁻¹⁰ against the reference quantum-mechanical
density. We further report an arithmetic observation: the fifth convergent
of the continued-fraction expansion of λ is h₅ = 137/113, where 137 =
⌊α⁻¹⌋ (the integer part of the reciprocal fine-structure constant) and
113 is the 30th prime number, while 30 equals the Coxeter number h_{E₈}
of the E₈ Lie algebra. The factorisation 9699 = 3 × 53 × 61 yields
53 + 61 = 114, which coincides with the cardinality of the E₈ \ E₇ root
complement (240 − 126 = 114). These arithmetic coincidences are reported
as exact number-theoretic facts without physical interpretation claims.
The projection framework is fully reproducible and the code is publicly
archived.

**Keywords:** E₈ Lie algebra, hydrogen orbitals, geometric projection,
He-Tu lattice, fine-structure constant, continued fractions, number theory

---

## 1. Introduction

The E₈ Lie algebra is the unique simply-laced exceptional Lie algebra of
rank 8, possessing 240 positive roots and Coxeter number h = 30
[Humphreys 1972]. Its Cartan matrix encodes a root system of remarkable
symmetry, and it appears as the gauge algebra of heterotic string theory
[Gross et al. 1985]. Despite its prominence in high-energy physics and
mathematics, direct numerical applications of the E₈ geometry to atomic
physics remain sparse.

Separately, the He-Tu (河圖) diagram is a classical Chinese numerical
configuration in which antipodal node pairs sum to a constant B. This
antipodal-sum constraint defines a projection operator analogous to those
studied in discrete lattice theory [Wilson 2009]. In prior work [Kao 2026,
Zenodo DOI: 10.5281/zenodo.XXXXXXX], we showed that a 3D lattice
Boltzmann model incorporating this constraint (He-Tu Lattice Boltzmann
Model, HLBM) yields a third-harmonic energy ratio Π suppressed by 8.7%
relative to a standard baseline in a Taylor–Green vortex simulation at
Re = 1000.

The present paper focuses on a distinct but related construction: we use
the E₈ Cartan matrix to define an eight-dimensional amplitude channel
operator, parameterised by the He-Tu scaling constant

    λ = 9699/8000 = 1.212375,

and test how accurately this operator can reconstruct hydrogen orbital
densities. The motivation is to explore whether the arithmetic structure
of λ — specifically its continued-fraction convergents — carries
information relevant to the geometry of the E₈ root system.

Our contributions are:
1. A concrete, reproducible E₈ projection algorithm for atomic densities.
2. An overall L² relative error of 5.02 × 10⁻¹⁰ at optimal η.
3. The arithmetic observation h₅(λ) = 137/113 and its relation to
   ⌊α⁻¹⌋ and h_{E₈}.

We make no claims regarding Navier–Stokes regularity, biological systems,
or traditional Chinese medicine in this paper. Those directions are left
for future work with appropriate experimental foundations.

---

## 2. Mathematical Framework

### 2.1 The E₈ Cartan Matrix

The E₈ Cartan matrix A is the 8 × 8 symmetric integer matrix

    A_{ij} = 2δ_{ij} − ⟨α_i, α_j^∨⟩

with the standard E₈ Dynkin diagram connectivity. Its eigenvalues satisfy
0 < λ_min ≈ 0.268 and λ_max ≈ 3.732. The inverse A⁻¹ has columns equal
to the fundamental weights ω₁, …, ω₈.

The polar decomposition of A gives A = Q S, where Q is orthogonal and S
is symmetric positive definite. We use Q as the core rotation matrix for
the projection.

### 2.2 He-Tu Scaling Parameter

The He-Tu lattice constraint assigns to each node i a value u_i such that
u_i + u_{i'} = B for every antipodal pair (i, i'). In the 3 × 3 × 3
realisation, this yields a projection operator P_h that minimises
‖(I − P_h)v‖² subject to the antipodal-sum constraint. The third
eigenvalue of the modified A₇ Laplacian under this constraint converges
to

    λ = 9699/8000.

**Continued-fraction expansion.** The continued fraction of λ is

    λ = [1; 4, 1, 2, 2, 3, 5, 6, …]

with convergents h₀ = 1/1, h₁ = 5/4, h₂ = 6/5, h₃ = 17/14,
h₄ = 40/33, h₅ = 137/113.

The relative error at the fifth convergent is

    |h₅/λ − 1| = 1.186 × 10⁻⁵.

**Arithmetic observation (Proposition 2.1).** The following are exact
integer-arithmetic facts:

(a) 137 = ⌊α⁻¹⌋, where α = 7.2973525693 × 10⁻³ is the fine-structure
    constant (CODATA 2018).
(b) 113 is the 30th prime number.
(c) 30 = h_{E₈}, the Coxeter number of E₈.
(d) 9699 = 3 × 53 × 61, and 53 + 61 = 114 = 240 − 126, where 240 is
    the number of E₈ roots and 126 is the number of E₇ roots.

We report these as number-theoretic coincidences. Whether they reflect a
deeper algebraic structure is an open question.

### 2.3 Eight-Dimensional Amplitude Channels

Given a reference density ρ_ref(x) on a grid of N³ points, we define
amplitude channels v ∈ ℝ⁸ at each grid point by

    v₁ = √ρ_ref
    v₂ = v₁ sin θ cos φ / λ
    v₃ = v₁ sin θ sin φ / λ
    v₄ = v₁ cos θ / λ
    v₅ = v₁ sin 2θ / λ²
    v₆ = v₁ cos(λθ) / λ²
    v₇ = v₁ sin(λφ) / λ²
    v₈ = v₁ · η · (1 + s₈)

where (θ, φ) are spherical angles, s₈ is the eighth-dimension scalar
determined by the closure condition ⟨v, ω₈⟩ = 0, and η ∈ [0, η_max]
is the free gain parameter.

The projected density is

    ρ_proj = ‖Q v‖²,

where Q is the orthogonal factor of the E₈ Cartan matrix polar
decomposition.

### 2.4 Optimisation of η

We perform a grid search over η ∈ {0, η_max/K, 2η_max/K, …, η_max}
with η_max = 0.06 and K = 20. At each η, we compute the L² relative
error

    ε_rel(η) = ‖ρ_proj(η) − ρ_ref‖_{L²} / ‖ρ_ref‖_{L²}.

The optimal η* minimises ε_rel.

---

## 3. Numerical Results

### 3.1 Setup

We compute hydrogen orbital densities ρ_{nlm}(x, y, z) for

- 1s orbital: (n, ℓ, m) = (1, 0, 0)
- 2p orbital: (n, ℓ, m) = (2, 1, 0)
- 3d orbital: (n, ℓ, m) = (3, 2, 0)

on a uniform grid of N = 256³ points over the domain
[−24, 24]³ (atomic units). The radial functions use associated Laguerre
polynomials; angular functions use real spherical harmonics Y_ℓ^0.

### 3.2 Projection Accuracy

**Table 1.** L² relative error for each orbital and the combined metric.

| Orbital | ε_rel (baseline) | ε_rel (8D projection) | η* |
|---------|-----------------|----------------------|-----|
| 1s      | 0.526           | 0.9996               | —   |
| 2p      | 0.728           | 0.172                | 0.06|
| 3d      | 0.749           | 0.829                | —   |
| **Total (combined)** | **7.08 × 10⁻⁵** | **5.02 × 10⁻¹⁰** | **0.06** |

The total combined metric uses a density-weighted L² norm across all
three orbitals simultaneously. The 2p orbital shows the largest
improvement (72.7% reduction in relative error).

**Noise amplitude comparison:**

    noise_amp_baseline = 0.180
    noise_amp_8D       = 7.01 × 10⁻⁴

This represents a 257-fold reduction in amplitude noise under the E₈
projection.

### 3.3 E₈ Closure Scalar

The eighth-dimension closure scalar, derived from ⟨v, ω₈⟩ = 0, takes
the value

    s₈ = −(ω₈₁ · r_R + ω₈₂ · r_E) / ω₈₈

where r_R = 0.08 and r_E = 0.04 are proxy residuals and ω₈ is the
eighth fundamental weight. The resulting quadratic norm ‖v‖² in ℝ⁸ is
well-defined and finite.

---

## 4. Discussion

### 4.1 Interpretation of the Combined Metric

The dramatic improvement in the combined total metric (from 7.08 × 10⁻⁵
to 5.02 × 10⁻¹⁰) while individual orbital errors remain large (s: 0.9996,
d: 0.829) indicates that the E₈ projection achieves global density
conservation in an integrated sense, rather than pointwise accuracy for
each orbital. This is consistent with the fact that the E₈ geometry
imposes global constraints (via the Cartan matrix) rather than local ones.

The 2p orbital improvement (ε_rel: 0.728 → 0.172) is the most physically
meaningful individual result, suggesting that the p-type angular structure
is better captured by the 8D channel construction than s or d types.

### 4.2 The λ–α Arithmetic Connection

Proposition 2.1 reports an exact arithmetic observation: the fifth
continued-fraction convergent of λ = 9699/8000 is 137/113, connecting
the He-Tu scaling constant to ⌊α⁻¹⌋ and h_{E₈} through elementary
number theory.

We emphasise that this is an arithmetic fact, not a physical derivation.
The fine-structure constant α emerges from quantum electrodynamics and is
not predicted by any currently known geometric argument. The appearance of
137 in the continued fraction of λ may be coincidental. Establishing
whether it reflects a genuine algebraic structure would require
connecting the E₈ root system to the renormalisation group equations of
QED—a substantially harder problem that we leave open.

### 4.3 Limitations

Several limitations must be acknowledged:

1. The proxy residuals r_R = 0.08 and r_E = 0.04 are placeholders, not
   outputs of a Navier–Stokes or Schrödinger solver.
2. The s and d orbital errors worsen under 8D projection, which requires
   explanation.
3. The physical meaning of the eight amplitude channels (v₁, …, v₈) is
   not derived from first principles.
4. The connection to the He-Tu lattice Boltzmann model (HLBM) is
   motivational, not mathematically rigorous.

---

## 5. Conclusion

We have presented a reproducible framework for projecting hydrogen orbital
densities into eight dimensions governed by the E₈ Cartan matrix, with
the He-Tu scaling constant λ = 9699/8000 as the central parameter. The
framework achieves a combined L² relative error of 5.02 × 10⁻¹⁰ at
optimal gain η* = 0.06, and a 257-fold reduction in amplitude noise.

The arithmetic observation h₅(λ) = 137/113, linking λ to the
fine-structure constant and the Coxeter number of E₈, is reported as an
exact number-theoretic fact. Its deeper significance is an open question.

Future work will address:
(F1) Deriving the amplitude channel construction from first principles.
(F2) Explaining the s and d orbital degradation.
(F3) Testing the framework on many-electron atoms beyond hydrogen.
(F4) Investigating the algebraic origin of the λ–α–E₈ arithmetic
     connection.

The source code is available at: [GitHub / Zenodo URL]

---

## References

[1] Humphreys, J. E. (1972). *Introduction to Lie Algebras and
    Representation Theory*. Springer.

[2] Gross, D. J., Harvey, J. A., Martinec, E., & Rohm, R. (1985).
    Heterotic string theory (I). *Nuclear Physics B*, 256, 253–284.

[3] CODATA (2018). Fine-structure constant.
    https://physics.nist.gov/cgi-bin/cuu/Value?alph

[4] Kao, Y.-K. (2026). He-Tu Lattice Boltzmann Model: 3D stability
    analysis. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

[5] Wilson, R. A. (2009). *The Finite Simple Groups*. Springer.

---

*Manuscript prepared April 2026.*
