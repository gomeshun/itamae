# Actual-paper citation verification: itamae

Audit date: 2026-10-15 (JST). This supersedes the equation-attribution conclusions of the 2026-10-14 citation index; the old index remains a historical snapshot.

Inspected source revision: `4aeb13ed4e3e3e87171fe800e08f74b95c1244a2`. All cited definitions include private helpers, nested functions and tests. A moved function is traced through its retained qualified name; a move does not imply a new physical prescription.

Comparison baseline: `none` at `not applicable` (new shared library).

Coverage: 16 definitions, 18 original reference occurrences; 18 after correction. 5 reference keys / 5 unique papers. Every original occurrence has a reviewed source comparison, including citations removed as inapplicable.

The linked full texts were actually opened and the claimed equations, coefficients or contextual values inspected. Publisher landing-page reachability is separate from full-text access. Equation existence is separate from algebraic derivation, implementation adaptation and physical validation. PDF page numbers refer to the inspected version, not necessarily printed page numbers.

Categories: `printed` = applicable relation/value appears in the paper; `derived` = a documented algebraic/numerical consequence; `adopted` = reuse with explicit assumptions; `adapted` = an implementation convention or coefficient set differs; `context` = source supplies background rather than an equation; `unresolved` = the exact attribution is not established. Counts are citation uses, not counts of independent physical validations.

adopted: 6, context: 2, derived: 7, printed: 3

## Source comparisons

| Key | Primary full text and inspected location | Result |
|---|---|---|
| planck | [Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.](https://arxiv.org/pdf/1807.06209) — 1807.06209v4 [astro-ph.CO] 9 Aug 2021; Table 2, p.16 (TT,TE,EE+lowE+lensing column) | **adopted**. H0=67.36, Omega_m=0.3153, omega_b=0.02237, tau=0.0544 and A_s=2.100e-9 support the rounded native/control inputs. Omega_b is obtained by dividing omega_b by h^2. The Colossus planck18 preset selects a different data combination; its name does not make it identical to the native parameter vector. |
| colossus | [Diemer (2018), COLOSSUS: A Python Toolkit for Cosmology, Large-scale Structure, and Dark Matter Halos.](https://arxiv.org/pdf/1712.04512) — 1712.04512v3 [astro-ph.CO] 18 Dec 2018; Table 1, p.2; Sections 2.1-2.2 | **context**. The v3 paper lists planck18 (H0=67.66, Omega_m=0.3111, Omega_b=0.0490) and planck18-only (67.36, 0.3153, 0.0493) separately. It describes named presets and global cosmology state. Saving/restoring that state is this wrapper's implementation, not a published physical equation. |
| growth | [Carroll, Press & Turner (1992), The Cosmological Constant, Section 3.5.](https://ned.ipac.caltech.edu/level5/Carroll/Carroll3_5.html) — NED full-text reprint; Section 3.5, Eq. (29), NED full-text reprint | **derived**. The displayed growth suppression is g=2.5*Omega_m/[Omega_m^(4/7)-Omega_Lambda+(1+Omega_m/2)*(1+Omega_Lambda/70)]. The normalized D(z)=g(z)/[(1+z)*g(0)] and code redshift derivatives are derived uses. The DOI reaches a publisher landing page; the NED equation image provides accessible full text. |
| ludlow | [Ludlow et al. (2016), The mass-concentration-redshift relation of cold and warm dark matter haloes, Appendix C.](https://arxiv.org/pdf/1601.02624) — 1601.02624v2 [astro-ph.CO] 10 May 2016; Eq. (3) and following text, p.5; Appendix C, pp.18-19, Eqs. (C8)-(C12) | **printed**. The paper explicitly uses delta_sc(z)=1.686/D(z). The variance fit is sigma=22.26*xi^0.292/(1+1.53*xi^0.275+3.36*xi^0.198), xi=(M/(1e10 h^-1 Msun))^-1. It reports about 2.5% fit accuracy over 1e-7 to 1e15 h^-1 Msun. Differentiation is derived. WDM host histories integrate their spectrum and do not use this CDM variance fit; their redundant citation was removed. |
| nfw | [Navarro, Frenk & White (1997), A Universal Density Profile from Hierarchical Clustering.](https://arxiv.org/pdf/astro-ph/9611107) — astro-ph/9611107v4 21 Oct 1997; Eq. (1), p.2; Eqs. (2)-(3), p.7 | **derived**. The density and circular-speed laws are printed. f(x)=ln(1+x)-x/(1+x) follows by integrating the density. Potential with zero at infinity, exact mass-conversion root finding, and c=-1/W0(-exp[-(y+1)])-1 are derived consequences; the paper does not print the implementation's inverse or solver. The 2.163 and 4.625 factors are rounded maximum-speed identities. |

The companion JSON records URL status, UTC retrieval time and the SHA-256 of each full-text file. Original identifiers were retained where an alternative full-text link was added. No third-party paper files are redistributed with this report.

## Definition coverage

| Definition at inspected revision | Original baseline definition | References before → after | Classification after correction |
|---|---|---|---|
| [src/itamae/cosmology/colossus.py::ColossusCosmology](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/colossus.py#L11) | no author baseline: new shared library | colossus, planck | colossus: context; planck: adopted |
| [src/itamae/cosmology/colossus.py::ColossusCosmology.__post_init__](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/colossus.py#L42) | no author baseline: new shared library | colossus, planck | colossus: context; planck: adopted |
| [src/itamae/cosmology/colossus.py::ColossusCosmology.collapse_threshold](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/colossus.py#L233) | no author baseline: new shared library | ludlow | ludlow: printed |
| [src/itamae/cosmology/native.py::NativeFlatLCDM](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/native.py#L75) | no author baseline: new shared library | planck | planck: adopted |
| [src/itamae/cosmology/native.py::NativeFlatLCDM.__post_init__](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/native.py#L104) | no author baseline: new shared library | planck | planck: adopted |
| [src/itamae/cosmology/native.py::NativeFlatLCDM.growth_factor](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/native.py#L291) | no author baseline: new shared library | growth | growth: derived |
| [src/itamae/cosmology/native.py::NativeFlatLCDM.collapse_threshold](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/cosmology/native.py#L330) | no author baseline: new shared library | ludlow | ludlow: printed |
| [src/itamae/halo/nfw.py::nfw_mass_function](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L41) | no author baseline: new shared library | nfw | nfw: derived |
| [src/itamae/halo/nfw.py::invert_nfw_mass_function](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L84) | no author baseline: new shared library | nfw | nfw: derived |
| [src/itamae/halo/nfw.py::NFWProfile](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L191) | no author baseline: new shared library | nfw | nfw: derived |
| [src/itamae/halo/nfw.py::NFWProfile.__post_init__](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L220) | no author baseline: new shared library | nfw | nfw: derived |
| [src/itamae/halo/nfw.py::NFWProfile.enclosed_mass](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L243) | no author baseline: new shared library | nfw | nfw: derived |
| [src/itamae/halo/nfw.py::NFWProfile.density](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L277) | no author baseline: new shared library | nfw | nfw: printed |
| [src/itamae/halo/nfw.py::NFWProfile.potential](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/src/itamae/halo/nfw.py#L313) | no author baseline: new shared library | nfw | nfw: derived |
| [tests/test_cosmology.py::test_colossus_backend_matches_order_of_magnitude](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/tests/test_cosmology.py#L28) | no author baseline: new shared library | planck | planck: adopted |
| [tests/test_public_input_contracts.py::test_native_cosmology_identifies_distinct_float_parameters](https://github.com/gomeshun/itamae/blob/4aeb13ed4e3e3e87171fe800e08f74b95c1244a2/tests/test_public_input_contracts.py#L382) | no author baseline: new shared library | planck | planck: adopted |

## Source and validation boundaries

Compared all 74 tracked Python files in this checkout after stripping only leading docstrings: executable syntax is unchanged. All 0 sealed historical Python files are byte-identical. No model coefficients, numerical controls, assertions or tolerances were altered.

Original qualified-name links record the author/owner baseline even when the current file or containing class moved. A missing one-to-one name mapping is reported as such; it is not evidence that a physical equation is new. The companion JSON preserves baseline source-file hashes.

The retained test results from the previous citation audit belong to their recorded source revisions. They are not relabelled as new physical runs. This audit checks documentation, paper correspondence and executable-source identity; it does not establish convergence, observational calibration or validity outside the published domains. Source-file hashes can change when documentation changes, so persisted output/restart identity checks must not be bypassed.
