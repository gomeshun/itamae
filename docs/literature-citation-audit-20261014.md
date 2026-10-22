# Literature citation audit: itamae — started 2026-10-14 UTC

Inventoried **476 active definitions**, including inherited/moved functions, private/nested helpers, constructors and tests. Updated **16 docstrings**, of which **16** receive or complete local paper references. The remaining edits clarify units, assumptions or NumPy structure.

The executable baseline before this citation edit is `c452c5bbbbf8276870a8ff1f61f77d8bff0d8d7c`. Comparison baseline: new shared library; no author-upstream baseline. Moving a definition does not make it new; this citation pass also includes inherited definitions.

The audit screens numerical literals and existing attributions, then matches adopted values to primary equations/tables. Generic algorithms, numerical controls, unit conversions and frozen regression outputs do not receive invented physical calibrations. The full [definition index and hash evidence](literature-citation-audit-20261014.json) distinguish edited definitions, unchanged references and archived evidence. An inventory entry is not proof of physical calibration.

## Findings and limits

- Native defaults use rounded Planck-2018 values; the Colossus preset is identified separately. A background choice does not recalibrate halo fits.
- Generic protocols, callback integration, polynomial expansions, grid controls and unit conversion identities do not acquire invented empirical calibrations.

## Verification

- All 16 changed docstrings parse as NumPy style, cover explicit arguments and contain their assigned paper identifiers.
- All 74 inventoried tracked Python files have identical executable ASTs after removing only docstrings; assertions, tolerances and equations are unchanged. 0 historical Python files remain byte-identical.
- 333 distinct existing local pytest cases pass on Python 3.11.15. This is source-checkout testing, not a new scientific catalogue or a new calibration.
- Full Ruff lint and formatting checks pass.
- Historical numerical results keep their original calculation revisions. The docstring edit changes source bytes and build identity; executable equivalence does not permit relabeling old outputs. No tests are deleted or weakened.

## Changed definitions and paper mapping

| Definition | References / change |
| --- | --- |
| [`ColossusCosmology`](../src/itamae/cosmology/colossus.py#L11) | [colossus](https://arxiv.org/abs/1712.04512), [planck](https://arxiv.org/abs/1807.06209) |
| [`ColossusCosmology.__post_init__`](../src/itamae/cosmology/colossus.py#L42) | [colossus](https://arxiv.org/abs/1712.04512), [planck](https://arxiv.org/abs/1807.06209) |
| [`ColossusCosmology.collapse_threshold`](../src/itamae/cosmology/colossus.py#L233) | [ludlow](https://arxiv.org/abs/1601.02624) |
| [`NativeFlatLCDM`](../src/itamae/cosmology/native.py#L75) | [planck](https://arxiv.org/abs/1807.06209) |
| [`NativeFlatLCDM.__post_init__`](../src/itamae/cosmology/native.py#L104) | [planck](https://arxiv.org/abs/1807.06209) |
| [`NativeFlatLCDM.growth_factor`](../src/itamae/cosmology/native.py#L291) | [growth](https://doi.org/10.1146/annurev.aa.30.090192.002435) |
| [`NativeFlatLCDM.collapse_threshold`](../src/itamae/cosmology/native.py#L330) | [ludlow](https://arxiv.org/abs/1601.02624) |
| [`nfw_mass_function`](../src/itamae/halo/nfw.py#L41) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`invert_nfw_mass_function`](../src/itamae/halo/nfw.py#L84) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`NFWProfile`](../src/itamae/halo/nfw.py#L191) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`NFWProfile.__post_init__`](../src/itamae/halo/nfw.py#L220) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`NFWProfile.enclosed_mass`](../src/itamae/halo/nfw.py#L243) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`NFWProfile.density`](../src/itamae/halo/nfw.py#L277) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`NFWProfile.potential`](../src/itamae/halo/nfw.py#L313) | [nfw](https://arxiv.org/abs/astro-ph/9611107) |
| [`test_colossus_backend_matches_order_of_magnitude`](../tests/test_cosmology.py#L28) | [planck](https://arxiv.org/abs/1807.06209) |
| [`test_native_cosmology_identifies_distinct_float_parameters`](../tests/test_public_input_contracts.py#L382) | [planck](https://arxiv.org/abs/1807.06209) |

## Bibliography

- **colossus**: [Diemer (2018), COLOSSUS: A Python Toolkit for Cosmology, Large-scale Structure, and Dark Matter Halos.](https://arxiv.org/abs/1712.04512)
- **growth**: [Carroll, Press & Turner (1992), The Cosmological Constant, Section 3.5.](https://doi.org/10.1146/annurev.aa.30.090192.002435)
- **ludlow**: [Ludlow et al. (2016), The mass-concentration-redshift relation of cold and warm dark matter haloes, Appendix C.](https://arxiv.org/abs/1601.02624)
- **nfw**: [Navarro, Frenk & White (1997), A Universal Density Profile from Hierarchical Clustering.](https://arxiv.org/abs/astro-ph/9611107)
- **planck**: [Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.](https://arxiv.org/abs/1807.06209)
