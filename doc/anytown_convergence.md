# Anytown initialization and relaxation investigation

This document records the investigation of the **original upstream input**, before the user requested model changes. That version has no available supply at t=0 and cannot be repaired through initialization or relaxation alone. The actual corpus INP has subsequently been adapted into a working physical model; see [the revised model and validation](anytown_physical_model.md). wolf-3 remains excluded because the pinned EPANET reference failed at both tested accuracy settings; its provenance is retained in `tests/public_networks/excluded_networks.json`.

## Supply deficit in the original input

Pumps 78, 79 and 80 connect reservoir 40 to the network, but their speed patterns 2, 3 and 4 contain only zeros. Both tanks 41 and 42 have initial level 10 ft and minimum level 10 ft. They have no usable stored water above the minimum and cannot drain. The 19 positive-demand nodes require a total of 9800 GPM, approximately 0.618284 m3/s (618.284 kg/s).

The reservoir is separated by stopped pumps; the tanks cannot supply this deficit. Opening a stopped pump or draining a tank below its minimum would change the specified problem. The new common warning `EPANET.NO_AVAILABLE_SUPPLY` reports the number of affected nodes, total demand, first node ID and likely settings to check. The audit excludes closed links and zero-speed pumps; negative prescribed demands are counted as external injections. It ignores directional restrictions, so it is a conservative topology warning, not a complete feasibility solver.

## Numerical experiments

All tests retained the original INP bytes and all physical parameters, controls and patterns. Each hydraulic solve allowed up to 1000 iterations; continuity tolerance was 1e-8 kg/s.

| Experiment | Settings | Runs | Converged |
| --- | --- | ---: | ---: |
| Existing adaptive relaxation | Initial pressure 0, 30, 70, 100 m; uniform initial mass flow 1e-8, 1e-4, 0.01, 1, 100 kg/s; initial relaxation 1 or 0.1; growth multiplier 1 or 1.2; head tolerance 1e-12 m | 80 | 0 |
| Fixed relaxation | Initial pressure 0 or 68.58 m; initial mass flow 1e-8, 0.01, 1, 100 kg/s; relaxation 1, 0.5, 0.1, 0.01, 0.001; existing Jacobian refresh policy or refresh every iteration; head tolerance 0.1 mm | 80 | 0 |
| EPANET warm start | Independent EPANET heads/flows loaded through the existing XML initialization interface; the same five fixed relaxation factors and two Jacobian policies; head tolerance 0.1 mm | 10 | 0 |

Fixed relaxation and unconditional Jacobian refresh were tested in a temporary instrumented build, not added as production defaults. The initial adaptive-relaxation experiments use the existing compiled hydraulic core. Among the warm-start trials, the lowest final reported RMS head residual was approximately 32.4652 m, still far above 0.1 mm. These results are experimental evidence, not a proof that every conceivable numerical strategy fails; the absence of available supply is the decisive physical obstruction.

The six pipes 110, 113, 114, 115, 116 and 125 have diameter 0.0001 inches (2.54e-6 m) and length 6000 ft. With forced demand through such pipes, Hazen–Williams losses reach approximately 1e25 m. At the returned head of node 5, adjacent double-precision values are about 2.147e9 m apart. Head residuals of 0.1 mm cannot be resolved at that scale by the current unscaled double-precision formulation. These pipes must not silently be treated as closed: they connect demanded nodes in the input.

## EPANET's nominal convergence is unreliable here

The pinned EPANET 2.2 library, with 1000 trials and ACCURACY 1e-6 or 1e-5, reports only a negative-pressure warning instead of an unbalanced warning. It returns heads near -9.906e24 m, stopped pump flows of zero, and no net supply from the reservoir/tanks. An independent audit finds a continuity deficit of approximately 0.370970207 m3/s at node 20 (limit 1e-6 m3/s). This is not a usable reference solution.

Reference extraction now independently checks finite results and junction conservation after conversion to SI. For each junction, the net incoming signed link flow minus demand must be within max(1e-6 m3/s, 1e-6 times local throughput/demand). A failed audit retries the second accuracy setting and then classifies the reference as `reference_numerical_failure`. Negative pressure alone is not rejected.

The same audit also finds previously unnoticed unreliable references in `conditional_controls_2`, `control_comb`, `fcv_open_no_downstream_sources`, `prv_closed_no_upstream_sources` and `2fcvs`. These files remain diagnostic fixtures; no numerical agreement is claimed for them. All 74 previously matching strict-profile networks still match.

## Historical conclusions and saved evidence

During this investigation, no initial-condition or relaxation default was changed. A usable Anytown benchmark needs an operational supply configuration and meaningful candidate-pipe diameters. That is an input/model decision, rather than a convergence setting, and the originals were not changed. The separate `Anytown_multipointcurves.inp` already matches EPANET under the strict reference profile.

At the end of that investigation, before the physical adaptation, the standard corpus had 97 files. Normal execution: 83 solved, 13 rejected, one nonconvergent (Anytown). Strict independent comparison: 74 matches, one mismatch (ky10), four STACI numerical failures with usable references, nine unsupported cases with usable references, three EPANET-invalid inputs, and six unreliable references. Default-profile comparison: 57 matches, 22 mismatches, nine unsupported cases, three invalid inputs and six unreliable references; three of the unreliable-reference cases are solved by STACI, so they cannot be counted as numerical matches.

At that stage, all 75 other CTests passed; ky10 still failed the public numerical reference test. Detailed experiment configurations, outputs, temporary harness/build scripts and the prior nominal EPANET output are retained under `tests/test-results/anytown-numerics/`. Current public reports are under `tests/test-results/public-networks-reference/` and `public-networks-default-reference/`. The original input is archived as `doc/model_sources/Anytown_upstream.inp` (also copied into `tests/test-results/anytown-numerics/original-input.inp`). Current counts and the intentional adapted hash are documented in [anytown_physical_model.md](anytown_physical_model.md).
