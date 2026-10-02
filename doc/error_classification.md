# Current public-network error classification

PSV, PBV, emitters and PDA are now implemented. Seven previously unsupported public networks solve and match EPANET. The `io` failure is a running POWER pump at a zero-demand dead end, explained by `EPANET.POWER_PUMP_DEAD_END`. `JEP5-13` supports negative PRV settings and now matches the reference. The three infeasible valve fixtures emit `HYDRAULICS.VALVE_CONSTRAINT`; `io` retains its power-pump dead-end diagnostic. No retained model is currently blocked by a missing PSV/PBV/emitter/PDA implementation.

See [pressure elements and validation](epanet_pressure_elements.md), [ky10 reference convergence](ky10_convergence.md), and the manifest for individual expectations. A supported model can still fail numerically or have no reliable independent reference; such cases are not equivalence passes.

## Historical development-stage report


The stages 1–4 counts below are historical. Current validation covers 97 retained files, with independent junction-conservation checks on the EPANET reference; see [current reference status](epanet_reference.md).

Initial-time numerical comparison with official EPANET is part of the standard suite. Of 98 unchanged original INP files, 58 have verified numerical agreement. STACI solves 59 files; the extra one is below EPANET's minimum junction-count requirement and cannot have a numerical EPANET reference. Full definitions and reproduction commands are in [epanet_reference.md](epanet_reference.md).

## Historical completed development stages 1–4

1. Independent double precision EPANET C API reference for all original files, standard CTest/portable-runner integration, explicit required-reference mode, numerical mismatch reports and unit fixtures.
2. Invalid inputs separated from numerical cases: bad_syntax (unknown section), bad_values (unknown endpoint), and cv_controls (EPANET minimum junction-count limitation; a solvable STACI fixture, not a physically invalid network).
3. Legacy SI means LPS; two/three-field TANKS records mean fixed-head reservoirs, exactly as in EPANET. Inputs remain unchanged; conversion is internal and warned. The 57,460-pipe case now converges and matches EPANET. JEP5-13 progresses to an explicit unsupported PRV diagnostic.
4. Corrected EPANET pipe models (H-W exponents/coefficient, low-flow regularization, minor-loss derivative, D-W viscosity and transition law), pump power units and power/piecewise head curves, default demand pattern 1, supported initial-time controls and tank limits. Replaced dimensionless ACCURACY-to-dimensional residual mapping with tight physical tolerances. System construction now uses indexed lookup, making the largest case practical. Machine-readable steady results are exported without rewriting INP inputs; stale hydraulic exports are removed before a new solve.

No new valve, emitter or pressure-dependent demand model was added.

## Historical first-error categories (before valve/emitter/PDA support)

Counts below are disjoint first-rejection categories; a network can contain further unsupported features after the first is resolved.

| Category | Cases | Required work |
| --- | ---: | --- |
| Unsupported PRV, PSV, PBV, FCV, GPV valves | 27 | New hydraulic models and state transitions; TCV already supported |
| Nonzero emitters | 3 | Pressure-dependent outlet model |
| Pressure-dependent demand (PDA) | 2 | Demand model / nonlinear solver integration |
| Initially unfed components | 2 | Control scheduling and connectivity investigation; do not open links silently |
| Unknown endpoint / unknown section | 2 | Correctly rejected negative fixtures; preserve the originals |
| STACI non-convergence | 3 | Anytown, NW_Model, NW_Model1: conditioning/check-valve/solver investigation |

The 36 compatibility/input rejections include the two negative fixtures and wolf-3. The official reference is unbalanced on wolf-3 under both configured convergence settings; this is separately reported and is not an equivalence pass. EPANET rejects the smaller cv_controls fixture with code 223; STACI's successful mathematical solve does not establish agreement with EPANET.

The ky16 EPS problem improved from 16 failed frames to one frame at the strict 1e-12 m RMS setting with the model corrections. The configurable default is now 0.1 mm RMS and ky16 passes all 25 frames, including the portable integration test. Resolving the remaining strict-precision failure would require further solver work. Original, strict and default metadata are retained in test-results. No test is counted as a solved/equivalent network solely because an expected error was produced.

See [diagnostics.md](diagnostics.md) for the common GUI error protocol and [the manifest](../tests/public_networks/manifest.json) for each network's independent-reference expectation.

## Historical PRV/FCV/GPV implementation stage

At this historical stage, PRV, FCV and GPV had been added. The then-remaining compatibility classes included PSV, PBV, emitter/PDA and disconnected or invalid inputs. Normal execution is now 84 solved / 13 rejected / zero nonconvergent; wolf-3 was removed and Anytown explicitly adapted into a physical model. Strict comparison is 75 matches and one solved mismatch (ky10), plus the classified numerical/input/reference failures. See [epanet_valves.md](epanet_valves.md) for exact scope; ky10 still failed at that stage; it now matches, as described at the top of this document.

The [physical Anytown model](anytown_physical_model.md) resolves the previously unavailable supply and candidate-pipe conditioning through intentional input changes. Zero-speed pumps remain supported and are independently tested.
