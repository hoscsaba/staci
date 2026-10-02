# EPANET numerical equivalence

Validation status updated 2026-10-02. The public initial-state corpus contains
97 retained INP files: 96 unchanged upstream files and one explicitly adapted
physical Anytown. Original/adapted hashes are recorded in the manifest.

The latest initial-state reference run has **88 numerical matches, three EPANET
input rejections and six unreliable EPANET references**. Standalone execution
under the strict profile has 89 solved cases, four compatibility/input rejections
and four failed hydraulic cases. Expected failures are regression passes, not
numerical matches. `ky10`, `JEP5-13`, `NW_Model`, `NW_Model1` and `GES4-9` now match;
older reports describing them as open failures are historical.

## Reproduce

```sh
python3 tests/setup_epanet_reference.py
cmake -S . -B build -DSTACI_EPANET_EXECUTABLE=/absolute/path/runepanet -DSTACI_EPANET_LIBRARY=/absolute/path/libepanet2.so
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 tests/test_public_networks.py --binary build/staci --reference-library /absolute/path/libepanet2.so --require-reference --output-dir build/public-reference
python3 tests/test_epanet_eps_corpus.py --binary build/staci --library /absolute/path/libepanet2.so --output build/eps-reference --timeout 300
```

Use `libepanet2.dylib` on macOS or `epanet2.dll` on Windows. The pinned official
OWA EPANET 2.2 revision is `4d8d82ddc260fce216af9321fc3d9a4646ac6827`.
Configure both reference paths to enable all comparisons; missing tools produce
skips, not agreement. The shared C library is queried in double precision, with
serialized calls because its input tokenizer has shared state.

## Quantities and tolerances

IDs are original EPANET IDs; internal boundary edges are excluded. Both outputs
are converted to SI. A numerical comparison passes against the larger absolute
or relative limit. Tolerances are common across the corpus.

| Quantity | Absolute tolerance | Relative tolerance |
|---|---|---|
| Total head / junction pressure head | 0.01 m | 1e-6 |
| Prescribed junction demand | 1e-8 m³/s | 1e-4 |
| Reservoir/tank net exchange | 1e-6 m³/s | 1e-4 |
| Signed link flow | 1e-6 m³/s | 1e-4 |
| Velocity magnitude | max(1e-5 m/s, flow absolute tolerance / area) | 1e-4 |
| Endpoint head difference | 0.02 m | 1e-5 |
| Open/closed state | exact | exact |
| EPS chemical concentration | 1e-5 kg/m³ (0.01 mg/L) | 3% |
| EPS water age | 60 s | 3% |

The steady reference harness requests STACI head residual 1e-12 m and continuity
1e-8 kg/s. The full-period EPS harness requests 1e-9 m and 1000 iterations.
Ordinary INP execution defaults to 0.1 mm RMS. Residual convergence and agreement
with another solver are different tests.

Reference acceptance additionally requires finite results, a head-equation error
below 0.1 mm, independent junction continuity and physical constant-power pump
operating points. Retried numerical settings are retained in report metadata.
`wolf-3` was removed after its reference failed; `cv_controls`, `bad_syntax` and
`bad_values` remain deliberately non-comparable input cases.

## Extended periods and chemistry

A steady match does not establish EPS agreement. The full-period runner records
hydraulic and quality outcomes separately. Remaining full-period hydraulic
mismatches include Net3 and Net6 variants; water-quality agreement there is not
proof of overall compatibility. TRACE and water-age tank mixing remain limited.

The chemical implementation now includes first-order bulk decay, flow-dependent
wall mass transfer (Schmidt/Sherwood correlations), complete-mix tank storage and
reaction, downstream initial pipe concentrations, concentration-tolerance segment
merging, and reactive stagnant-node concentrations. A zero input minimum tank
volume means the cylindrical volume at minimum level, as in EPANET. Tank storage
therefore includes water below the operating minimum.

`epanet_eps_chemical_reference` compares **Net1 (25 report frames), Net2
(56 frames) and Net2-CL2 (56 frames)** with unchanged inputs and tolerances. Both hydraulics and chemistry
pass. Net1 formerly had 188 mismatching chemical values; Net2 had 580. Net1 covers
bulk/wall decay and tank cycling; Net2 isolates reaction-free tracer transport,
a negative-demand concentration source and tank dilution.

Chemical EPS supports MIXED tanks and first-order reactions. Unsupported mixing
models (2COMP/FIFO/LIFO), reaction orders, limiting potential and roughness
correlations produce an explicit compatibility error in chemical EPS. They do
not block an initial hydraulic solve. MSX multi-species reactions are outside
this single-species INP comparison, even for an INP located in an MSX example.

Numeric `[TIMES]` values honor explicit seconds/minutes/hours/days. The
Net2-CL2 regression covers `Quality Timestep 5.00 min`, previously interpreted
as five hours.

Remaining chemical discrepancies are reported explicitly. The stagnant branch
of `prv_open_no_upstream_sources` retains a node-concentration discrepancy:
small accepted hydraulic residual circulations can select different incoming
segments even though hydraulic quantities meet their comparison tolerances.
Net6_plus and Net3 variants also have hydraulic discrepancies, so their chemical
mismatches cannot be interpreted independently of those flows.

The Batch-NH2Cl INP exposes a pinned-reference edge case: with all flows stagnant,
EPANET can leave reported tank quality unchanged despite nonzero bulk decay.
The reference harness now marks this chemical reference unavailable and records
`QUALITY.STAGNANT_TANK`, while retaining its valid hydraulic comparison. STACI's
isolated-tank reaction is checked separately against the discrete mass balance.
An earlier claim of Batch chemical equivalence used unchanged tank quality on
both sides and did not validate reaction physics.

`ky16` still has one failed report frame with the extreme 1e-12 m/50-iteration
profile; all 25 frames pass at the normal 0.1 mm setting. `solver_options` verifies
this bounded partial-result behavior (or a fully converged improvement), the
reported strict limit and exit code 3. Passing that regression does not mean the
strict ky16 simulation is fully solved.

The Anytown tank-limit correction closes the incident inward/outward links while
retaining the tank head. Whole-second event times are rounded upward so an event
cannot occur just before a level threshold and then be missed until the next
grid point. Earlier measurements reduced a 14.4 mm tank discrepancy to about
0.42 mm; a tiny flow difference can still exceed the strict comparison limit.

## Physical failures and integration

The four failed steady cases are documented in [valve validation](epanet_valves.md)
and [pressure elements](epanet_pressure_elements.md). They have incompatible
flow/supply constraints or a constant-power dead end. Relaxed residual tolerances
must not make numerical regularization count as a physical source.

Reports retain per-network settings, mismatches, warnings and unavailable-reference
reasons. See [README_DEV](../README_DEV.md) for external calls, result handling and
JSONL errors. Historical investigations are [Anytown](anytown_convergence.md),
[its explicit physical adaptation](anytown_physical_model.md), and [ky10](ky10_convergence.md).
