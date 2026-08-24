# Channel reference test suite

`run_channel_tests.py` exercises only the circular `channel1` element. Each
network is generated from `channel.spr`; the original file is never modified.
All calculations and reported values use SI units.

The suite covers:

- positive uniform flow at normal depth;
- positive profiles below normal depth, between normal and critical depth,
  and above critical depth;
- reverse profiles below and above critical depth;
- a supercritical-to-subcritical hydraulic jump with circular-section momentum
  matching;
- open-to-full mixed profiles in both flow directions;
- fully pressurised flow in both directions;
- a hydrostatic zero-flow state.

## Independent references

The Python oracle does not call STACI hydraulic functions. It implements:

- exact circular-segment area, wetted perimeter and top width;
- Manning's equation for normal depth;
- `Fr² = Q² T / (g A³) = 1` for critical depth;
- `dy/dx = (S0 - Sf) / (1 - Fr²)` with an independent RK4 integration for
  gradually varied profiles;
- Darcy-Weisbach for completely pressurised cases;
- `M = Q²/(gA) + I₁`, with `I₁ = ∫(y-z)T(z)dz`, independently evaluated on
  both sides of the hydraulic jump.

These equations and assumptions follow the following public engineering
references:

- [FHWA HDS-3: Design Charts for Open-Channel Flow](https://www.fhwa.dot.gov/engineering/hydraulics/pubs/hds3.pdf)
- [USACE HEC-RAS: equations for basic profile calculations](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.2/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-steady-flow-water-surface-profiles/equations-for-basic-profile-calculations)
- [USACE HEC-RAS: critical-depth determination](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/latest/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-steady-flow-water-surface-profiles/critical-depth-determination)
- [FHWA HEC-22, fourth edition](https://www.fhwa.dot.gov/engineering/hydraulics/pubs/hif24006.pdf)
- [EPA SWMM Hydraulics Addendum: open/full transition and Preissmann slot](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1014FBE.txt)

The mixed open/full checks are explicitly marked as *model-consistency* tests.
The present STACI steady model has no Preissmann slot or pressure-wave model,
so these cases cannot be treated as independent validation of the transition's
physics.

When two open ends require a supercritical-to-subcritical transition, a smooth
GVF profile is impossible. STACI uses the signed diffusive-wave relation as the
additional steady discharge closure, integrates the supercritical and
subcritical profiles from their controlling boundaries, and locates their
discontinuity by equal specific momentum. The `hydraulic-jump-positive` case
checks both the discharge and the momentum residual independently. End depths
alone still do not define a unique discharge and jump location, which is why
the discharge closure is stated explicitly.

An exact gradually-varied profile cannot cross `Fr = 1`, because the governing
ODE is singular there. The suite therefore computes critical depth independently
and tests profiles on both sides of it. On the unchanged descending bed, a
reverse flow has an adverse bed slope in its own flow direction, so a Manning
normal depth is not physically defined for that direction; reverse cases are
classified relative to critical depth instead.

## Running

Build STACI first, then run only this suite:

```bash
python3 tests/run_channel_tests.py --binary build/staci
```

Results are retained in `tests/test-results/channel/`. The human-readable
`run_channel_tests.log` explains every case, while `results.csv` is convenient
for automated comparison. Every case directory also contains two longitudinal
sections, `channel-profile.svg` and `channel-profile.pdf`, with the calculated water surface,
channel invert and crown in absolute SI elevation. Reverse-flow cases are
automatically plotted in their actual flow direction. A labeled topology
overview precedes the profiles and uses the calculated flow signs for its arrows.
Every channel label contains the calculated SI discharge, and node labels list
the flow-oriented endpoint invert elevation and water depth (`z_e`, `h_e`,
`z_v`, `h_v`). Every profile includes a red dashed energy grade line and red
markers for the upstream `z_e+h_e+v_e^2/(2g)` and downstream `z_v+h_v` rest
levels.

All cases are mandatory reference checks; any flow or profile discrepancy above
the documented tolerance makes the command fail. After configuring CMake, the
suite can also be run with:

```bash
ctest --test-dir build -R channel_reference_suite --output-on-failure
```

## Stationary multi-channel network

`channel_network_merge_split.spr` contains six open-surface channel elements.
At `JUNCTION`, two channels arrive and two channels leave; the two downstream
branches then continue through one additional channel each. Four pool elements
prescribe the outer water levels, while STACI solves the three internal water
levels and all channel flows. The two inlet reaches have deliberately unequal
bed slopes: 0.5% in `CHANNEL_1` and 1.5% in `CHANNEL_2`; their source pools keep
the same 0.5 m inlet depth. The `CHANNEL_5` invert rises from 1.00 m to 1.05 m
along the actual `BRANCH_A` to `SINK_A` flow direction, giving a -0.05%
adverse slope; the parallel `CHANNEL_6` retains a conventional falling bed.

`run_channel_network_test.py` requires at least five channel elements, locates
the 2-in/2-out junction from the SPR topology, runs the stationary solver, and
checks convergence, nonzero finite flows, open depths at both ends of every
channel, and the junction mass balance. It additionally verifies that
`CHANNEL_5` carries flow from `BRANCH_A` to `SINK_A` while its bed rises, and
checks that the calculated adverse-slope profile is finite and positive. It
also asserts a low positive slope in `CHANNEL_1`, a steep
slope in `CHANNEL_2`, and at least a 2.5:1 slope ratio. Run it directly or
through CTest:

```bash
python3 tests/run_channel_network_test.py --binary build/staci
ctest --test-dir build -R channel_network_merge_split --output-on-failure
```

The network test retains `channel-network-longitudinal-profile.svg` and a
multi-page `channel-network-longitudinal-profile.pdf`. Both start with a
topology sketch naming all `CHANNEL_*` elements and junctions; the test verifies
that every expected identifier is present. It also reports `Q` on every channel
and all distinct `e`/`v` endpoint elevation-depth pairs at each node. Since the
network branches, the drawing contains a separate panel for each complete
source-to-outlet flow path; shared reaches are repeated so every panel remains
a continuous longitudinal section. Each panel uses the same water-surface,
energy-line and rest-level conventions as the one-channel plots. The renderer has no third-party Python
dependencies and can also be invoked directly:

```bash
python3 tests/plot_channel_profiles.py --network solved-network.spr \
  --output profile.svg --pdf-output profile.pdf
```
