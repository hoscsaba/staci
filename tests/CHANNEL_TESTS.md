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
for automated comparison.

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
levels and all channel flows.

`run_channel_network_test.py` requires at least five channel elements, locates
the 2-in/2-out junction from the SPR topology, runs the stationary solver, and
checks convergence, nonzero finite flows, open depths at both ends of every
channel, and the junction mass balance. Run it directly or through CTest:

```bash
python3 tests/run_channel_network_test.py --binary build/staci
ctest --test-dir build -R channel_network_merge_split --output-on-failure
```
