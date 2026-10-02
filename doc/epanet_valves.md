# EPANET valves

The shared hydraulic core now imports PRV, FCV and GPV valves in addition to TCV. All four applications use this core. PSV and PBV are now supported too; emitter/PDA support is described in [the pressure-element update](epanet_pressure_elements.md).

* PRV fixes downstream pressure head when regulation is feasible, opens when upstream head is insufficient, and closes against reverse flow. Its setting is converted from the input pressure units, including specific gravity.
* FCV fixes the requested forward volume flow when feasible and behaves as an open valve when pressure cannot sustain it. `EPANET.FCV_UNATTAINABLE` explains the requested and achieved flows in the common diagnostics stream.
* GPV uses its named flow/headloss curve, with piecewise linear interpolation and endpoint extrapolation. Flat headloss segments are supported. Missing, short, nonfinite, negative or decreasing curves receive a section/element/line compatibility diagnostic; flow coordinates must strictly increase.

Fixed OPEN/CLOSED states, numeric PRV/FCV settings and EPS controls use the existing valve control interface. Numeric GPV settings are rejected because its setting names a curve. Automatic hydraulic closure remains distinct from external disabling, so tank-limit controls cannot be undone by valve status updates. Imported settings and GPV curves are included in SI INP export; automatic initial guesses are not exported as forced OPEN states.

The solver rebuilds the Jacobian for valve networks, checks statuses before declaring convergence and backtracks Newton steps. Tiny coupling/closed-link resistance prevents floating or competing regulators from immediately producing a singular Jacobian. Limits remain configurable with the common CLI options. Native SPR pipe/throttle formulas are preserved.

Reference testing also exposed and corrected the existing EPANET TCV minor-loss conversion and HEAD-pump reverse-flow closure. EPANET's published 0.02517 valve coefficient is converted to SI consistently with imported pipes. Reservoir/tank net exchange is compared at the same 1e-6 m3/s absolute resolution as link flow; prescribed junction demand retains 1e-8 m3/s. This distinction accounts for EPANET's closed-link penalty leakage and is common to all networks, not tuned per case.

## Current verification (2026-10-02)

`epanet_reference_valves` covers 62 independent steady cases and nine controlled
EPS snapshots, including all ten flow units, active/open/reverse/closed states,
minor losses and export round trips. `epanet_pressure_elements` extends this to
PSV/PBV/emitter/PDA. The current public initial-state reference run has 88 matches;
ky10 and JEP5-13 are included. Negative PRV/PSV pressure settings are supported.

## Infeasible regulation and physical diagnostics

The three formerly labelled strict solver failures have physical constraint
conflicts in the computed regularized state:

| Network | Conflict |
|---|---|
| `fcv_open_no_downstream_sources` | FCV `VALVE` is set to 25 GPM while downstream continuity needs 250 GPM; the alternative tank connection is closed |
| `prv_closed_no_upstream_sources` | The upstream pump is closed; supplying upstream demand requires about 250 GPM of reverse flow through the closed PRV |
| `2fcvs` | Active FCV `44` is set to 6.5 GPM while the computed continuity solution requires 17 GPM |

Small diagonal/coupling terms help the Newton matrix but can mathematically allow
these conflicts at enormous pressure differences. `HYDRAULICS.VALVE_CONSTRAINT`
now reports the valve, its setpoint/state and computed flow, and rejects this
state even if the RMS residual meets the normal 0.1 mm limit. The independent
EPANET references also fail physical acceptance checks. These models remain
unchanged diagnostic fixtures; success would require a feasible supply/demand
configuration or a deliberate change to valve operation.

`epanet_infeasible_valve_diagnostics` verifies the three cases and `io` at default
and strict tolerances, preserving input hashes and requiring failure diagnostics.
This is error-handling coverage, not a numerical-equivalence claim.

Algorithms follow the pinned EPANET 2.2 hydraulic-status and hydraulic-coefficient
sources. See [current reference status](epanet_reference.md) for the validation
scope, tolerances and reproduction commands.
