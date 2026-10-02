# PSV, PBV, emitters and pressure-dependent demand

The shared STACI core now supports all six EPANET valve types: TCV, PRV, PSV, PBV, FCV and GPV. PSV/PBV extend the existing `EpanetValve` class, derived from `JelleggorbesFojtas`; they reuse settings, status, controls, EPS and export support.

* PSV holds upstream pressure when feasible, opens when downstream pressure prevents regulation, and closes against reverse flow. Pressure units and specific gravity use the same conversion as PRV.
* PBV imposes a head drop. When the open-valve minor loss exceeds its setting, it uses that loss instead. Explicit OPEN/CLOSED overrides follow EPANET.
* `EpanetEmitter` is a new single-junction outlet element derived from `Agelem`. Its signed discharge is `Q = C sign(p) |p|^exponent`, with coefficient conversion from input flow/pressure units to SI. The inverse relation supplies its Newton equation and derivative. Negative-pressure atmospheric inflow follows EPANET and produces `EPANET.EMITTER_BACKFLOW` in the common diagnostics stream.
* PDA is a junction demand model: positive nominal demands are zero below minimum pressure, proportional to `((p-pmin)/(preq-pmin))^exponent` between minimum/required pressure, and fully delivered above required pressure. Negative demands remain fixed injections. The pressure derivative enters the continuity Jacobian; demand patterns continue to update nominal demand. Results report delivered demand plus emitter discharge.

The solver initializes emitter flow consistently with its pressure guess, refreshes the nonlinear Jacobian, and backtracks pressure-outlet Newton steps using fixed physical merit scales. Running constant-power pumps cannot cross zero flow during these safeguarded steps. Final RMS tolerances and explicit iteration limits are unchanged. Original public INP files are unchanged.

## Verification

`epanet_pressure_elements` adds 54 independent initial-state EPANET comparisons in LPS/GPM/CFS, 54 lossless INP export round trips, and six emitter/PDA EPS report-frame comparisons with varying source pressure. Fixtures cover active/open/closed/reverse/minor-loss valves and negative/zero/partial/full pressure states. Existing malformed-input tests cover unknown junctions, nonfinite/negative emitter coefficients, invalid exponents and incompatible PDA pressure limits.

Seven previously rejected public networks now solve and match EPANET: `Net6_plus`, `epanet_leaks`, `NET1emit`, `NET1negemit`, `psv_open_no_downstream_sources`, `NET1-PBV`, and `cheung`. Targeted tests and the full corpus check initial hydraulics; no full-corpus EPS/quality equivalence is claimed.

## Remaining issues

`io.inp` is no longer rejected for PSV/PDA. Its running POWER pump `pump2` feeds terminal zero-demand junction `j4`. Continuity requires zero flow there, incompatible with finite head at nonzero constant power. `EPANET.POWER_PUMP_DEAD_END` describes this model issue and suggests stopping the pump or providing an outlet. The input is retained. A nominal EPANET convergence result at effectively zero pump flow does not establish a physical operating point.

`JEP5-13.inp` now supports its negative PRV pressure setting and matches the
independent EPANET reference. `ky10`, `NW_Model`, `NW_Model1` and `GES4-9` also
match. These are no longer open compatibility failures.

The strict public initial-state run has 88 matches, three EPANET input rejections
and six unreliable references. Standalone outcomes are 89 solved, four
compatibility/input rejections and four hydraulic failures. The latter are `io`
and the three physical valve conflicts detailed in [valve validation](epanet_valves.md).
No retained network is blocked by a missing PSV/PBV/emitter/PDA implementation.
These counts do not establish full-corpus EPS or chemical equivalence; see the
[current EPS and chemical status](epanet_reference.md).
