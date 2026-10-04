# STACI flushing reference

[Back to the README](../README.md)

Commands below assume the repository root as the working directory unless stated otherwise.

## Single and simultaneous hydrant flushing analysis

`staci_flush` solves a baseline with every added hydrant closed, then opens one
hydrant at a time, restoring the baseline before each independent steady-state
scenario. It lists original pipes where **absolute velocity is strictly greater
than the supplied threshold**, including pipes flowing in reverse. Existing
pump states stay fixed. It ranks hydrants by additional pipe volume and estimates
opening times from directed flow travel times. It does not simulate sediment removal.

```bash
./build/staci_flush \
  --inp /path/to/network.inp \
  --config /path/to/flushing_config.json
```

The default `mode` is `single`, preserving independent one-at-a-time scenarios.
Set `"mode": "multi"` to open all listed junctions in one hydraulic solve.
`MODE` is accepted as an alias; specifying both keys is an error.

Put the flushing settings in `flushing_config.json`:

```json
{
  "mode": "single",
  "hydrant_area_m2": 0.002,
  "total_loss_coefficient": 2.0,
  "velocity_threshold_mps": 0.5,
  "output_dir": "results/flushing",
  "min_pressure_head_m": 0.0,
  "hydrant_node_ids": ["J1", "J2"],
  "write_network_files": true
}
```

`output_dir` is resolved relative to the **config file's directory**, regardless
of the current working directory; absolute paths also work. The pressure setting
is optional and defaults to zero. `write_network_files` defaults to false. Put
exact, unique junction IDs in `hydrant_node_ids`.
Alternatively, omit this array and supply the legacy `--hydrants` input.
Supplying both is rejected. Area, coefficient, velocity and output directory
are required. Unknown keys, wrong JSON types and settings duplicated on the command line are rejected. The original
individual command-line options remain supported when no config is supplied.

See the [worked example](../examples/flushing/README.md) for a complete three-hydrant
network, JSON config, expected plan and step-by-step explanation. Model-specific
projects can provide their own launcher that calls the built executable. Choose
a new or empty output directory for each run.

These numbers are illustrative; provide the outlet characteristics and criteria
for your project. Area is in **m²**, velocity in **m/s**, pressure head in **m**.
The global outlet law is `Q = A * sqrt(2*g*max(h,0)/K)`, with `g = 9.81 m/s²`,
at the junction elevation. **K is the total resistance coefficient, including
outlet kinetic head**: if a supplied local loss coefficient is `zeta`, use
`K = 1 + zeta`. Hydrant discharge is solved together with network pressure;
it is not imposed as a constant withdrawal. Separate outlet elevations and
individual hydrant parameters are not supported yet.

The optional legacy text hydrant list contains one exact junction ID per line; blank lines and
lines beginning with `#` are ignored. Alternatively, use the existing Pusztavacs
JSON format:

```json
{"hydrants": [{"dxf_handle": "H1", "node_id": "J1", "status": "matched"}]}
```

All entries must be matched to existing junctions. Asset IDs must be unique;
distinct assets may share a junction and are evaluated separately in `single`
mode. In `multi` mode, repeated junctions are rejected to avoid ambiguous outlet areas. Without
`dxf_handle`, the node ID is used as the asset ID. Unknown nodes, duplicate
assets and unmatched entries are errors.

The output directory must be new or empty. Inputs are copied there before
loading, so STACI's temporary `.ros`/`.rps` files cannot affect the source model.
Outputs are:

| File | Content |
| --- | --- |
| `run.json` | Parameters, input fingerprints, initial pump states, assumptions and run status |
| `baseline_pipes.csv` | Original pipe geometry, signed flow/velocity and baseline threshold flag |
| `scenarios.csv` | One row per hydrant: status, discharge, pressure, residuals and qualifying pipe count/length/volume |
| `summary.txt` | Human-readable qualifying pipe volume (m³) for each hydrant, plus unique covered volume across valid scenarios |
| `networks.csv`, `networks/<original_name>_hydrant_<nodeID>.inp` | Optional scenario-to-file index and EPANET inspection networks |
| `config.json` | Copy of the supplied configuration |
| `flushing_plan.csv`, `flushing_plan.txt` | Greedy activation sequence, additional/cumulative volume, coverage percentage and both opening-time estimates |
| `pipe_travel_times.csv` | Per-pipe transport diagnostics, in seconds |
| `scenario_hydrants.csv` | Individual active-outlet flow and pressure for each scenario |
| `scenario_pipes.csv` | Every original pipe for each valid scenario, signed/absolute velocity and threshold flags |
| `pipes_above_threshold.csv` | Only qualifying pipe–hydrant pairs from valid scenarios |
| `pipe_coverage.csv` | Every pipe, covering hydrants, maximum valid absolute velocity and corresponding hydrant |
| `work/network.inp`, `hydrants.*` | Copies of the inputs, alongside solver working files |

### Simultaneous opening (`multi`)

The area and loss coefficient apply **to each outlet**, and all outlet flows
are solved together with the network. This is not a sum of single-hydrant runs.
`scenarios.csv` has one combined row, with `hydrant_id` equal to `multi` and
`node_id` containing a JSON array of active junction IDs. Discharge is the sum
of outlet discharges; `hydrant_pressure_head_m` is their minimum pressure head.
`scenario_hydrants.csv` gives each outlet's individual flow and pressure.
`run.json` records the mode and `scenario_count: 1`.

For a valid scenario, `flushing_plan.csv` and the numbered plan in
`flushing_plan.txt` contain one result row. Each qualifying pipe's volume is
counted once, so total, additional and cumulative volumes coincide. The
volume/flow estimate divides this volume by **total** outlet discharge.
The travel-time estimate takes the longest reachable route to any active
outlet. It also considers routes passing through an open junction toward
another outlet: the first outlet may withdraw only part of the incoming flow.
Pipes unable to reach any active outlet through above-threshold pipes are
ignored in the travel-time estimate.

Invalid simultaneous hydraulics produce one diagnostic row in `scenarios.csv`,
exit code 3, and no actionable plan rows, following the existing invalid-scenario
convention. With `write_network_files: true`, one `<original_name>_multi.inp`
contains all open-hydrant emitters. In coverage tables, `multi` identifies the
combined scenario; `hydrant_count` consequently counts covering scenarios.

### Flushing sequence and opening times

Each run also writes `flushing_plan.csv`, `flushing_plan.txt` and
`pipe_travel_times.csv`. The console prints the ordered plan.

In `single` mode, the first hydrant covers the largest qualifying pipe volume. After each
selection, the remaining hydrants are rescored using **only pipe volume not
covered by any earlier selected hydrant**. Exact ties use lexicographic hydrant
ID order. Every hydraulically valid hydrant is ranked; zero-additional-volume
entries appear last and are marked `redundant`. Invalid hydraulic scenarios are
excluded. Coverage here is a bookkeeping assumption, not a measured cleaning state.

For each scenario, flow direction defines a directed graph. A pipe's transit time
is `length / abs(velocity)`. Starting at the upstream end of **every qualifying
pipe**, the code finds the longest route to the active hydrant and takes the
maximum travel time. **Every pipe on the route must satisfy `abs(v) > threshold`.**
A below-threshold (or exactly-at-threshold) connector breaks the route. Remote
qualifying pipes without another fully qualifying route are ignored for timing.
Previously covered pipes are still included when calculating a later hydrant's
opening time. Pump transit is approximated as zero; tank and reservoir nodes
terminate transport routes. Links with `abs(Q) <= 1e-12 m³/s` are treated as stagnant.

The reported value is a **maximum route time** over the remaining paths, not
first arrival along the shortest route or a flow-weighted mean. Slow connectors
are excluded rather than assigned an artificial higher velocity. If another
fully qualifying path exists, that path can still contribute. This rule applies
to both `single` and `multi` modes.

`pipe_travel_times.csv` appends `route_json`: an ordered array of link IDs,
flow-directed `from`/`to` nodes and `travel_time_s` for the controlling route of
each qualifying pipe. Summing those segment times reproduces the row's total.
Use the plan's `critical_pipe_id` to locate the route responsible for a large
time. Ignored pipes have status `ignored_no_qualifying_path`, a blank time and
an empty route array. If none of the qualifying pipes remains connected, the
plan reports zero time with `no_connected_qualifying_pipes`. This is distinct
from a scenario having no qualifying pipes at all (`no_qualifying_pipes`).

The plan includes total, additional and cumulative pipe volumes, followed by
cumulative volume as a percentage of **all original network pipe volume**
(including pipes no hydrant covers). Volume outputs in m³ use two decimal places;
ranking and percentage calculations retain full precision internally. It also includes newly covered
pipe IDs, opening time in minutes (one decimal place), timing status, and the controlling
pipe ID. `pipe_travel_times.csv` gives the underlying per-pipe times. A directed cycle along an otherwise eligible route to the hydrant makes that
hydrant's opening time **undetermined**, left blank in the CSV. The volume rank
is retained. A scenario with no qualifying pipes gets zero time. `run.json`
records the number of undetermined times separately from hydraulic failures.

**Only the route-based time is filtered this way.** Velocity coverage, pipe-volume
ranking and the separate volume/flow estimate continue to include all qualifying
pipes, including those ignored for route timing.

A second timing is reported alongside the travel-time estimate:
`volume_over_flow_time_min = total qualifying pipe volume / hydrant discharge / 60`.
Both times in the plan use minutes with one decimal place. The CSV also includes
`hydrant_flow_m3s` and `volume_over_flow_status`.
The numerator is the full qualifying volume for that individual hydrant, including
previously covered pipes; it is neither the additional volume nor the cumulative
network volume. Calculations use unrounded volumes and discharge in m³/s. Zero
qualifying volume gives zero time; a nonpositive discharge with positive volume
leaves this estimate undetermined. Both estimates remain separate: the program
does not automatically choose between them. Volume/flow time represents one
nominal volume exchange, not a guarantee of contaminant clearance.

Finite times are **advective transport estimates**, assuming the material is
mobilized at the start of a steady scenario. They do not model continued sediment
release, dispersion, storage mixing or changing demands. Flow splitting can send
part of the material toward customer demands or other outlets; arrival at the
hydrant is not a guarantee that all contamination exits there. An undetermined
time requires review before treating the ranked list as an executable schedule.

Pipe volume is the sum of `pi * diameter² / 4 * length` for qualifying original
pipes, in m³; it is reported on the console, in `summary.txt`, and as
`qualifying_volume_m3` in `scenarios.csv`. This is contained pipe volume, not the
volume discharged through a hydrant over time. The unique coverage total counts
each pipe once across valid scenarios. Invalid scenarios have no qualifying volume.

Set `write_network_files` to `true` to export one INP per hydrant (including
clearly marked invalid scenarios). `networks.csv` maps filenames to node IDs. Files use
`<original_name>_hydrant_<nodeID>.inp`, for example
`Pusztavacs_production_hydrant_J184.inp`. Open a file in EPANET and inspect a pipe's **Tag**:
`FLUSHED`, `BELOW_THRESHOLD`, or `INVALID_SCENARIO`. These tags replace original
pipe tags in exported copies only; original tags remain in the input model.
The original coordinates, vertices and other drawing sections are retained.
Tags are static STACI classifications; EPANET does not automatically colour them.
To highlight velocity coverage in EPANET, run the analysis and select
**View → Query → Links → Velocity → Above**, then enter your threshold. This
shows EPANET's recalculated result, which can differ slightly from STACI.
Double-click a pipe to read its saved Tag. Re-running EPANET does not update tags.

The active hydrant is represented by an emitter with exponent 0.5, using the
same area/K law and the input model's flow/pressure units. See the
[EPANET emitter definition](https://wateranalytics.org/EPANET/_emits_page.html).
Exports freeze customer demands and reservoir heads at the evaluated snapshot,
remove time patterns and controls/rules, and set duration to zero. Pump statuses
and initial tank levels remain as supplied. Hydraulic results recalculated by
EPANET can differ from STACI; the tags always describe the STACI result.
These inspection exports contain emitters and cannot be fed back to the current
`staci_flush` importer, which rejects pre-existing emitters.

`newly_above_threshold` distinguishes pipes that cross the threshold only after
opening the hydrant from those already above it in the baseline. Coverage includes
both. Flows are in m³/s, with an additional L/s column in the scenario summary.
The coverage CSV stores hydrant lists as JSON arrays inside quoted CSV cells.

`--min-pressure-head-m` optionally sets a service criterion at **every junction**;
it defaults to zero. Negative-pressure, unconverged, nonfinite, or excessive
residual results are excluded from pipe coverage. Results below a positive
service criterion remain diagnostic rows in `scenarios.csv`. The baseline must
converge with nonnegative junction pressures, but need not meet that optional
service criterion. Customer demands remain pressure-independent. A numerically
valid result alone does not demonstrate acceptable service pressure.

Only `--snapshot initial` is currently supported (and is the default): initial
tank levels are fixed, first demand/reservoir pattern multipliers are used, and
rules/controls are frozen. OPEN/CLOSED initial pump statuses are applied explicitly.
No network valve is operated. This first flushing implementation rejects
INP files containing `[VALVES]`, `[EMITTERS]`, closed
or check-valve pipes, nonzero pipe minor losses, pump speed/pattern settings,
or pressure-dependent demands are rejected. H-W and D-W headloss are supported
with positive roughness; non-default viscosity is rejected. Models requiring
these features need additional flushing-specific support before use here.

Exit codes: `0` all scenarios valid; `1` calculation/baseline failure;
`2` invalid input/configuration; `3` one or more invalid scenarios
(valid scenarios are still exported). Run
`staci_flush --help` for the argument list.

Implementation lives in `src/staci_flush.cpp`, `src/flushing.cpp`,
`src/flushing_plan.cpp` and `src/HydrantOutlet.cpp`, sharing the internal `staci_core` solver target with
`staci`. Model-specific folders only need their inputs and a command or launcher
that calls this executable.

The CTest suite includes the worked example, outlet-law/Jacobian checks and
end-to-end tests for analytic flow, reverse-flow coverage, scenario independence, JSON input, pump
status, invalid inputs, pressure filtering and source-file preservation. Tests
run automatically on pull requests and pushes to `master` on Linux, macOS and Windows:

```bash
ctest --test-dir build --output-on-failure
# Run only the standard flushing regressions and worked example:
ctest --test-dir build -R flushing --output-on-failure
```
