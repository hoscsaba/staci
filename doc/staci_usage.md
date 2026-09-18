# STACI usage reference

[Back to the README](../README.md)

Commands below assume the repository root as the working directory unless stated otherwise.

## Quick start

The commands below use the included `tests/anytown_1med.spr` network. From the
repository root, use `./build/staci` on Linux or macOS and
`.\build\Release\staci.exe` on Windows.

### Import an EPANET network

STACI recognizes the `.inp` extension automatically. For example:

```bash
./build/staci -l model.inp
./build/staci -g model.inp -e J-1 -p demand
./build/staci -s model.inp
```

The importer reads every EPANET section and reports a warning for data that
does not have a STACI equivalent. It currently transfers:

- junction elevations, separate demand categories, category labels, pattern
  references and complete demand-pattern series, demand multipliers, and
  per-node initial-quality metadata; the steady snapshot applies the first
  multiplier while retaining the remaining values;
- reservoir base heads, complete head-pattern references and multiplier series,
  and tank initial levels as fixed hydraulic boundaries;
- pipes, with EPANET-to-SI conversion for flow, length, diameter, and
  Darcy-Weisbach roughness, dimensionless minor-loss coefficients,
  `Open`/`Closed` status overrides, and one-way `CV` operation;
- Hazen-Williams and Darcy-Weisbach head-loss selection;
- pump head curves and constant-power pumps, including initial status,
  relative speed, speed patterns, original curve IDs, and retained energy and
  efficiency references;
- throttle control valves (`TCV`) as constant, 50%-position
  `JelleggorbesFojtas` elements. The importer retains the dimensionless TCV
  setting, diameter, separate minor-loss coefficient, and
  `ACTIVE`/`OPEN`/`CLOSED` state while converting the active hydraulic
  coefficient to STACI's SI mass-flow form;
- relevant solver settings such as trial count and accuracy (the latter is an
  approximate mapping because the convergence criteria differ);
- a lossless EPANET document containing the original text, section and record
  order, line endings, inline comments and their record IDs, tags, map data,
  and non-hydraulic run configuration.

Rule-based controls are executed in EPANET EPS mode but do not affect a single
steady-state import. TCV settings and states are supported; other regulated
valve types and emitters are not simulated. Chemical EPS applies initial
quality, all four EPANET source types and their patterns, plus first-order global
and pipe-specific bulk and wall reactions. Tank chemical storage/mixing and
non-first-order reactions remain unsupported and produce explicit warnings.
Descriptive metadata, map sections, and non-hydraulic run configuration do not
affect the hydraulic snapshot, but they remain attached to the imported
EPANET document for lossless export. Each incompatible calculation feature is
identified on standard error with its section, element where available, and
source line.

Hydraulic results are not written into `.inp` files because EPANET input files
describe network definitions rather than simulation results. Hydraulic runs
can still create the normal `.ros`, `.rps`, and `.rrs` sidecar files. Network
properties can be changed safely with `-m`; STACI writes a new `.inp` file and
does not alter the source.

### Export a STACI network to EPANET

```bash
./build/staci --export-epanet network.spr -o network.inp
# Equivalent short form:
./build/staci -y network.spr -o network.inp
```

The generated file uses metric `LPS` units. Junctions, fixed-head boundaries,
tanks, pipes, pump curves, constant-power pumps, constant non-negative throttle
valve curves, demands, initial chemical quality, and the selected pipe-friction
model are exported. A compatible throttle curve is written as an EPANET `TCV`.
Tank operating limits are inferred because native STACI pools do not store all
EPANET tank constraints. STACI-only elements are skipped with explicit
warnings.

When the source is an EPANET INP file, STACI exports its retained
`EpanetDocument` representation. `[TITLE]`, inline and full-line comments,
`[TAGS]`, `[COORDINATES]`, `[VERTICES]`, `[LABELS]`, `[BACKDROP]`, `[TIMES]`,
`[REPORT]`, `[ENERGY]`, and all `[OPTIONS]` records retain their original text,
association, order, whitespace, and line endings. An unchanged INP therefore
round-trips byte-for-byte. The `-m` INP modification path changes only the
requested field and preserves all other lines.

### EPANET–STACI element mapping (easiest first)

The table below compares the hydraulic node and link types in EPANET with the
current STACI model. The order estimates the effort required to achieve a
complete, editable and hydraulically meaningful mapping, starting with the
closest matches. `Direct` means that a corresponding STACI type exists and is
already created by the INP importer; it does not imply that every auxiliary
EPANET field is supported yet.

| Order | EPANET element | Closest STACI representation | Current mapping | Work needed for complete compatibility |
|---:|---|---|---|---|
| 1 | Junction | `Csomopont` | **Direct; structured metadata retained** | |
| 2 | Reservoir | `Csomopont` + `KonstNyomas` | **Direct; complete head pattern retained and applied** | |
| 3 | Pipe | `Cso`, including one-way `CV` mode | **Direct; hydraulic and status fields retained and applied** | |
| 4 | Pump, `POWER` form | `EpanetPowerPump` | **Direct; definition, status, speed, patterns and energy metadata retained** | |
| 5 | Pump, `HEAD` curve form | `Szivattyu` | **Direct; original curve and all pump attributes retained** | |
| 6 | Tank | `Csomopont` + `Vegakna` | **Supported in EPS; approximate in steady mode** | EPS applies operating limits, cylindrical storage, and referenced piecewise-linear volume curves. The steady importer uses initial level as a fixed boundary; native SPR export must infer missing operating limits. Tank quality mixing remains. |
| 7 | Throttle control valve (`TCV`) | `JelleggorbesFojtas` | **Direct; constant curve, setting, minor loss and status retained and applied** | |
| 8 | General-purpose valve (`GPV`) | `JelleggorbesFojtas` | **Similar STACI type exists; not mapped** | Adapt the referenced GPV head-loss curve without losing its ID or points and validate that its sign and interpolation semantics match STACI. |
| 9 | Emitter attached to a junction | None | **No equivalent** | Add a pressure-dependent outlet representation with emitter coefficient and pressure exponent, plus INP import/export support. |
| 10 | Pressure breaker valve (`PBV`) | None | **No exact equivalent** | Add a valve with a prescribed pressure drop and open/closed state; an ordinary throttling curve is not an exact substitute. |
| 11 | Flow control valve (`FCV`) | None | **No equivalent** | Add a flow-setpoint valve and its active/open/closed operating states. |
| 12 | Pressure reducing valve (`PRV`) | None | **No equivalent** | Add downstream-pressure regulation and the EPANET active/open/closed state model. |
| 13 | Pressure sustaining valve (`PSV`) | None | **No equivalent** | Add upstream-pressure regulation and the EPANET active/open/closed state model. |

EPANET also has data objects that are not standalone hydraulic elements. Their
current STACI counterparts, ordered approximately from easier to harder to
represent structurally, are:

| EPANET data object | STACI support |
|---|---|
| Base demands and `[DEMANDS]` categories | **Supported:** retained as separate typed node metadata with category labels, independent pattern references and complete multiplier series. The hydraulic snapshot receives their time-zero aggregate. |
| Curves | **Partial:** HEAD pump curves and pump-efficiency curves are retained with their original IDs and SI points. Referenced tank-volume curves are applied by EPS using EPANET piecewise-linear volume-depth interpolation. A general editable tank-curve model and GPV curve mapping remain. |
| Patterns | **Partial:** complete demand, reservoir-head, pump-speed and pump energy-price patterns are retained by their corresponding STACI elements; hydraulic patterns are used by EPS. A reusable network-level model for other pattern types is still missing. |
| Simple controls and rule-based controls | **Supported in EPS:** simple controls execute `OPEN`, `CLOSED`, `ACTIVE`, numeric pump speeds and numeric TCV loss settings for node pressure/tank level, elapsed-time, and clock-time triggers. The EPS rule engine supports `IF`/`AND`/`OR`, `THEN`, multiple actions, `ELSE`, priorities, node/tank/link/system premises, TCV setting/status premises, and the separate rule timestep. Both sections remain losslessly retained in the imported INP document; neither yet has a network-level editable STACI object outside EPS. |
| Water age, chemical sources, reactions and tank mixing | **Partial:** EPS solves `QUALITY AGE` and `QUALITY CHEMICAL` with segment-based plug flow, timestep-volume node mixing, reservoirs, pumps and valves. Initial quality, `CONCEN`/`MASS`/`SETPOINT`/`FLOWPACED` sources with patterns, and first-order global or pipe-specific bulk/wall coefficients are applied. Tank storage/mixing, trace mode, and non-first-order kinetics remain. |

STACI also contains `Csatorna` (open channel) and `BukoMutargy` (overflow/weir),
for which EPANET has no native hydraulic element. These elements are therefore
skipped with a warning when a native SPR network is exported to INP.

An imported TCV and any native `JelleggorbesFojtas` with a constant,
non-negative loss curve are exported as EPANET `TCV` links without a warning.
The TCV setting is converted through `h_L = K v^2/(2g)` rather than copied into
STACI's mass-flow coefficient. A genuinely position-dependent native STACI
curve still cannot fit into one TCV setting; it is omitted with a warning that
identifies the missing connection.

TCV properties can be inspected or changed explicitly:

```bash
./build/staci -g model.inp -e TCV1 -p tcv_setting
./build/staci -m model.inp -e TCV1 -p tcv_setting -n 8 -o modified.inp
```

The mapping status above describes the editable STACI calculation model. It is
separate from lossless INP preservation: unsupported records and sections in an
imported `EpanetDocument` can still be re-exported unchanged even when STACI
cannot simulate them.

### EPANET–STACI compatibility TODO (easiest first)

The following work is ordered by expected implementation difficulty. It is
limited to INP data coverage, lossless import/export, and the set of available
network elements; changes to the hydraulic solution method are intentionally
out of scope. Each completed item should include a minimal INP fixture and a
field-level `INP -> STACI -> INP` round-trip assertion.

1. [x] **Preserve descriptive metadata.** Store and re-export `[TITLE]`, inline
   comments, and `[TAGS]` without changing their text or association with
   nodes and links.
2. [x] **Preserve map metadata.** Round-trip `[COORDINATES]`, `[VERTICES]`,
   `[LABELS]`, and `[BACKDROP]`, including numeric values, backdrop declarations,
   associations, and ordering.
3. [x] **Preserve non-hydraulic run configuration.** Store and re-export all
   recognized `[TIMES]`, `[REPORT]`, `[ENERGY]`, and currently unused
   `[OPTIONS]` entries instead of reducing them to warnings.
4. [ ] **Retain the source unit system.** When an INP network is imported and
   exported again, preserve its original EPANET flow units and convert every
   affected field back consistently; keep `LPS` as the default for native SPR
   exports.
5. [x] **Complete pipe field coverage.** Import and export minor-loss
   coefficients, `Open`/`Closed`/`CV` pipe status, and matching `[STATUS]`
   overrides. `CV` uses a one-way `Cso` mode so the original pipe geometry,
   friction model, and minor-loss coefficient remain active.
6. [ ] **Round-trip all tank fields.** Retain minimum and maximum level,
   minimum volume, diameter, initial level, and the referenced volume-curve ID
   instead of inferring missing values during export.
7. [x] **Preserve multiple demand categories.** Keep every `[DEMANDS]` row,
   category label, pattern reference, and junction base demand separately
   instead of storing only one summed STACI demand value. The solver-facing
   aggregate remains available without replacing the retained components.
8. [ ] **Add a native pattern data model.** Store complete `[PATTERNS]` series,
   the default pattern, and all object-to-pattern references so an
   INP–SPR–INP round trip does not collapse a pattern to its first multiplier.
9. [x] **Complete pump data coverage.** Preserve pump curve IDs, `POWER`/`HEAD`
   form, initial status, relative speed, speed-pattern reference, and associated
   efficiency and energy records during import and export.
10. [ ] **Classify and preserve every curve.** Distinguish pump, efficiency,
    volume, and general-purpose curves, validate their references, and export
    them with their original IDs and numeric data.
11. [ ] **Represent emitters as network elements.** Add an import/export mapping
    for `[EMITTERS]`, including the emitter coefficient and the relevant
    pressure-exponent option, with an explicit placeholder representation when
    no equivalent STACI element is available.
12. [ ] **Map general-purpose valves.** TCV records are now mapped to constant
    `JelleggorbesFojtas` curves with unit-correct setting, minor-loss and status
    handling. Implement the remaining EPANET `GPV` mapping using its referenced
    head-loss curve without changing its interpolation semantics.
13. [ ] **Represent regulated EPANET valves.** Add native element records and
    lossless INP mappings for `PRV`, `PSV`, `PBV`, and `FCV`, including setting,
    minor-loss coefficient, and initial status.
14. [ ] **Expose controls and rules as editable network objects.** EPS now
    parses and executes `[CONTROLS]` and `[RULES]` as typed runtime objects, and
    imported INP text is re-exported losslessly. A reusable network-level API
    is still needed to edit their references, thresholds, actions, priorities,
    and ordering outside an EPS run.
15. [x] **Preserve water-quality configuration.** Round-trip `[QUALITY]`,
    `[SOURCES]`, `[REACTIONS]`, and `[MIXING]`, including units, source types,
    coefficients, tank mixing models, and pattern references. EPS now applies
    chemical initial quality, sources and first-order reactions; unsupported
    tank mixing and reaction orders remain losslessly preserved and warned.
16. [ ] **Introduce a lossless fallback for unknown sections.** Keep unknown or
    newer EPANET sections and unsupported records as raw, ordered INP data so
    STACI can modify known properties without silently discarding future
    EPANET extensions.

### Run an EPANET extended-period simulation

```bash
./build/staci --epanet-eps network.inp -o results/network
# Equivalent short form:
./build/staci -z network.inp -o results/network
```

EPS mode reads `Duration`, hydraulic/pattern/report timesteps, pattern start,
junction demand patterns, demand multiplier, reservoir head patterns, pump
speed patterns, tank geometry and operating levels, initial link status, and
EPANET simple controls from the input file. Supported controls can open or close
imported pipes and pumps, or assign a numeric relative speed to an imported
`POWER` or `HEAD` pump. Their triggers can use tank level, junction pressure,
elapsed simulation time (`AT TIME`), or time of day (`AT CLOCKTIME`), including
`START CLOCKTIME` and daily clock-event recurrence. Pressure thresholds are
converted from the declared EPANET pressure units to SI pressure head.

Both canonical `HH:MM` and WNTR-generated `HH:MM:SS` time values are accepted.
`START CLOCKTIME 00:00:00 AM` is interpreted as midnight as well as the
canonical `12:00 AM` spelling. Time-based events shorten the current hydraulic step so off-grid event times
are solved exactly without shifting later pattern and report boundaries.
Node-based controls are re-evaluated after the hydraulic solve and the state is
resolved until the control settings stabilize. These runtime control records
belong to the EPS orchestration layer; no new hydraulic `Agelem` subclass is
needed. Pump affinity laws are applied at every hydraulic state. Tank levels
are advanced from the solved boundary flow in volume space. Cylindrical tanks
use their area, while tanks referencing a `[CURVES]` volume curve use
piecewise-linear volume-depth interpolation and inverse interpolation. Flow is
blocked in the prohibited direction at the minimum or maximum level.

The `[RULES]` engine uses the separately configurable `RULE TIMESTEP` (default:
one tenth of the hydraulic timestep). It supports node demand/head/pressure,
tank level/fill time/drain time, absolute link flow/status/setting, total system
demand, elapsed time, and clock time. EPANET's OR-before-AND grouping,
interval-aware time equality, `THEN`/`ELSE` action lists, strict priority
replacement, and first-rule precedence for equal priorities are retained.
Rule actions accept both EPANET's `IS` spelling and WNTR's equivalent `=`
spelling for `STATUS` and `SETTING` assignments. Tank levels are projected with
the last solved flow while the rule interval is scanned. Volume-curve tanks are
projected in volume space, and their fill/drain times use the actual remaining
curve volume. The first rule check that causes a real link change becomes the
next hydraulic event.

The primary output is `PREFIX.h5` using the `STACI EPS OUTPUT v1` schema. A
small `PREFIX.meta.json` sidecar contains run metadata and precomputed value
ranges for plot color bars. Four analysis-friendly SI CSV files are retained
for interoperability and manual inspection:

- `PREFIX-nodes.csv`: pressure head, total head, elevation, demand, water age
  in seconds, and chemical concentration in `kg/m3` for each node and report time;
- `PREFIX-links.csv`: flow, head loss, endpoints, type, enabled state,
  volume-averaged water age, and chemical concentration for each link;
- `PREFIX-tanks.csv`: tank level, volume, inflow, and limits;
- `PREFIX-summary.csv`: duration, timestep, state count, failed states, and
  warning count.

All physical quantities in the HDF5, JSON, and CSV outputs use SI units:
seconds, metres, cubic metres, cubic metres per second, and metres per second.
The CSV files use one observation per row and can be read directly by Excel,
LibreOffice, pandas, R, MATLAB, or plotting tools.

The HDF5 dynamic arrays use `[time, element]` order:

```text
/time                       int64   [T]    s
/nodes/head                 float64 [T,N]  m
/nodes/pressure_head        float64 [T,N]  m
/nodes/demand               float64 [T,N]  m3/s
/nodes/water_age            float64 [T,N]  s
/nodes/chlorine             float64 [T,N]  kg/m3
/links/flow_rate            float64 [T,L]  m3/s
/links/velocity             float64 [T,L]  m/s
/links/headloss             float64 [T,L]  m
/links/status               uint8   [T,L]
/links/water_age            float64 [T,L]  s
/links/chlorine             float64 [T,L]  kg/m3
/tanks/level                float64 [T,K]  m
/tanks/volume               float64 [T,K]  m3
/tanks/inflow               float64 [T,K]  m3/s
/simulation/converged       uint8   [T]
/simulation/iterations      uint32  [T]
```

Static `nodes`, `links`, and `tanks` datasets provide ID-to-index mapping,
topology, elevations, coordinates, lengths, diameters, and operating limits.
Coordinates from an EPANET file are converted using that network's length-unit
family. Missing coordinates and non-applicable link dimensions are stored as
`NaN`.

Dynamic datasets are extensible along the time dimension and are written in
16-frame chunks. For 10,000 nodes a `float64 [16,10000]` chunk is about 1.28 MB.
Shuffle plus gzip level 1 is used when the linked HDF5 library provides the
deflate filter. The file is flushed after every chunk and its `frames_written`
and `simulation_status` attributes allow interrupted runs to be recognized.

If HDF5 was not found while configuring the build, CMake prints a warning and
the EPS command writes only the SI CSV and JSON files. Set `HDF5_ROOT` for a
non-standard installation, for example:

```bash
cmake -S . -B build -DHDF5_ROOT=/path/to/hdf5
```

`QUALITY AGE` and `QUALITY CHEMICAL` are supported for EPS pipe networks with node
mixing, reservoir boundaries, zero-volume pumps/valves, and plug-flow pipe
storage. Chemical concentrations are converted from the INP declaration and
written as SI `kg/m3`; bulk coefficients are converted from 1/day to 1/s and
wall coefficients from length/day to m/s. Tank quality storage/mixing,
non-first-order kinetics, and trace simulation are not yet represented. TCV
`ACTIVE`/`OPEN`/`CLOSED` actions and numeric loss settings are supported;
regulated non-TCV valve types and exact interpolation of a
simple node-pressure or tank-level threshold crossing inside a hydraulic
timestep are not yet implemented. Rule conditions are sampled at `RULE
TIMESTEP` instants, matching EPANET's discrete rule-engine model; unsupported
valve operations emit warnings.

### Steady water-quality solution

For constant demands, link states, reservoir concentrations and reaction
coefficients, STACI can solve the asymptotic water-quality state directly:

```bash
./build/staci --steady-quality network.inp -o results/steady
# Equivalent short form, overriding the INP quality mode:
./build/staci -q network.inp -o results/steady --quality-mode both
```

The direct solver first computes one steady hydraulic state. It directs every
link according to its solved signed flow, uses `tau = volume / abs(flow)` for
pipe travel time, and applies the exact first-order transfer factor
`exp(k*tau)`. Complete instantaneous mixing gives one sparse linear system for
node water age and another for node chemical concentration. Pumps and valves
are represented as zero-volume links. Tanks are fixed quality boundaries at
their imported initial values for this snapshot. This is an asymptotic,
fixed-boundary calculation; EPS must be used when demands, controls, tank levels, source
strengths or flow directions change with time. Non-first-order reactions are
rejected explicitly.

All output values use SI units. The command writes:

- `PREFIX-steady-nodes.csv`: node water age in seconds and concentration in
  `kg/m3`;
- `PREFIX-steady-links.csv`: signed flow in `m3/s`, travel time, volume-average
  age and concentration, and the first-order transfer factor;
- `PREFIX-steady-summary.csv`: input, selected mode and matrix dimensions.

The selected hydraulic parameter can be differentiated without rerunning the
quality solution:

```bash
./build/staci -q network.inp -o results/steady \
  --quality-mode both --quality-sensitivity \
  -e P1 -p diameter
```

`diameter`, `friction_coeff`, and junction `demand` are supported. The existing
hydraulic implicit sensitivity supplies all link-flow derivatives; STACI then
differentiates travel time, wall reaction, pipe transfer and node mixing and
solves the resulting right-hand side with the already factorized quality
matrix. `PREFIX-steady-sensitivity.csv` reports derivatives per SI parameter:
per metre for diameter, per dimensionless roughness value for
`friction_coeff`, and per `m3/s` for demand.

### Native C++ examples

User-facing, MATLAB-independent programs are provided under
[`examples/cpp`](../examples/cpp/README.md). They use the same C++17 STACI sources
as the command-line tools and are built by default; they can be disabled with
`-DSTACI_BUILD_CPP_EXAMPLES=OFF`.

```bash
cmake -S . -B build -DSTACI_BUILD_CPP_EXAMPLES=ON
cmake --build build --target \
  staci_example_hydraulics \
  staci_example_epanet_inp --parallel

./build/staci_example_hydraulics
./build/staci_example_epanet_inp tests/epanet_eps_smoke.inp
```

`staci_example_hydraulics` demonstrates native SPR loading, steady hydraulics
and SI node/link results. `staci_example_epanet_inp` focuses on EPANET INP
import, reports the mapped STACI element types, solves the imported snapshot and
prints the physical EPANET links represented by STACI.

When the official OWA EPANET 2.2 toolkit is available, CMake additionally builds
`staci_example_epanet_eps_compare`. It runs the official toolkit and STACI EPS
on the same INP file, matches time/element identifiers and compares total head,
junction pressure/demand, link flow and velocity in SI units:

```bash
python3 tests/setup_epanet_reference.py
cmake -S . -B build \
  -DSTACI_BUILD_CPP_EXAMPLES=ON \
  -DSTACI_EPANET_TOOLKIT_ROOT="$PWD/build/epanet-reference"
cmake --build build --target staci_example_epanet_eps_compare --parallel

./build/staci_example_epanet_eps_compare \
  tests/epanet_eps_smoke.inp cpp-eps-results
```

The EPS example calls the EPANET toolkit directly rather than parsing its text
report. It prints maximum absolute differences and returns a failing exit code
when the included benchmark tolerances are exceeded. STACI HDF5, JSON and SI
CSV outputs and EPANET's diagnostic report are retained in the selected output
directory.

### General in-memory MATLAB interface

The optional `StaciModel` API keeps each parsed STACI network behind a locked
`staci_mex` handle. After the initial SPR or INP read, property updates,
hydraulic calculations, steady-quality calculations, sensitivities and result
transfers happen directly in memory without temporary data or result files.
Multiple independent `StaciModel` objects may be open in one MATLAB process.
All public numerical properties use SI units.

```matlab
cd('/path/to/staci')
addpath('matlab')

network = StaciModel('tests/anytown_1med.spr');
pipes = network.linkIds("pipe");
diameters = network.getLinkProperty(pipes, "diameter_m");
network.setLinkProperty(pipes, "diameter_m", 1.05 * diameters);

status = network.solveHydraulics();
nodes = network.nodeTable();
links = network.linkTable();
quality = network.solveSteadyWaterAge();
sensitivity = network.hydraulicSensitivity(pipes(1), "diameter_m");
```

Important operations are:

- network lifecycle: constructor, `release`, and `resetHydraulicState`;
- introspection: `nodeIds`, `linkIds`, `linkTypes` and the immutable `Info`
  structure;
- vectorized properties: `getNodeProperty`, `setNodeProperty`,
  `getLinkProperty`, and `setLinkProperty`;
- hydraulics and results: `solveHydraulics`, `nodeResults`, `linkResults`,
  `nodeTable`, and `linkTable`;
- steady quality: `solveSteadyWaterAge` or `solveSteadyQuality` with `age`,
  `chemical`, `chlorine`, or `both` mode;
- implicit hydraulic sensitivity: `hydraulicSensitivity` for `diameter_m`,
  `friction_coeff`, and `demand_m3s` parameters.

Common node properties include `elevation_m`, `pressure_head_m`, `total_head_m`,
`demand_m3s`, `water_age_s`, `concentration_kgm3`, and
`source_concentration_kgm3`. Common link properties include `flow_rate_m3s`,
`velocity_ms`, `status`, and pipe-specific `diameter_m`, `length_m`,
`friction_coeff`, and `minor_loss`. Pump speed/power, valve position/setting,
pool water level and fixed boundary head are also exposed when supported by the
selected STACI element type. Unsupported property/type combinations raise a
MATLAB exception instead of entering a legacy interactive error path.

On first use, `build_staci_mex.m` configures `build-matlab/` and builds the
platform-specific module. It recognizes `CMAKE_COMMAND` and the usual
Homebrew/CMake.app locations when MATLAB does not inherit the login-shell PATH:

```matlab
addpath('matlab')
mexFile = build_staci_mex();
```

User-facing programs are under [`examples/`](../examples/README.md), separately
from developer regression tests. Start with:

```matlab
addpath('matlab')
addpath('examples/matlab')
getting_started_hydraulics();
modify_network_in_memory();
steady_water_quality_demo();
hydraulic_sensitivity_demo();
```

The Anytown Global Optimization Toolbox demonstration is now an example built
entirely on the same general API:

```matlab
result = optimize_anytown_residence_time();
```

It optimizes all 41 pipe diameters continuously between `0.01` and `1.00 m`,
minimizes demand-weighted steady water age and enforces `0...60 m` nodal
pressure head. Feasible seed networks prevent GA from stopping after one
generation when a uniformly random population violates the hydraulic limits.
Its final CSV, MAT and grouped diameter/pressure plot are written to
`matlab/matlab-results/`; GA evaluations remain file-free.

MATLAB developer tests are deliberately kept under `tests/matlab`:

```matlab
addpath('matlab')
addpath('tests/matlab')
results = run_matlab_tests();
```

They cover handle lifecycle, introspection, vector property round trips,
hydraulics, steady water age, central-difference sensitivity validation and all
short user examples. The optimization is exercised separately with reduced GA
settings in development because its default run has 50 generations.

For manual CMake configuration:

```bash
cmake -S . -B build-matlab \
  -DSTACI_BUILD_MATLAB_MEX=ON \
  -DSTACI_BUILD_OPTIMIZERS=OFF \
  -DSTACI_ENABLE_HDF5=OFF \
  -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build build-matlab --target staci_mex --parallel
```

The normal executables continue to use UMFPACK. The MEX target uses Eigen
`SparseLU`, avoiding a second OpenMP runtime inside MATLAB when a package-manager
SuiteSparse build was linked with OpenMP.

### List network elements

```bash
./build/staci -l tests/anytown_1med.spr
```

This also writes `element_list.txt` in the current working directory.

### Run a steady-state hydraulic simulation

```bash
./build/staci -s tests/anytown_1med.spr
```

The solved values are written back to the input network file. Work on a copy if
the original input must remain unchanged:

```bash
cp tests/anytown_1med.spr network.spr
./build/staci -s network.spr
```

`network.spr.rrs` contains `OK` when the solver converges or `ERROR!` when it
does not.

An earlier result file can be used for initialization:

```bash
./build/staci -s network.spr -i previous-result.spr
```

### Read a property

```bash
./build/staci -g network.spr -e NODE210 -p pressure
./build/staci -g network.spr -e PIPE74 -p mass_flow_rate
```

The value is printed and also written to `tmp.dat`. Readable edge properties
are `diameter`, `mass_flow_rate`, `bottom_level`, `water_level`, and `position`.
Readable node properties are `pressure`, `head`, and `demand`.

### Modify a property and save a new network

```bash
./build/staci -m network.spr \
  -e PIPE74 -p diameter -n 0.8 -o network-modified.spr
```

Writable edge properties are `diameter`, `bottom_level`, `water_level`, and
`position`. The writable node property is `demand`. The original file is copied
to the path supplied with `-o`, then the requested value is saved there.

The same command works directly with EPANET input files:

```bash
./build/staci -m network.inp \
  -e PIPE74 -p diameter -n 0.8 -o network-modified.inp
```

Modification values use STACI SI units, so pipe diameter and tank levels are
specified in metres and demand in m³/h. STACI converts values back to the units
declared by the EPANET file. All other `.inp` sections—including controls,
patterns, coordinates, and metadata—are preserved. Supported EPANET write-back
properties are pipe `diameter`, junction `demand`, and tank
`bottom_level`/`water_level` (using the imported `EPANET_TANK_ID` element ID).

### Transport calculations

```bash
# Residence time
./build/staci -t network.spr

# Concentration distribution
./build/staci -c network.spr
```

For native SPR chemical calculations, `cl_input` and stored concentrations use
`kg/m3`, `cl_k` is a bulk first-order coefficient in `1/s`, and `cl_w` is a
first-order wall coefficient in `m/s`. `cl_length` is the chemical simulation
duration in hours. The legacy SPR command holds the solved hydraulic state
constant; use EPANET EPS mode when demands or operating states vary with time.

The `-t` command is the legacy fixed-hydraulic-state SPR transport calculation.
For an EPANET extended-period water-age simulation, set `QUALITY AGE` in the
INP `[OPTIONS]` section and run `--epanet-eps`; node and link ages are then
written to CSV and HDF5 in SI seconds.

For chlorine, set `QUALITY CHEMICAL Chlorine mg/L` and use `[QUALITY]`,
`[SOURCES]`, and `[REACTIONS]` normally. The same EPS command writes node and
link concentrations to CSV and HDF5 as SI `kg/m3`.

### Sensitivity calculations

```bash
# Sensitivity to one element property
./build/staci -r network.spr -e PIPE74 -p diameter

# Hydraulic solution, residence time, and demand sensitivity
./build/staci -d network.spr
```

### Export connectivity

```bash
./build/staci -x network.spr
```

This produces `nodelist.txt` and `connected_nodes.txt`. The scripts in
`python_tools/` can be used for additional connectivity analysis.

## Command-line reference

| Option | Purpose |
| --- | --- |
| `-s`, `--stac FILE` | Steady-state hydraulic simulation |
| `-i`, `--ini FILE` | Initialization/result file used with a hydraulic run |
| `-t`, `--travel_time FILE` | Residence-time calculation |
| `-c`, `--conc_transp FILE` | Concentration transport calculation |
| `-m`, `--mod_prop FILE` | Modify a property; also requires `-e`, `-p`, `-n`, and `-o` |
| `-l`, `--list_all_elements FILE` | List all nodes and edges |
| `-g`, `--get_data FILE` | Read a property; also requires `-e` and `-p` |
| `-r`, `--sensitivity FILE` | Parameter sensitivity; also requires `-e` and `-p` |
| `-d`, `--demand_sensitivity FILE` | Hydraulic, residence-time, and demand-sensitivity calculation |
| `-x`, `--export_for_connectivity_check FILE` | Export node connectivity |
| `-y`, `--export-epanet FILE` | Export a network to the `.inp` path supplied with `-o` |
| `-z`, `--epanet-eps FILE` | Run EPANET extended-period hydraulics and write chunked HDF5, metadata JSON, and SI CSV files using the `-o` prefix |
| `-q`, `--steady-quality FILE` | Solve asymptotic AGE and/or CHEMICAL quality for fixed hydraulics; use `--quality-mode`, optionally `--quality-sensitivity`, `-e`, and `-p` |
| `-e`, `--element_ID ID` | Select a node or edge |
| `-p`, `--property_ID NAME` | Select a property |
| `-n`, `--newValue VALUE` | New numeric property value |
| `-o`, `--outfile FILE` | Output network file for property modification |

## Generated files

STACI writes files relative to the current working directory or beside the
input filename, depending on the command:

| File | Meaning |
| --- | --- |
| `INPUT.ros` | Steady-state log/output |
| `INPUT.rps` | Steady-state progress |
| `INPUT.rrs` | Steady-state convergence marker |
| `INPUT.rot`, `INPUT.rpt`, `INPUT.rrt` | Residence-time output, progress, and completion marker |
| `INPUT.roc`, `INPUT.rpc`, `INPUT.rrc` | Concentration output, progress, and completion marker |
| `element_list.txt` | Element list created by `-l` |
| `tmp.dat` | Numeric value created by `-g` |
| `nodelist.txt`, `connected_nodes.txt` | Connectivity exports created by `-x` |
| `PREFIX-nodes.csv`, `PREFIX-links.csv` | EPS node and link time series |
| `PREFIX-tanks.csv`, `PREFIX-summary.csv` | EPS tank time series and run summary |
| `PREFIX.h5` | Chunked `STACI EPS OUTPUT v1` data for visualization |
| `PREFIX.meta.json` | EPS run metadata, dimensions, status codes, warnings, and value ranges |
| `PREFIX-steady-nodes.csv`, `PREFIX-steady-links.csv` | Steady-quality SI node and link results |
| `PREFIX-steady-summary.csv` | Steady-quality mode and dimensions |
| `PREFIX-steady-sensitivity.csv` | Optional water-quality and hydraulic-flow derivatives per SI parameter |
