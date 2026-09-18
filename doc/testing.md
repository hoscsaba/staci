# STACI testing reference

[Back to the README](../README.md)

Commands below assume the repository root as the working directory unless stated otherwise.

## Run the integration tests

The CTest suite includes deterministic, fixed-seed pagmo tests for both
optimizer programs. The calibration test uses the small `anytown_1med.spr`
network, while the splitter test uses `LOV-LOVOTV-2-input_mod.spr`, which has
2486 edge elements (2479 elements participate in splitting after excluding
fixed-pressure and pool boundary elements):

```bash
ctest --test-dir build --output-on-failure
```

The test configurations use very small populations and generation counts; they
verify integration and file output, not optimization quality.

The focused steady-quality checks can be run separately:

```bash
ctest --test-dir build -R epanet_steady_quality --output-on-failure
```

They verify the sparse equations against closed-form pipe/mixing results,
compare the asymptotic result with the existing 24-hour time-marching AGE and
CHEMICAL models, and compare analytic diameter sensitivities with centered
finite differences. The readable benchmark log is
`build/steady-quality-test/steady-quality-test.log`.

The portable Python test runner recursively tests every `.spr` and `.inp` file
under `tests/`. It uses copies below `tests/test-results/`, so hydraulic
calculations never overwrite the original test networks. The runner also tests
`staci_calibrate` and runs `staci_split` on the large LOV network with both two
and three requested segments. On Linux and macOS:

```bash
python3 tests/run_tests.py --binary ./build/staci
```

On Windows:

```powershell
py tests\run_tests.py --binary build\Release\staci.exe
```

SPR files up to 1 MB receive a steady-state hydraulic test, convergence-marker
check, EPANET export, and round-trip import. Larger SPR files receive the export
and round-trip test by default, avoiding multi-minute hydraulic runs. To force
hydraulics for every SPR file, add `--full-spr-hydraulics` and optionally raise
the per-command limit with `--timeout SECONDS`.

Every INP file receives a byte-for-byte metadata round-trip, import, and
extended-period simulation test. The runner verifies the SI CSV schema, JSON
metadata, and convergence counts. A
separate HDF5 test generates an EPS result and verifies its file signature,
required datasets, SI unit attributes, extendible time dimension, and
16-frame chunking. This check requires either the Python `h5py` package or the
`h5dump` command-line utility; use the `H5DUMP` environment variable if the
utility is not on `PATH`. The HDF5 test output is always retained in a
directory below `tests/test-results/hdf5/`.

The runner continues after individual failures and writes the full
human-readable report to `tests/test-results/run_tests.log`. Before every run,
the complete `tests/test-results/` directory is deleted and recreated, so it
contains only the latest results. Generated artifacts are grouped as follows:

- `networks/`: SPR and INP integration-test copies and their outputs;
- `hdf5/`: the retained chunked EPS HDF5 test file and related SI outputs;
- `epanet-reference/`: official EPANET reports, STACI SI results, and a
  machine- and human-readable `comparison.json` for each validation network;
- `calibration/`: calibration settings, logs, and best result;
- `split/2-segments/` and `split/3-segments/`: splitter settings, membership,
  logs, and `network-2-segments.png` or `network-3-segments.png`.

The two PNG files are generated directly from the SPR coordinates, links, and
the splitter's `membership.txt`; no plotting package is required. The splitter
uses one graph seed per segment and assigns nodes with a multi-source graph
traversal, so every segment is connected by construction. The test independently
checks this property and fails if any requested segment contains multiple graph
components. Old generated STACI sidecar logs (`.ros`, `.rps`, `.rrs`, and
related variants) are also removed at startup.

### Validate STACI against official EPANET 2.2

The test suite can build the official
[OpenWaterAnalytics EPANET 2.2](https://github.com/OpenWaterAnalytics/EPANET)
solver locally and compare its hydraulic results with STACI. The helper pins
the source to commit `4d8d82ddc260fce216af9321fc3d9a4646ac6827`, so the
reference calculation is reproducible. It uses CMake and the platform's C/C++
compiler; on Apple Silicon this creates a native `arm64` executable and needs
neither Rosetta nor a system-wide EPANET installation.

```bash
cmake --build build --target setup_epanet_reference
ctest --test-dir build -L epanet --output-on-failure
```

The same setup and tests work with a Visual Studio CMake build on Windows and
with GCC or Clang on Linux. A separately installed solver can be selected at
configure time:

```bash
cmake -S . -B build -DSTACI_EPANET_EXECUTABLE=/path/to/runepanet
```

The portable Python suite discovers this local solver automatically. Its
location can also be supplied explicitly; `--require-epanet-reference` makes a
missing reference solver an error instead of a documented skip:

```bash
python3 tests/setup_epanet_reference.py
python3 tests/run_tests.py --binary ./build/staci \
  --epanet-binary ./build/epanet-reference/build/bin/runepanet \
  --require-epanet-reference
```

Thirteen targeted fixtures validate every EPANET hydraulic element that currently
has a calculation-model counterpart in STACI, together with the supported
simple-control and rule paths:

| Fixture | Direct EPANET comparison coverage |
|---|---|
| `epanet_reference_pipe.inp` | Junction, reservoir, and two Hazen-Williams pipes |
| `epanet_reference_tcv.inp` | Pipe minor loss and an active TCV |
| `epanet_reference_pattern.inp` | Multi-state junction demand pattern |
| `epanet_reference_anytown.inp` | 21 junctions, two reservoirs, 41 Darcy-Weisbach pipes, one HEAD pump, seven EPS states, and direct `QUALITY AGE` comparison |
| `epanet_reference_anytown_chlorine.inp` | The same AnyTown EPS hydraulics with a 1 mg/L reservoir boundary, first-order bulk chlorine decay, and direct SI chemical-quality comparison |
| `epanet_node_metadata.inp` | Separate demand categories, demand multiplier, initially closed pipe, minor loss, and CV pipe |
| `epanet_pumps.inp` | POWER pump, multi-point HEAD pump, one-point HEAD pump, initial speed, and speed patterns |
| `epanet_eps_smoke.inp` | Tank storage and level, reservoir head pattern, demand pattern, POWER pump, and tank-level simple control |
| `epanet_controls.inp` | Elapsed-time, clock-time, junction-pressure, OPEN/CLOSED, and numeric pump-setting controls |
| `epanet_rules.inp` | Time, clock, demand, pressure, flow, link status, tank level/fill/drain time, AND/OR/ELSE, equality, priority conflicts, pump setting, WNTR-style `=` actions, and `HH:MM:SS` times |
| `epanet_volume_curve.inp` | Non-cylindrical tank storage, CMH input conversion, piecewise-linear EPANET volume-depth interpolation, SI volume output, and direct EPANET reference comparison |
| `epanet_reference_tcv_controls.inp` | Numeric and OPEN TCV simple controls |
| `epanet_reference_tcv_rules.inp` | TCV setting/status premises and actions plus multiple-action THEN/ELSE branches |

At every report time the tests compare SI node head, junction pressure and
demand, link flow, velocity magnitude, signed total head loss, and convergence.
For control and rule fixtures, STACI's enabled/closed link state is also checked
exactly against EPANET's `STATUS FULL` event history. The comparison JSON stores
maximum errors both globally and separately for every node and link, so a
regression is localized to its element and result quantity.

The Anytown reference also runs the official EPANET 2.2 water-quality solver in
`AGE` mode. STACI advances plug-flow pipe parcels with perfect instantaneous
junction mixing and zero-age reservoir boundaries using the EPS hydraulic state
from each interval. Node ages are compared in SI seconds with a 600 s absolute
acceptance limit; the current maximum difference is about 525 s over the
six-hour simulation. This tolerance includes the effect of the independently
computed hydraulic flows, whose largest Anytown difference is about 5%.

The separate AnyTown chlorine benchmark starts with zero concentration in the
network, applies a 1 mg/L `CONCEN` source at reservoir `NODE210`, and uses a global first-order
bulk coefficient of -0.50/day. EPANET report values in mg/L are converted before
comparison with STACI's SI `kg/m3` output. With a shared 60 s quality timestep,
the current maximum node-concentration difference is `6.15e-5 kg/m3` over six
hours; the automated acceptance limit is `8e-5 kg/m3`.

The three original small fixtures retain strict absolute tolerances: 0.002 m
for head and head loss, `1e-9 m3/s` for flow, and `5e-6 m/s` for velocity.
Element, pump, tank, control, and rule fixtures have explicit case-level limits
based on their current observed solver differences; exact demand, convergence,
and enabled/closed state checks remain independent of those numeric limits.
Anytown uses large-network acceptance limits of 1% for head, 0.016 m3/s or 3%
for flow, 0.08 m/s or 3% for velocity, and 0.4 m or 5% for link head loss. Its
current observed maxima are 0.566 m, 0.0153 m3/s, 0.073 m/s, and 0.377 m,
respectively. These separate limits document the larger difference between
STACI's and EPANET's independently implemented Darcy-Weisbach and pump models
while still detecting hydraulic, EPS-pattern, control/rule, status, or
unit-conversion regressions. CTest marks these cases skipped if `runepanet` has
not been built; the rest of the suite remains usable offline.

The command `ctest --test-dir build -L epanet --output-on-failure` runs only
these targeted compatibility checks. It does not run the complete STACI
regression suite.

### Channel-only hydraulic reference tests

The circular open-channel element has a separate reference suite derived from
`tests/channel.spr`. It compares STACI against independently implemented
circular geometry, Manning normal depth, critical depth, gradually varied flow,
and Darcy-Weisbach pressure-flow calculations. It covers both flow directions,
open/full transitions, fully pressurized cases, hydrostatic equilibrium, and a
momentum-matched hydraulic jump where a smooth GVF profile is impossible:

```bash
python3 tests/run_channel_tests.py --binary ./build/staci
```

On Windows with a Visual Studio Release build:

```powershell
py tests\run_channel_tests.py --binary build\Release\staci.exe
```

Results are retained in `tests/test-results/channel/`. All channel cases are
mandatory and any reference discrepancy makes the command fail. See
[tests/CHANNEL_TESTS.md](../tests/CHANNEL_TESTS.md) for the case matrix, governing
equations, and engineering references. The same suite is registered in CTest
as `channel_reference_suite`. Each case also produces two longitudinal-section
files, `channel-profile.svg` and `channel-profile.pdf`. Each visualization
starts with a labeled topology sketch and calculated flow arrows. Channel
labels include the calculated SI discharge, while node labels show the
flow-oriented channel-end invert elevations and water depths (`z_e`, `h_e`,
`z_v`, `h_v`). The longitudinal diagrams use absolute SI elevations, show the
channel invert, crown, calculated water surface and red dashed energy grade
line, and always orient the horizontal axis in the calculated flow direction.
Red endpoint markers identify the upstream rest level
`z_e + h_e + v_e^2/(2g)` and downstream rest level `z_v + h_v`.

The standard suite also contains
`tests/channel_network_merge_split.spr`, a stationary network of six
open-surface channels. Two channels merge at the central manhole and two
channels leave it, after which both branches continue through another channel.
The two source channels deliberately have different bed slopes: `CHANNEL_1`
is a mild 0.5% reach and `CHANNEL_2` is a steep 1.5% reach. Their source water
levels preserve the same 0.5 m inlet depth, isolating the effect of slope.
One parallel branch is adverse: the `CHANNEL_5` invert rises from 1.00 m to
1.05 m along the actual `BRANCH_A` to `SINK_A` flow direction, corresponding
to a -0.05% bed slope. `CHANNEL_6` retains a conventional falling bed and the
branch discharges become asymmetric.
The dedicated test verifies the topology, convergence, nonzero channel flows,
open water depths at every channel end, and mass conservation at the 2-in/2-out
junction. It also verifies the actual flow orientation and negative
flow-direction bed slope in `CHANNEL_5`, and checks that its calculated profile
remains finite, open and positive-depth. The inlet slope values and their
minimum ratio are checked explicitly as well:

```bash
python3 tests/run_channel_network_test.py --binary ./build/staci
```

It is registered in CTest as `channel_network_merge_split` and is also executed
by `tests/run_tests.py` as the `CHANNEL-NETWORK` case. Its retained
`channel-network-longitudinal-profile.svg` and the multi-page
`channel-network-longitudinal-profile.pdf` begin with a topology overview that
labels every channel and junction, includes every calculated channel discharge,
and lists the applicable `z_e`, `h_e`, `z_v`, `h_v` endpoint data at each node.
Nodes shared by incoming and outgoing channels show both endpoint roles. The
following elevation-correct panels contain the red dashed energy grade line and
rest-level markers for every source-to-outlet route, so the channels are shown
consecutively even through the merge and split. The standalone renderer can be
used for any solved SPR network containing `channel1` elements:

```bash
python3 tests/plot_channel_profiles.py \
  --network tests/test-results/channel-network/channel_network_merge_split.spr \
  --output tests/test-results/channel-network/channel-profile.svg \
  --pdf-output tests/test-results/channel-network/channel-profile.pdf
```
