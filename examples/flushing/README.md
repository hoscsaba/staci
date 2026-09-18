# Worked example: three hydrants on a branched network

This small, synthetic network demonstrates pressure-dependent hydrant discharge,
reverse flow, overlapping coverage, greedy ranking and both opening-time estimates.
It requires no project-specific files or external Python libraries.

## Run

Build STACI using the repository [build instructions](../../README.md), then run
from the repository root:

```sh
./build/staci_flush --inp examples/flushing/network.inp --config examples/flushing/flushing_config.json
```

On Windows with the Visual Studio generator:

```powershell
.\build\Release\staci_flush.exe --inp examples/flushing/network.inp --config examples/flushing/flushing_config.json
```

Results appear in `examples/flushing/results/`. This path is relative to the JSON
config, even when you launch the executable from another directory. Before another
run, move the previous results aside or choose a different `output_dir` in the
config; the executable refuses to overwrite a nonempty output directory.

## Network and settings

```text
R (40 m head) --- MAIN --- A --- BRANCH_B --- B
                         |
                         +---- BRANCH_C --- C
```

R is a fixed-head reservoir. Junctions A, B and C are at zero elevation, have no
customer demand and each have a hydrant. MAIN is 100 m long with 150 mm diameter.
BRANCH_B is 100 m long and BRANCH_C is 150 m long, both with 100 mm diameter.
All pipes use Hazen–Williams C = 120. BRANCH_C is stored from C to A in the INP,
so its flushing flow is negative when hydrant C is open. Absolute velocity is
used for coverage.

The config uses outlet area 0.002 m², total K = 2, velocity threshold 0.5 m/s,
zero minimum pressure head, and enables tagged EPANET network exports. K includes
outlet kinetic head. These are demonstration parameters.

## Expected results

The closed-hydrant baseline has zero flow. Opening a hydrant causes flow only
along its supply path. Both branches and MAIN exceed the threshold when active.
The network contains **3.73 m³** of pipe volume.

| Rank | Hydrant | Total qualifying volume (m³) | Additional volume (m³) | Cumulative volume (m³) | Coverage (%) | Travel time (min) | Volume/flow time (min) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | C | 2.95 | 2.95 | 2.95 | 78.95 | 1.9 | 1.9 |
| 2 | B | 2.55 | 0.79 | 3.73 | 100.00 | 1.5 | 1.5 |
| 3 | A | 1.77 | 0.00 | 3.73 | 100.00 | 0.8 | 0.8 |

C ranks first because MAIN + BRANCH_C has the largest volume. B comes next:
MAIN is already counted, so only BRANCH_B adds coverage. A remains in the full
ranking with `redundant = 1` because its only pipe, MAIN, is already covered.
Opening times still describe A's complete scenario, even though its additional
coverage is zero. Calculations retain full precision before formatting, so
summing the displayed rounded volumes can differ by 0.01 m³.

The two time estimates coincide here: along each serial supply path, every pipe
carries the same discharge, so `sum(length / velocity) = sum(pipe volume) / Q`.
They generally differ in networks with flow splitting or customer demands.
Approximate discharges are A: 37.774 L/s, B: 29.100 L/s, C: 26.435 L/s.

The machine-readable reference is [expected_plan.csv](expected_plan.csv).
Inspect `results/flushing_plan.txt` or `results/flushing_plan.csv` after running.
`results/pipe_travel_times.csv` contains the underlying transport calculations.

## Inspect in EPANET

Open `results/networks/network_hydrant_C.inp`. Only C has an emitter. MAIN and
BRANCH_C have the `FLUSHED` pipe tag; BRANCH_B has `BELOW_THRESHOLD`. Geometry is
included for graphical inspection. To highlight EPANET's recalculated coverage,
run its analysis, then use **View → Query → Links → Velocity → Above → 0.5**.
The saved tags represent STACI results and remain unchanged by EPANET analysis.

## Automated regression

```sh
ctest --test-dir build -R flushing --output-on-failure
```

Add `-C Release` for Visual Studio builds. The `flushing_worked_example` test runs
this config in a temporary directory, checks the expected plan, independently
solves the serial-pipe/outlet head balance, verifies reverse flow and emitter/tag
exports, and checks that original example inputs are unchanged. Other flushing
tests cover overlap reranking, slow connecting paths, cycles, unreachable pipes,
pressure filtering and malformed inputs. GitHub runs them on all three platforms.
