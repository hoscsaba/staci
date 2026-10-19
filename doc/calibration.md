# Calibration with SPR or EPANET INP networks

`staci_calibrate` varies selected **pipe diameters** using pagmo's genetic
algorithm. Every candidate uses the same diameters across all supplied periods.
Bounds are 0.5 to 2 times each original diameter. It does not fit roughness,
demands, pump curves or chlorine parameters.

## Inputs and invocation

```sh
staci_calibrate --settings settings.json --seed 12345
```

Settings may be JSON or legacy XML. See the runnable
[two-pipe INP example](../examples/calibration/README.md), the larger
[Anytown example](../examples/calibration/anytown/README.md), and the
[general settings template](../examples/config/staci_calibrate_settings.json).
The JSON keys are also XML child names under `<settings>`.

For a numbered INP snapshot series, use:

```json
{"dir_name": "./", "fname_prefix": "network_", "network_extension": ".inp",
 "Start_of_Periods": 0, "Num_of_Periods": 2, "dt": 1.0}
```

These are the network-related fields of a complete settings file. They load
`network_0.inp` and `network_1.inp`. The default extension remains `.spr`.
Supported values for `network_extension` are `.spr` and `.inp`.
`Start_of_Periods` also selects the first measurement column (zero-based).
All snapshots must retain identical node/link IDs, order, types and endpoints;
period-specific operating data such as demands may differ.

For one explicitly named snapshot use `"network_file": "network.inp"` and
`"Num_of_Periods": 1`. `fname_prefix` is then unused and may be omitted;
`network_extension` must be omitted. The extension selects the importer.
The file is relative to `dir_name`, or may be absolute. Settings paths do not
change the process working directory. Measurement paths retain the historical
`dir_name` + `sollwert_dfile` convention, so include a trailing separator in
`dir_name` when supplying a directory.

Required measurements are specified by `sollwert_dfile` (JSON or semicolon CSV):

```json
{"measurements": [
  {"id": "J1", "type": "node", "values": [30.2, 29.8]},
  {"id": "T1", "type": "pool", "values": [3.0, 3.1]}
]}
```

These illustrative values are **pressure head in metres**, not bar, and tank
water depth in metres above the bottom. Use the original INP junction/tank IDs;
the internal `EPANET_TANK_` prefix is not required. Either measurement type may
be used; choose weights accordingly. Each row needs at least
`Start_of_Periods + Num_of_Periods` finite values. CSV equivalents are
`J1;node;30.2;29.8;` and `T1;pool;3.0;3.1;`.

## Explicit pipe selection

To calibrate specified pipes, use these fields in the settings JSON:

```json
"type_of_pipe_selection": "list",
"pipe_ids": ["Small", "Large"]
```

To calibrate every pipe:

```json
"type_of_pipe_selection": "all"
```

These modes do not require `num_of_active_pipes` or `Dmin`; any such existing
fields are ignored. `pipe_ids` is only valid in `list` mode. IDs are exact and
case-sensitive. Empty lists, duplicate IDs, unknown IDs and non-pipe elements
are input errors (exit 2, `INPUT.CONFIG`). Results use the network's pipe order,
regardless of the order in the supplied list. Pumps, valves and tank boundaries
are never included by `all`. Both SPR and INP networks are supported.

The equivalent legacy XML uses a JSON array as the element text:

```xml
<type_of_pipe_selection>list</type_of_pipe_selection>
<pipe_ids>["Small", "Large"]</pipe_ids>
```

For all pipes use `<type_of_pipe_selection>all</type_of_pipe_selection>` and
omit `pipe_ids`. Other settings remain the same.

Other required settings are:

- `type_of_pipe_selection`: `list`/`all` as above, or `largest_diameter` or `most_sensitive` with
  `num_of_active_pipes`, or `Dmin` with a diameter threshold in metres.
- `weight_p_err`: pressure-error weight, between 0 and 1; tank weight is its complement.
- `popsize`, `ngen`, `pmut`, `pcross`: population, generation count, mutation and crossover probabilities.
- `Spoil_Active_Pipes`: normally false; enables the legacy diameter-perturbation experiment.
- `logfilename`, `best_logfilename`, `global_debug_level`, `Staci_debug_level`.

## Objective and period handling

The pressure term is the sum of absolute differences between **mean-centered
pressure time series**, separately for each measured node. The tank term is
the sum of absolute water-level errors. Their weighted sum is minimized.
Constant pressure offsets are therefore ignored. **One period always gives
zero pressure error**, and each measured tank's first level initializes the
model, so a one-period run is a loading/smoke check rather than an informative
fit. `CALIBRATION.UNINFORMATIVE_PRESSURE` warns about the one-period pressure
case. Use at least two different operating conditions and sufficient independent
measurements; multiple periods alone do not guarantee identifiable parameters.

Each candidate resets measured initial tank levels. Later levels use the
preceding solved tank flow, constant surface area and `dt` (hours) before the
next hydraulic solve. This is the existing calibration period model, **not an
EPANET EPS run**. INP duration does not automatically create periods or sample
patterns/controls; supply separate snapshots with the intended operating data.
The common steady INP importer applies its usual compatibility rules and
warnings. Variable-area tank storage and full EPS controls are not implemented
by this calibration loop.

## Results and errors

`best.dat` contains active pipe IDs, fitted diameters in metres, diameter ratios,
and measured/computed comparisons. `bog.dat` contains available optimizer
progress records; configured logs track settings and best candidates.
`sensitivity_MassFlowRates.dat` and `sensitivity_Pressures.dat` contain final
sensitivities. Fixed-name outputs use the process working directory: use separate
job directories for parallel runs. No fitted network is automatically exported.

Common JSONL diagnostics apply (`--diagnostics-file PATH`). Missing networks,
invalid format settings or incompatible period topology return exit code 2;
initial hydraulic nonconvergence returns 1. Review warnings and the fitted
parameters even after a successful exit. A genetic algorithm result is not a
guarantee of a global optimum or a uniquely identified physical model.
