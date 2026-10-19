# Start a calibration with your own network

`staci_calibrate` fits **pipe diameters**, using separate steady network snapshots
and measured pressure/level time series. It does not fit pipe roughness or run an
INP extended-period simulation automatically. Use at least two different
operating conditions: the pressure objective removes each sensor's mean, so one
pressure sample provides no fitting information.

## 1. Run a supplied example first

Install/build the programs with optimizer support (see the
[installation guide](../README.md#installation)). Verify `staci_calibrate --help`.
The examples are in the source repository; the binary archive alone contains
executables, not these example inputs.

From the repository root on macOS/Linux:

```sh
cp -R examples/calibration/anytown /tmp/my-calibration
cd /tmp/my-calibration
/absolute/path/to/staci_calibrate --settings settings.json --seed 12345
```

In Windows PowerShell, choose a new, nonexistent job folder:

```powershell
Copy-Item -Recurse examples\calibration\anytown C:\jobs\my-calibration
Set-Location C:\jobs\my-calibration
& 'C:\staci\bin\staci_calibrate.exe' --settings settings.json --seed 12345
$LASTEXITCODE
```

Keep the executable's dependency DLLs beside it on Windows. The prepared
[Anytown example](../examples/calibration/anytown/README.md) needs no Python;
Python is only required to regenerate its synthetic measurements. It has three
periods, two selected pipes, six pressure sensors and two tanks. A smaller
[two-pipe example](../examples/calibration/README.md) is also provided.
The examples demonstrate formats and recovery from synthetic data, not validation
against field measurements.

## 2. Prepare a separate directory for your real job

Use a new directory so that fixed-name outputs cannot overwrite another run:

```text
my-real-job/
  settings.json
  measurements.json
  network_0.inp
  network_1.inp
```

Export one snapshot for each measured operating condition. All snapshots must
have the same node/link IDs, order, types and connections. Set demands, pump
speeds, source heads and valve settings to the actual conditions for that period.
Do not expect patterns, controls or INP duration to generate the snapshots for
you. Use matching measurement times and a uniform interval `dt` in hours.

First run each snapshot hydraulically from this directory:

```sh
/absolute/path/to/staci -s network_0.inp
/absolute/path/to/staci -s network_1.inp
```

Resolve initial hydraulic failures before fitting. Review warnings and check
that the modeled operating conditions are plausible. Calibration cannot
compensate reliably for incorrect topology, source heads or demand assumptions.

## 3. Supply measured data in SI units

Create `measurements.json`. This is an illustrative two-sensor, two-period file;
**replace IDs and values with your own measurements**:

```json
{
  "measurements": [
    {"id": "J1", "type": "node", "values": [30.2, 28.9]},
    {"id": "J2", "type": "node", "values": [25.4, 23.8]}
  ]
}
```

- `node` values are **pressure heads in metres**, not total hydraulic head or bar.
  For water at density 1000 kg/m³, 1 bar is approximately 10.197 m pressure head.
  Convert gauge pressure consistently; do not add ground elevation.
- `pool` values are **water depths above tank bottom in metres**, not elevations.
  For example, a tank row can be `{"id":"T1","type":"pool","values":[3.0,2.9]}`.
- Use original INP IDs exactly, including case. Preserve string IDs such as `"01"`.
- Each row needs a finite number for every selected period; missing/blank values
  and NaN are rejected. The first array entry corresponds to `network_0.inp`
  when `Start_of_Periods` is zero.
- All measurement values and diameter bounds use SI even if the INP uses US units.

## 4. Choose pipes and configure the fit

This is a complete pressure-only, two-period settings file. Replace `P1`, `P2`,
the bounds and the time step with values appropriate for your model:

```json
{
  "global_debug_level": 1,
  "Staci_debug_level": 0,
  "dir_name": "./",
  "fname_prefix": "network_",
  "network_extension": ".inp",
  "logfilename": "calibrate.log",
  "best_logfilename": "calibrate-best.log",
  "sollwert_dfile": "measurements.json",
  "Start_of_Periods": 0,
  "Num_of_Periods": 2,
  "Spoil_Active_Pipes": false,
  "dt": 1.0,
  "weight_p_err": 1.0,
  "type_of_pipe_selection": "list",
  "pipe_ids": ["P1", "P2"],
  "min_diameter_m": 0.08,
  "max_diameter_m": 0.60,
  "popsize": 20,
  "ngen": 25,
  "pmut": 0.2,
  "pcross": 0.8
}
```

Start with a small number of uncertain pipes and enough independent measurements.
To fit every pipe, use `"type_of_pipe_selection": "all"` and remove `pipe_ids`;
a larger problem needs more data and optimization effort.

`min_diameter_m` and `max_diameter_m` are absolute bounds for all selected pipes.
Delete either key to use that pipe's original diameter times 0.5 or 2 for the
missing bound. The effective lower bound must be positive and smaller than the
upper bound. `Dmin` is a separate pipe-selection option, not a fitted lower bound.

Use `weight_p_err: 1.0` for pressure-only data, `0.0` for tank-only data, or a
value between them for both. Measured tank levels initialize each candidate at
the first period; subsequent levels use a constant-area tank balance and `dt`.
The pressure objective fits changes over time, ignoring constant sensor offsets.
Read the [objective and tank limitations](calibration.md#objective-and-period-handling)
before interpreting the fitted diameters physically.

The population/generation values above are starting values, not a convergence
guarantee. Keep `Spoil_Active_Pipes` false for real jobs. Increase `Num_of_Periods`
and add matching snapshots/measurement entries for additional periods.

## 5. Run, inspect and repeat

Change into your job directory (passing an absolute settings path does not change
the working directory), then run:

```sh
/absolute/path/to/staci_calibrate --settings settings.json --seed 12345 --diagnostics-file diagnostics.jsonl
```

On Windows use the PowerShell `& 'C:\...\staci_calibrate.exe'` invocation above
with the same arguments. A zero exit code means execution completed, not that
the parameters are uniquely identified or the global optimum was found.

| File/result | What to check |
| --- | --- |
| `diagnostics.jsonl` | Warnings and readable errors; exit 2 means input/configuration error, exit 1 means calculation failure |
| `calibrate.log` | Selected pipe IDs, effective diameter intervals, initial solves and progress |
| `best.dat` | Pipe IDs in network order; fitted diameters in metres; ratios to original diameters; measured/computed comparisons |
| `calibrate-best.log`, `bog.dat` | Best-candidate history and available optimization progress |
| `sensitivity_Pressures.dat`, `sensitivity_MassFlowRates.dat` | Final sensitivity information |

If a diameter sits on a bound, check model/data consistency and identifiability
before widening the interval. Repeat with other seeds and increased population
or generation counts in **separate job directories**. Compare residuals and
parameter stability, and validate against operating conditions not used for fitting.

The program does **not** write a calibrated INP automatically. Apply reviewed
fitted diameters to copies of every relevant snapshot, converting metres back
to the INP diameter convention (millimetres for metric INP, inches for US INP).
Then solve and validate those copies. Preserve original models and measurements.

For all options and the legacy XML/CSV equivalents see the
[calibration reference](calibration.md); for Python/MATLAB/backend invocation see
[README_DEV](../README_DEV.md).
