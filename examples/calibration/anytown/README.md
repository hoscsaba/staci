# Anytown: three-period diameter calibration

This larger example supplements the [two-pipe example](../README.md).
It has 22 junctions, one reservoir, two tanks, 43 pipes and three pumps.

## Run

Copy this directory into a writable job directory, change into the copy, then:

```sh
/absolute/path/to/staci_calibrate --settings settings.json --seed 12345
```

The prepared inputs are self-contained; Python is only needed to regenerate
synthetic measurements. `settings.json` uses `type_of_pipe_selection: "list"`
and `pipe_ids: ["1", "2"]`. The three snapshots use demand factors 0.75, 1.0
and 1.15, separated by 0.01 h (36 s). Measurements contain six pressure-head
series and both tank-level series, all in metres. The INP files retain Anytown's
GPM/feet/inches convention; the importer converts these to SI internally.

The known diameters of pipes 1 and 2 are both 0.6096 m. Their input diameters
are deliberately changed to 0.67056 m and 0.54864 m (+10% and -10%). The example
uses population 20 and 25 generations. A macOS run with seed 12345 returned
0.609500 m and 0.609443 m, respectively. Numerical/optimizer results can vary
across platforms; these are illustrative, not exact golden values.

Read `best.dat` for fitted diameters, ratios to the perturbed input values, and
measurement comparisons. Configured logs and sensitivity files are written in
the job directory. The program does not automatically export a fitted INP.
For all 43 pipes, set `type_of_pipe_selection` to `all` and remove `pipe_ids`;
this substantially enlarges the fitting problem and does not guarantee unique
identification from these measurements.

## Source and synthetic-data construction

`anytown.inp` is an exact copy of the repository's physically adjusted
`tests/public_networks/wntr/wntr/tests/networks_for_testing/Anytown.inp`, not the
unmodified upstream model. The original corpus file is untouched. Its source
hash, model counts and perturbations are recorded in `provenance.json`.
See [the physical model history](../../../doc/anytown_physical_model.md).

`prepare.py` derives three steady operating snapshots. It removes time patterns,
controls/rules and quality settings, fixes all pumps at nominal speed, and sets
initial tank depths to 20 ft (6.096 m). Demands are scaled explicitly in each
snapshot. The initial levels are inside the original 10–35 ft bounds.

**Measurements are synthetic STACI results, not field measurements or an
independent EPANET validation.** The generator solves the unperturbed-diameter
model for each period, updating tank levels from the preceding tank flow and
constant area, exactly as the calibration period model does. It then writes
calibration inputs with the two perturbed diameters. It does not fit an EPS
simulation or claim contaminant/water-quality agreement.

Regenerate into a new/empty directory (the source file is read-only):

```sh
python3 prepare.py --staci /absolute/path/to/staci --output /absolute/path/to/new-job
```

The standard CTest case `calibration_anytown_example` regenerates the data in
isolation, verifies topology and source preservation, runs the prepared example,
and checks recovered diameters within 1% and pressure/level agreement within
0.1 m. See the [calibration contract](../../../doc/calibration.md).
