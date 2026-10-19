# INP calibration example

Copy this directory to a writable job directory, change into that copy, then run:

```sh
/absolute/path/to/staci_calibrate --settings settings.json --seed 12345
```

The two INP files have the same topology and different demands. The JSON
measurements are synthetic pressure heads and tank levels in metres, computed
from the H-W pipe equations for the supplied diameters. The period length is
0.01 hours (36 seconds). The largest pipe (`Large`) is calibrated.

This small example demonstrates execution and result formats; it does not
establish identifiability or guarantee an optimum. See the
[calibration guide](../../doc/calibration.md) for the objective, required data,
INP snapshot limitations and output files. Results are written in the job
working directory; the input INP files are not replaced by fitted networks.

For a larger looped network with pumps, two tanks and three operating periods,
see the [Anytown example](anytown/README.md), which selects pipes explicitly.
