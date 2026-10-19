# STACI examples

This directory contains user-facing examples. Automated regression and
developer checks belong under `tests/`.

- [`cpp/`](cpp/README.md) demonstrates the native C++ API, EPANET INP import,
  and a direct official-EPANET-versus-STACI EPS comparison;
- [`matlab/`](matlab/README.md) demonstrates the in-memory MATLAB interface,
  hydraulic calculation, vectorized model changes, steady water quality,
  hydraulic sensitivities and diameter optimization;
- [calibration/](calibration/README.md) provides two-pipe and Anytown INP
  calibration examples, including explicit pipe selection.

API values and calibration measurements use SI units. INP network files retain
the units declared in their `[OPTIONS]` section (Anytown uses US units).
