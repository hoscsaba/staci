# STACI native C++ examples

These programs demonstrate the reusable C++ API without MATLAB. Build them with:

```bash
cmake -S . -B build -DSTACI_BUILD_CPP_EXAMPLES=ON
cmake --build build --target \
  staci_example_hydraulics \
  staci_example_epanet_inp --parallel
```

## Basic SPR hydraulics

`basic_hydraulics.cpp` loads a native SPR network through `StaciSession`, solves
the steady hydraulic state and prints SI node/link results:

```bash
./build/staci_example_hydraulics
./build/staci_example_hydraulics path/to/network.spr
```

## EPANET INP import

`epanet_inp_hydraulics.cpp` separately demonstrates EPANET INP parsing, the
resulting STACI element types, a steady hydraulic snapshot and SI results:

```bash
./build/staci_example_epanet_inp
./build/staci_example_epanet_inp path/to/network.inp
```

The default `epanet_eps_smoke.inp` intentionally contains a demand pattern,
reservoir head pattern, tank, power pump and simple control, so import warnings
also demonstrate which features belong to EPS rather than one steady snapshot.

## Official EPANET versus STACI EPS

`epanet_staci_eps_compare.cpp` runs STACI extended-period simulation and the
official OWA EPANET 2.2 toolkit on the same INP file. It matches time/element
identifiers and reports maximum absolute SI differences for total head,
junction pressure, junction demand, flow and velocity. The program returns 0
when the included benchmark tolerances pass and 2 when compared values exceed
them.

The target is enabled only when CMake finds the official EPANET headers and
library. The repository's pinned reference build can be prepared with:

```bash
python3 tests/setup_epanet_reference.py
cmake -S . -B build \
  -DSTACI_BUILD_CPP_EXAMPLES=ON \
  -DSTACI_EPANET_TOOLKIT_ROOT="$PWD/build/epanet-reference"
cmake --build build --target staci_example_epanet_eps_compare --parallel
```

Run the default three-state EPS benchmark:

```bash
./build/staci_example_epanet_eps_compare \
  tests/epanet_eps_smoke.inp \
  cpp-eps-results
```

The second argument is the output directory for STACI's chunked HDF5 (when
available), metadata JSON and SI CSV result files and EPANET's diagnostic
`epanet.rpt`. The comparison itself calls the EPANET toolkit directly and does
not parse the formatted report.

All three examples are also registered as CTest cases when their dependencies
are available.
