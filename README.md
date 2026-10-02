# STACI

STACI is a command-line simulator for steady-state hydraulic and transport
calculations in pressurized and open-channel networks. Its native network and
result format is XML-based `.spr`; EPANET `.inp` network files can also be
imported, modified, exported, and used for extended-period simulations.

The project is written in C++17 and built with CMake. It supports GCC and Clang
on Linux and macOS, and MSVC on Windows. The numerical solver uses UMFPACK from
SuiteSparse.

Compared with EPANET's water-distribution focus, STACI additionally provides:

- native stationary open-channel elements with rectangular, circular and
  arbitrary cross-sections, plus overflow/weir elements;
- hydraulic parameter and demand-sensitivity calculations;
- pagmo2-based pipe calibration and graph-aware network partitioning for
  sensor-placement studies;
- chunked, incrementally written HDF5 EPS results with SI arrays and explicit
  node/link index mappings.

## Table of contents

- [Standalone programs](#standalone-programs)
- [Installation](#installation)
  - [Quick start and build options](#quick-start-and-build-options)
  - [Precompiled executables](#precompiled-executables)
  - [Requirements](#requirements)
  - [Linux](#build-on-linux)
  - [macOS](#build-on-macos)
  - [Windows](#build-on-windows)
  - [Install to a staging directory](#install-to-a-staging-directory)
  - [Non-standard SuiteSparse installations](#non-standard-suitesparse-installations)
- [Usage](#usage)
  - [staci: hydraulic and transport calculations](#staci-hydraulic-solver)
  - [staci_split: partitioning and sensor placement](#staci_split)
  - [staci_calibrate: pipe diameter calibration](#staci_calibrate)
  - [staci_flush: single-hydrant flushing plans](#staci_flush)
- [Project directories and launching programs](#project-directories-and-launching-programs)
- [Installation and execution troubleshooting](#installation-and-execution-troubleshooting)
- [Tests](#tests)
- [Further documentation and examples](#further-documentation-and-examples)
- [EPANET support, GUI diagnostics, solver settings and validation](#epanet-hydraulic-support-and-current-validation)

## Standalone programs

| Executable | Purpose | Main inputs |
| --- | --- | --- |
| `staci` | Hydraulics, transport, sensitivity, network inspection and conversion | Native `.spr` or EPANET `.inp` network and command-line options |
| `staci_split` | Network partitioning or sensor-placement optimization | XML/JSON settings (`--settings`) and a network |
| `staci_calibrate` | Fit selected pipe diameters to measured pressures and pool levels | XML/JSON settings (`--settings`), period-specific `.spr` networks and measured data |
| `staci_flush` | Evaluate hydrants individually and rank a flushing sequence | EPANET `.inp` network and JSON flushing config |

Each program runs independently. A model-specific project needs only its inputs
and a command or launcher that calls the appropriate executable. The optimizer
tools use pagmo2; `staci_split` also uses Eigen3 and igraph. The build also provides the internal `staci_core` target and C++ examples,
but does not install a separate C++ link library or promise a stable public API.

## Installation

### Quick start and build options

Open a terminal in the repository root (the directory containing this README
and `CMakeLists.txt`). Install the dependencies for your operating system using
the sections below, then build and install:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake --install build --prefix "$HOME/.local"
```

The default build includes all four applications. Add `$HOME/.local/bin` to
`PATH`, or launch executables by their full path. For example, on Linux/macOS:

```sh
export PATH="$HOME/.local/bin:$PATH"
staci --help
staci_split --help
staci_calibrate --help
staci_flush --help
```

This `export` affects the current shell; put it in your shell startup file if
needed permanently. On Windows use the PowerShell/Visual Studio commands below
and the `.exe` files in `build\Release` or the installed `bin` directory.

| CMake option | Default | Effect |
| --- | --- | --- |
| `STACI_BUILD_OPTIMIZERS` | `ON` | Build `staci_split` and `staci_calibrate`; requires pagmo2, Eigen3 and igraph |
| `STACI_ENABLE_HDF5` | `ON` | Look for HDF5 and enable chunked EPS `.h5` output when found |
| `STACI_BUILD_CPP_EXAMPLES` | `ON` | Build the native C++ examples |
| `STACI_BUILD_MATLAB_MEX` | `OFF` | Build the optional MATLAB MEX interface |
| `BUILD_TESTING` | `ON` | Configure CTest checks and test examples; requires Python |
| `STACI_EPANET_LIBRARY` | Unset | Official EPANET shared-library path for independent reference tests |
| `STACI_EPANET_EXECUTABLE` | Unset | Official EPANET executable path for applicable reference tests |

For a smaller build containing `staci` and `staci_flush`:

```sh
cmake -S . -B build-minimal -DCMAKE_BUILD_TYPE=Release -DSTACI_BUILD_OPTIMIZERS=OFF -DSTACI_BUILD_CPP_EXAMPLES=OFF -DBUILD_TESTING=OFF
cmake --build build-minimal --parallel
```

SuiteSparse and nlohmann/json are still required. HDF5 is optional: without it,
EPS still supplies CSV/JSON outputs, but does not produce the chunked `.h5`
file. EPANET is a test/reference dependency, not a required runtime backend for
STACI's own hydraulic solver. Check the CMake configure output to see which
optional features were enabled. Reconfigure after changing options.


### Precompiled executables

The [Build platform binaries](.github/workflows/build-binaries.yml) workflow
builds Release executables on native GitHub-hosted runners for:

- Windows x64;
- Linux x64;
- macOS using the architecture of the current `macos-latest` runner.

Run the workflow manually from the repository's **Actions** tab, or push a tag
matching `v*`. Download the resulting `staci-Windows-*`, `staci-Linux-*`, or
`staci-macOS-*` artifact. See [dist/README.md](dist/README.md) for packaging and
runtime dependency details.

Precompiled binaries are intentionally generated by CI instead of being
committed to Git: executables are platform- and architecture-specific, and
SuiteSparse runtime requirements can change between releases.

### Requirements

- CMake 3.16 or newer;
- a C++17 compiler;
- nlohmann/json 3.2 or newer (shared diagnostics, JSON results and flushing);
- Python 3 for tests (or configure with `-DBUILD_TESTING=OFF`);
- SuiteSparse with UMFPACK, CHOLMOD, AMD, COLAMD, and SuiteSparseConfig.
- pagmo2 and its Boost/TBB dependencies;
- Eigen3 and igraph (used by `staci_split`);
- HDF5 development files (optional for general STACI use, required for the
  chunked EPS `.h5` output).
- MATLAB with Global Optimization Toolbox, CMake-visible Eigen3, and a
  configured supported C/C++ compiler (optional, only for the in-memory MEX
  optimization demonstration).

The optimizer targets are enabled by default. A minimal hydraulic-only build
can omit them with `-DSTACI_BUILD_OPTIMIZERS=OFF`.

CMake first uses SuiteSparse's official imported `SuiteSparse::UMFPACK` target.
For older system packages without CMake metadata, it falls back to finding the
headers and libraries in standard installation prefixes.

### Build on Linux

#### Debian or Ubuntu

```bash
sudo apt update
sudo apt install build-essential cmake libsuitesparse-dev libhdf5-dev \
  libpagmo-dev libeigen3-dev libigraph-dev nlohmann-json3-dev python3

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executables are `build/staci`, `build/staci_flush`, `build/staci_split`
and `build/staci_calibrate`; the last two require optimizers to be enabled.

#### Fedora

This builds the hydraulic and flushing programs. To build the optimizers too,
install pagmo2, Eigen3 and igraph development packages and omit
`-DSTACI_BUILD_OPTIMIZERS=OFF`.

```bash
sudo dnf install cmake gcc-c++ suitesparse-devel hdf5-devel json-devel python3

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTACI_BUILD_OPTIMIZERS=OFF
cmake --build build --parallel
```

To select Clang explicitly:

```bash
CC=clang CXX=clang++ cmake -S . -B build-clang -DCMAKE_BUILD_TYPE=Release
cmake --build build-clang --parallel
```

### Build on macOS

Install the command-line tools and dependencies:

```bash
xcode-select --install
brew install cmake suite-sparse hdf5 pagmo eigen igraph nlohmann-json python
```

Configure and build with Apple Clang:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executables are `build/staci`, `build/staci_flush`, `build/staci_split`
and `build/staci_calibrate`; the last two require optimizers to be enabled.

To use Homebrew GCC instead, adjust the version suffix to the installed
compiler:

```bash
CC=gcc-16 CXX=g++-16 cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release
cmake --build build-gcc --parallel
```

### Build on Windows

The recommended setup is Visual Studio 2022 with the **Desktop development
with C++** workload, CMake, Git, and vcpkg.

In PowerShell:

```powershell
vcpkg install suitesparse-umfpack:x64-windows hdf5:x64-windows pagmo2:x64-windows `
  eigen3:x64-windows igraph:x64-windows nlohmann-json:x64-windows

cmake -S . -B build `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"

cmake --build build --config Release --parallel
```

Install `suitesparse-umfpack` explicitly: the `suitesparse` umbrella package
does not include UMFPACK by default. Its required SuiteSparse dependencies are
installed automatically.

If `VCPKG_ROOT` is not defined, replace it with the absolute path to the vcpkg
checkout. The executables are normally in `build\Release\`, with `.exe` suffixes.

To build with Ninja instead of the Visual Studio generator:

```powershell
cmake -S . -B build -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"

cmake --build build --parallel
```

### Install to a staging directory

CMake can copy the built executables into a conventional `bin` directory:

```bash
cmake --install build --prefix install
```

For a multi-configuration generator such as Visual Studio:

```powershell
cmake --install build --config Release --prefix install
```

The installed executables are `staci`, `staci_flush`, `staci_calibrate`, and
`staci_split` (with `.exe` suffixes on Windows). The two optimizer programs are
omitted when `STACI_BUILD_OPTIMIZERS=OFF`. Installation copies the enabled
application executables; it does not package every shared dependency, project
settings file or example. Keep the required runtime libraries available, and
prepare project inputs separately. A staging directory can be removed to
uninstall its copied executables; remove its `bin` entry from `PATH` if present.

### Non-standard SuiteSparse installations

For a custom SuiteSparse prefix, pass it through CMake's standard search path:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/suitesparse
```

The fallback finder also recognizes the `SUITESPARSE_ROOT` environment
variable. Ensure the matching SuiteSparse shared libraries are on the runtime
library search path when using a dynamically linked build.

## Usage

The examples below run from the repository root with a single-configuration
build in `build/`, unless another working directory is shown. On Windows with
Visual Studio, use `build\Release\PROGRAM.exe` in place of `build/PROGRAM`.
Installed executables can be called directly when their `bin` directory is on
`PATH`.

### staci hydraulic solver

Use `staci` to solve steady hydraulics, calculate transport and sensitivities,
inspect or modify network properties, and convert networks. It also supports
EPANET extended-period hydraulics and steady water-quality calculations.

```sh
./build/staci --help

# Copy a sample: a steady solve writes results back into the native SPR file.
cp tests/anytown_1med.spr network.spr
./build/staci -s network.spr

# List nodes and edges, then export the network to EPANET.
./build/staci -l network.spr
./build/staci --export-epanet network.spr -o network.inp
```

Common commands:

| Command | Operation |
| --- | --- |
| `staci -s FILE` | Solve steady hydraulics |
| `staci -l FILE` | List elements; also write `element_list.txt` |
| `staci -g FILE -e ID -p PROPERTY` | Read a property; write its value to `tmp.dat` |
| `staci -m FILE -e ID -p PROPERTY -n VALUE -o NEW.spr` | Change a property and save a new network |
| `staci --export-epanet FILE -o NEW.inp` | Export an EPANET network |
| `staci --epanet-eps FILE.inp -o PREFIX` | Run extended-period hydraulics; write HDF5, metadata and CSV results |
| `staci -t FILE` | Calculate residence time |
| `staci -c FILE` | Calculate concentration transport |
| `staci -r FILE -e ID -p PROPERTY` | Calculate parameter sensitivity |
| `staci --steady-quality FILE.inp -o PREFIX` | Solve asymptotic water quality for fixed hydraulics |

An INP initial-state solve and an EPS run can be launched as follows:

```sh
./build/staci -s tests/epanet_eps_smoke.inp
./build/staci --epanet-eps tests/epanet_eps_smoke.inp -o /absolute/path/results/smoke
```

Create the output parent directory first. INP hydraulic solves preserve input
bytes and write SI results in `<network>.inp.hydraulics.json`. Native SPR
solves can update results inside the input file; work on a copy when preserving
an original model. EPS uses a result prefix and writes CSV/JSON metadata plus
HDF5 when enabled. Use network IDs exactly as listed by `-l`. Property units
are command/property-specific: native demand properties use m³/h, while
machine-readable hydraulic exports use SI m³/s.

Hydraulic runs produce `.ros` logs, `.rps` progress files and `.rrs` convergence
markers beside the input. See the [full usage reference](doc/staci_usage.md) for
property names and units, transport settings, EPS controls, quality calculations,
EPANET compatibility, all options and output files.

### staci_split

Use `staci_split` to divide a network into communities or select sensor locations.
Use **`--settings file.xml` or `--settings file.json`** to choose settings; the
network path is set by `fname` inside that file. The optional `--seed` argument
makes the optimizer's random choices reproducible. Without `--settings`, it
loads `staci_split_settings.xml`, or `staci_split_settings.json` if XML is absent.
See [JSON input formats](doc/json_inputs.md) for all auxiliary input schemas.

```sh
# Run inside a project directory containing staci_split_settings.xml.
/path/to/staci_split --seed 12345
```

A small partitioning config is:

```xml
<settings>
  <global_debug_level>1</global_debug_level>
  <n_comm>2</n_comm>
  <weight_type>topology</weight_type>
  <weight_type_mod>diameter</weight_type_mod>
  <fname>network.spr</fname>
  <logfilename>split.log</logfilename>
  <obj_type>modularity</obj_type>
  <popsize>4</popsize>
  <ngen>2</ngen>
  <pmut>0.25</pmut>
  <pcross>0.8</pcross>
</settings>
```

`modularity` partitions the graph into `n_comm` communities. Its weights can be
`topology`, `dp` or `sensitivity`; sensitivity weighting additionally uses
`weight_type_mod`. The alternative objectives `A-optimality` and `D-optimality`
select sensor locations using sensitivity-based criteria, with `weight_type` set
to `friction_coeff`, `diameter` or `demand`. The population size,
generation count and mutation/crossover probabilities above are demonstration
settings, not a convergence recommendation.

Results include `membership.txt` and the configured optimization log;
sensitivity modes also write matrix files. Use a separate working
directory for each project because these output names are reused.

With optimizers and testing enabled, CMake prepares a runnable example:

```sh
(cd build/optimizer-tests/split && ../../staci_split --seed 12345)
```

See the [settings template](tests/staci_split_settings.xml.in) and
[splitter manual (Hungarian)](doc/staci_split_manual.md) for further details.
CMake resolves the `@...@` paths in the template when generating the example.

### staci_calibrate

Use `staci_calibrate` to fit selected **pipe diameters** to measured node pressures
and pool levels over one or more periods. Each active diameter varies between
0.5 and 2 times its original value. It reads
**`--settings file.xml` or `--settings file.json`**; without this option it searches the current directory for `staci_calibrate_settings.xml`, then `.json`.

```sh
# Run inside a project directory containing staci_calibrate_settings.xml.
/path/to/staci_calibrate --seed 12345
```

A minimal single-period example is:

```xml
<settings>
  <global_debug_level>1</global_debug_level>
  <Staci_debug_level>0</Staci_debug_level>
  <dir_name>./</dir_name>
  <fname_prefix>calibrate_network_</fname_prefix>
  <logfilename>calibrate.log</logfilename>
  <best_logfilename>calibrate-best.log</best_logfilename>
  <sollwert_dfile>calibrate_targets.csv</sollwert_dfile>
  <Start_of_Periods>0</Start_of_Periods>
  <Num_of_Periods>1</Num_of_Periods>
  <Spoil_Active_Pipes>no</Spoil_Active_Pipes>
  <dt>1.0</dt>
  <weight_p_err>1.0</weight_p_err>
  <type_of_pipe_selection>largest_diameter</type_of_pipe_selection>
  <num_of_active_pipes>1</num_of_active_pipes>
  <popsize>4</popsize>
  <ngen>1</ngen>
  <pmut>0.2</pmut>
  <pcross>0.8</pcross>
</settings>
```

This loads `./calibrate_network_0.spr` and `./calibrate_targets.csv`. Subsequent
periods increment the network filename suffix. Keep a trailing separator in
`dir_name`; `dt` is the period length in hours. Target rows have the form
`NODE_ID;node;VALUE;` or `POOL_ID;pool;VALUE;`, with one value per period and a
trailing semicolon. Node pressures are in bar and pool levels in metres.

`weight_p_err` weights the pressure error against the pool-level error (1 selects
pressure only). Active pipes can be selected by `largest_diameter`,
`most_sensitive`, or `Dmin`; the first two use `num_of_active_pipes`, while `Dmin`
uses a diameter cutoff named `Dmin`. Use larger optimization runs and assess
convergence for a real calibration.

Outputs include `best.dat` with fitted diameter ratios and target comparisons,
`bog.dat` with optimization progress, the configured logs, and final sensitivity
files. Treat these as calibration results to review; the program does not
provide a command-line option to save a calibrated network under a new filename.

With optimizers and testing enabled, run the generated example:

```sh
(cd build/optimizer-tests/calibrate && ../../staci_calibrate --seed 12345)
```

See the [settings template](tests/staci_calibrate_settings.xml.in) and
[sample targets](tests/calibrate_targets.csv). JSON settings and measurement examples are in
[examples/config](examples/config). Select a configuration with `--settings`; relative paths inside it use the working directory.

### staci_flush

Use `staci_flush` to assess opening **one hydrant at a time**, without changing
valve positions, and construct a flushing sequence. Inputs are an EPANET network
and a JSON config containing the hydrant junction IDs and global outlet settings.

```sh
./build/staci_flush --inp /path/to/network.inp --config /path/to/flushing_config.json
```

Example `flushing_config.json`:

```json
{
  "hydrant_node_ids": ["J1", "J2"],
  "hydrant_area_m2": 0.002,
  "total_loss_coefficient": 2.0,
  "velocity_threshold_mps": 0.5,
  "min_pressure_head_m": 0.0,
  "output_dir": "results/flushing",
  "write_network_files": true
}
```

Use exact junction IDs from the INP file. Area is in m², velocity in m/s and
pressure head in m. The total coefficient includes outlet kinetic head:
`K = 1 + zeta` for a local loss coefficient `zeta`. Discharge is calculated from
junction pressure using `Q = A * sqrt(2*g*max(h,0)/K)`.
`output_dir` is relative to the config file and must be new or empty.

For each valid scenario, pipes qualify when their absolute velocity is strictly
above the threshold. The plan first selects the hydrant covering the largest
pipe volume, then repeatedly selects the largest **additional uncovered volume**.
It reports cumulative coverage as a percentage of all original pipe volume.
Each hydrant has two time estimates: the longest directed travel time from a
qualifying pipe to the hydrant using actual velocities (including slower
connecting pipes), and total qualifying volume divided by hydrant discharge.
Undetermined travel times are flagged. Plan volumes use two decimal places and
times use minutes with one decimal place.

Main results are `flushing_plan.txt`, `flushing_plan.csv`, `summary.txt` and
per-pipe/scenario CSV tables. With `write_network_files: true`, exported networks
are named `<original_name>_hydrant_<nodeID>.inp`, contain the open-hydrant emitter,
and tag qualifying pipes `FLUSHED` for inspection. These are independent steady
scenarios; the tool does not simulate sediment removal.

Run the complete [worked example](examples/flushing/README.md):

```sh
./build/staci_flush --inp examples/flushing/network.inp --config examples/flushing/flushing_config.json
```

Without a config file, a junction list can be supplied with `--hydrants`
using a text or JSON file and explicit outlet parameters:

```sh
./build/staci_flush --inp network.inp --hydrants hydrants.txt --hydrant-area-m2 0.002 --loss-coefficient 2 --velocity-threshold-mps 0.5 --min-pressure-head-m 0 --output-dir /absolute/path/new-flushing-results
```

The destination must be new or empty. Use `--help` for the exact supported
argument names. All four executables accept the common diagnostics and solver
options described below.

See the [flushing reference](doc/flushing.md) for the full output schema,
pressure filtering, timing assumptions, alternative inputs and EPANET inspection.

## Project directories and launching programs

Use a separate working directory for each model, optimization or flushing job.
The optimizer settings filenames are fixed and selected by the process working
directory; several legacy outputs also have fixed names. A typical layout is:

```text
project/
  network.inp                  # or a copied native network.spr
  staci_split_settings.xml     # only for a split job
  staci_calibrate_settings.xml # only for a calibration job
  calibrate_network_0.spr      # calibration period 0, when configured
  calibrate_targets.csv        # measured targets, when configured
  flushing_config.json        # only for a flushing job
  results/
```

| Program | How inputs are selected | Main outputs |
| --- | --- | --- |
| `staci` | Network path and mode on the command line | Hydraulic JSON for INP, native results, `.ros`/`.rps`/`.rrs`; EPS/quality files according to mode |
| `staci_split` | XML/JSON settings selected by `--settings` (or the default filename); `fname` selects the network | `membership.txt`, configured log and sensitivity matrices when applicable |
| `staci_calibrate` | XML/JSON settings selected by `--settings` (or the default filename), period-specific SPR files and measured targets | `best.dat`, `bog.dat`, configured logs and sensitivity files |
| `staci_flush` | `--inp` plus JSON config or explicit hydrant/outlet arguments | Plan TXT/CSV, summary, per-pipe/scenario tables and optional inspection INP files |

`staci_split` and `staci_calibrate` accept `--settings path.xml|path.json`.
Without it, the program uses `staci_<name>_settings.xml`, or the `.json` sibling
if XML is absent. Use the example configurations, replace
all model paths and IDs, then launch from the project directory. CMake's `.in`
templates contain substitution placeholders; use the generated examples or
replace the placeholders before running. Calibration currently assembles
period filenames with `.spr`, so shared INP support in the hydraulic core is
not a promise that this calibration input convention accepts INP periods.

For a GUI, configure the executable path, working directory and arguments as
separate process-launch fields. Supply the same absolute diagnostics path to
all applications and a separate result directory/prefix per job. For example:

```sh
cd /absolute/path/project
/absolute/path/bin/staci_split --seed 12345 --max-iterations 100 --diagnostics-file /absolute/path/project/results/staci-diagnostics.jsonl
```

Wait for process completion, inspect its exit code and the matching JSONL
`run_end`, then load the output files. Warnings can accompany successful
results. Exit code 3 indicates incomplete flushing/EPS results that require
inspection. Retain diagnostics when presenting failures to the user. The
complete protocol and exit-code table are included later in this README.

## Installation and execution troubleshooting

| Symptom | Check/action |
| --- | --- |
| Executable not found | Use its full path, or add the installed `bin` directory to `PATH`; on Windows check `build\Release` and `.exe` |
| Split/calibrate executable absent | Configure with `STACI_BUILD_OPTIMIZERS=ON` and install the optimizer dependencies |
| CMake cannot find a package | Install development packages; use the matching compiler architecture and `CMAKE_PREFIX_PATH`/package paths where needed |
| Shared-library/DLL load failure | Keep the required runtime libraries installed; CMake installation alone does not bundle them all |
| Settings not found | Supply `--settings path.xml` or `--settings path.json`, or use the default filename in the working directory |
| Missing model/measurement file | Verify settings paths, period suffixes and IDs; use absolute paths when launched by a GUI |
| Flushing refuses its output directory | Select a new or empty directory |
| EPS does not produce HDF5 | Check configure output for HDF5 detection; enable `STACI_ENABLE_HDF5` and install its development package |
| Solver does not converge | Read the common diagnostics, worst residuals and network-specific `.ros` log; check physical feasibility before changing limits |
| Reference tests skip | Provision official EPANET and set the library/executable paths; use `--require-reference` when comparison is mandatory |

On macOS, a compiler/SDK selection problem should be checked with
`xcode-select -p`; use a valid installed Command Line Tools or Xcode toolchain.
For an Xcode installation, configure/build with the appropriate
`DEVELOPER_DIR` if necessary. If CMake still selects a different SDK, set
`-DCMAKE_OSX_SYSROOT="$(DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer xcrun --sdk macosx --show-sdk-path)"`
at configuration time, with the actual installed Xcode path. Build and runtime architecture must match the
installed dependency architecture.

## Tests

The [2026-10-02 package audit](doc/package_validation.md) records the complete
CTest, portable-runner, MATLAB and installation checks, along with known limits.

Build first, then run:

```sh
ctest --test-dir build --output-on-failure
# Only the flushing regression tests:
ctest --test-dir build -R flushing --output-on-failure
```

Add `-C Release` for Visual Studio builds. The flushing tests cover the outlet
law, ranking and timing, input validation, network exports and the worked example.
See the [testing reference](doc/testing.md) for the full suite, optional official
EPANET comparisons and channel reference cases.
The suite also includes [97 public EPANET networks](tests/public_networks/README.md)
for standalone hydraulic snapshots and explicit compatibility diagnostics.

## Further documentation and examples

- [Developer integration guide: web backends, Python, MATLAB and error handling](README_DEV.md)

- [Full staci usage and command-line reference](doc/staci_usage.md)
- [Flushing configuration, planning and output reference](doc/flushing.md)
- [Worked flushing example with expected results](examples/flushing/README.md)
- [Testing and reference comparisons](doc/testing.md)
- [Native C++ examples](examples/cpp/README.md)
- [Examples and MATLAB demonstrations](examples/README.md)
- [Binary packaging and runtime dependencies](dist/README.md)

## EPANET hydraulic support and current validation

This section collects the implementation, GUI integration, numerical settings,
public-network validation and remaining issues in one place. Results below are
from the 2026-10-01 validation on macOS. The shared hydraulic core is used by
`staci`, `staci_split`, `staci_calibrate` and `staci_flush`; the public-network
corpus exercises standalone `staci`.

### Supported hydraulic models

| Feature | Implementation and behavior |
| --- | --- |
| Pipes | EPANET H-W and D-W equations, minor losses, viscosity/transition treatment, low-flow regularization and check valves |
| HEAD pumps | Head curves, operating speed/patterns and reverse-flow closure |
| POWER pumps | Constant hydraulic power, EPANET unit conversion and cubic speed scaling |
| Reservoirs and tanks | Fixed boundary at the initial level in steady mode; tank storage/levels and limits in EPS |
| TCV | Throttle valve with setting and minor-loss coefficient |
| PRV | Regulates downstream pressure, opens when supply pressure is insufficient and closes against reverse flow |
| PSV | Regulates upstream pressure, opens when downstream pressure prevents regulation and closes against reverse flow |
| PBV | Imposes a head drop; uses open-valve minor loss if it exceeds the setting |
| FCV | Regulates forward flow; opens and warns when the requested flow cannot be maintained |
| GPV | Named flow/headloss curve with piecewise interpolation and endpoint extrapolation |
| Emitters | New single-junction atmospheric outlet element with pressure-dependent discharge |
| PDA | Pressure-dependent positive demand, integrated into the continuity equation and its Jacobian |

PSV and PBV extend the existing `EpanetValve` class, which derives from
`JelleggorbesFojtas`. All control valves share the existing setting, status,
control, EPS and export interfaces. Explicit OPEN/CLOSED overrides are supported.
Pressure settings include input pressure-unit and specific-gravity conversion.
GPV settings identify curves; arbitrary numeric GPV settings are rejected.
Malformed curves and invalid coefficients produce informative diagnostics.

`EpanetEmitter` derives from `Agelem`. Its signed discharge law is
`Q = C * sign(p) * abs(p)^exponent`, where `p` is junction pressure head and the
coefficient is converted to SI. The inverse relation and analytic derivative
supply the branch equation. Negative-pressure atmospheric inflow reproduces
EPANET behavior and generates `EPANET.EMITTER_BACKFLOW`; the model author should
check whether that inflow is physical for the intended network.

PDA delivers zero positive demand below minimum pressure, full demand at or
above required pressure, and a fraction
`((p - pmin) / (preq - pmin))^exponent` between them. Negative demands remain
fixed injections. Demand patterns update nominal demand. Hydraulic outputs
report delivered demand plus emitter discharge. Exponents must be positive;
required pressure must exceed minimum pressure. Input validation follows
EPANET's minimum 0.1 difference in input pressure units.

For the new pressure-dependent models, the solver initializes emitter flow
consistently with pressure, rebuilds the nonlinear Jacobian and backtracks
Newton steps. These steps prevent running POWER pumps from crossing zero flow
and selecting an unphysical reverse-flow root. Final convergence tolerances
and explicit iteration limits remain unchanged.

Initial-state demand and pump patterns, supported initial controls and tank
limits are applied before steady comparison. Time-dependent tank behavior,
patterns and controls require EPS mode. Importing a network for steady mode
does not establish equivalence of its entire time-dependent simulation or
water-quality calculation.

### Common GUI diagnostics

All four applications append warnings and errors to the same JSON Lines file.
The default is `staci-diagnostics.jsonl` in the process working directory. A GUI
should pass an absolute `--diagnostics-file PATH` on every invocation. The
`STACI_DIAGNOSTICS_FILE` environment variable is an alternative; the CLI option
has precedence.

Each line is one JSON object with `schema_version`, UTC `timestamp`, `program`,
`run_id`, `sequence`, `event`, `severity`, `code` and `message`. Where available,
import errors also identify the network, section, element and input line. A run
has `run_start` and `run_end`; the latter includes `exit_code`, `error_count`
and `warning_count`. Group records by `run_id`, show the human-readable
messages and check the process exit code. A killed process may have no
`run_end`. Writes are locked across processes; complete records can interleave.
Consumers should tolerate additional fields.

| Exit code | Meaning |
| --- | --- |
| 0 | Successful operation, possibly with warnings |
| 1 | Calculation or execution failure |
| 2 | Invalid input/configuration or unsupported feature |
| 3 | Partial result: some frames or scenarios failed |

Errors also appear on stderr. A diagnostics-file open failure reports
`DIAGNOSTICS_OPEN` and exits with code 1. Numerical failures of discarded
calibration candidates are warnings; failure of the initial calibration
baseline is an error.

Useful hydraulic codes include `EPANET.COMPATIBILITY`,
`HYDRAULICS.NONCONVERGENCE`, `HYDRAULICS.ILL_CONDITIONED`,
`EPANET.FCV_UNATTAINABLE`, `EPANET.NO_AVAILABLE_SUPPLY`,
`EPANET.EMITTER_BACKFLOW` and `EPANET.POWER_PUMP_DEAD_END`. Nonconvergence
messages include RMS residuals with units/limits and the worst node/link.
Unsupported or malformed physical definitions are rejected instead of being
silently dropped from the network.

Successful steady INP runs write `<network>.inp.hydraulics.json`, containing
convergence and SI node/link results. A new steady solve removes a previous
hydraulic export before validating the input, preventing stale GUI results.
Existing console output and network-specific result files remain available.
See [the full diagnostic protocol](doc/diagnostics.md).

### Solver settings

All four applications accept these common options, including `--option=value`:

| Option | Meaning |
| --- | --- |
| `--head-tolerance-m VALUE` | RMS head-equation residual limit in metres |
| `--mass-tolerance-kg-s VALUE` | RMS junction-continuity residual limit in kg/s |
| `--max-iterations INTEGER` | Maximum nonlinear iterations |
| `--diagnostics-file PATH` | Shared machine-readable warning/error log |

CLI settings override input/application settings. Numeric values must be
positive and finite, and iteration counts must be integers. Invalid values
produce `CLI_ARGUMENT` and exit code 2.

The INP defaults are **0.0001 m (0.1 mm)** RMS head residual and **1e-8 kg/s**
RMS continuity residual. The iteration limit comes from the input setting;
without a specified INP trial limit it defaults to 100. Native SPR settings
remain in effect unless overridden. RMS convergence is not a bound on the
largest individual residual or on the difference from EPANET: nearly stagnant
branches can amplify small head errors into flow differences.

```sh
staci -s network.inp --head-tolerance-m 0.0001 --mass-tolerance-kg-s 1e-8 --max-iterations 100 --diagnostics-file /absolute/path/staci-diagnostics.jsonl
staci --epanet-eps ky16.inp -o ky16 --head-tolerance-m 0.0001 --max-iterations 100
```

The ky16 regression passes all 25 EPS report frames at the default 0.1 mm
setting. A failure at much tighter precision is a separate numerical issue.

### Public networks and independent comparison

The corpus contains **97 retained INP files**: 96 unchanged upstream files and
one intentionally adapted physical Anytown model. Inputs span small targeted
fixtures through large networks, including a 57,460-pipe case. Sources,
revisions, hashes, element counts and expected outcomes are recorded in
[the manifest](tests/public_networks/manifest.json). Tests verify input hashes
and do not rewrite the corpus inputs during execution.

The independent reference is the official OWA EPANET 2.2 C API, pinned to
revision `4d8d82ddc260fce216af9321fc3d9a4646ac6827`. It supplies double-precision
initial-state (`t = 0`) results, avoiding text-report rounding. Reference
settings use 1000 trials and relative accuracy 1e-6, with a fresh 1e-5 retry
for unreliable convergence. Nonfinite results and failed independent junction
continuity audits are rejected. An unavailable reference is a skip, never an
agreement pass; `--require-reference` makes availability mandatory.

The comparison checks heads, junction pressure, delivered demand, source/tank
exchange, signed flow, velocity, endpoint head difference and enabled state.
Limits are common to every model, not tuned per network:

| Quantity | Absolute tolerance | Relative tolerance |
| --- | --- | --- |
| Total head / junction pressure head | 0.01 m | 1e-6 |
| Junction demand | 1e-8 m³/s | 1e-4 |
| Reservoir/tank net exchange | 1e-6 m³/s | 1e-4 |
| Signed link flow | 1e-6 m³/s | 1e-4 |
| Velocity magnitude | max(1e-5 m/s, absolute flow limit / area) | 1e-4 |
| Endpoint head difference | 0.02 m | 1e-5 |
| Enabled/closed state | Exact | Exact |

A numeric value passes against the larger of the absolute and relative limits.
The 1e-12 m head / 1e-8 kg/s continuity profile used by the standard reference
corpus is deliberately stricter than the normal runtime default.

Current initial-state results: **88 numerical matches** across 97 retained
networks; three inputs are rejected by EPANET and six references are physically
unreliable. Standalone strict execution is 89 solved, four compatibility/input
rejections and four hydraulic failures. Expected diagnostic cases are regression
passes, not solved/equivalent networks. See [current validation](doc/epanet_reference.md)
for the full scope and reproducible commands.

Net1 and Net2 now pass full-period hydraulic and chemical comparison, including
MIXED tank storage, first-order reactions and wall mass transfer. Full-corpus
EPS equivalence remains incomplete. Unsupported chemical mixing/reaction laws
produce explicit errors; steady hydraulics remains available for those inputs.

Seven formerly rejected models now solve and match EPANET: `Net6_plus`,
`epanet_leaks`, `NET1emit`, `NET1negemit`,
`psv_open_no_downstream_sources`, `NET1-PBV` and `cheung`.

### Remaining issues and model adaptations

- **io.inp:** all its element types are supported. Running POWER pump `pump2`
  feeds terminal zero-demand junction `j4`, with no outlet. Continuity requires
  zero flow, incompatible with finite head at nonzero constant power.
  `EPANET.POWER_PUMP_DEAD_END` explains this and suggests stopping the pump or
  providing a physical discharge path. The input was not changed.
- **JEP5-13.inp:** negative PRV settings are supported and the reference matches.
- **bad_syntax.inp / bad_values.inp:** retained negative fixtures for an unknown
  section and an unknown endpoint.
- **conditional_controls_2.inp / control_comb.inp:** initial source/connectivity
  compatibility failures; their EPANET references are also unreliable.
- **NW_Model.inp / NW_Model1.inp / GES4-9.inp:** now match under the strict reference profile.
- **fcv_open_no_downstream_sources.inp / prv_closed_no_upstream_sources.inp /
  2fcvs.inp:** physical regulation/supply conflicts, now identified by
  `HYDRAULICS.VALVE_CONSTRAINT`; their EPANET references are also unreliable.
- **cv_controls.inp:** solves in STACI, but EPANET rejects the fixture under its
  minimum junction-count requirement, so no independent agreement is claimed.
- **ky10.inp:** now matches the reliable reference operating point; the earlier
  false-convergence investigation is retained in `doc/ky10_convergence.md`.

Anytown was intentionally converted into a physical model with user-approved
input changes. All three pump speeds and their 24-value speed patterns are
1.0. Six 0.0001-inch candidate pipes (`110`, `113`, `114`, `115`, `116`, `125`)
were resized to 16 inches, trunks `1`, `2`, `3` to 24 inches, and the other
8/10/12/16-inch distribution pipes to 12/16/18/24 inches respectively; existing
30-inch pipes were retained. Topology, demands, demand patterns, pump curves,
roughness, lengths, elevations, tank definitions and time steps were preserved.
Both the original and adapted hashes are recorded, and the original is retained
at [doc/model_sources/Anytown_upstream.inp](doc/model_sources/Anytown_upstream.inp).
The adapted model converges at normal/strict tolerances and matches the initial
EPANET state. Its 24-hour run passes physical pressure, velocity and tank-limit
checks; this is not a claim of full EPS numerical equivalence.

Zero-speed HEAD and POWER pumps are supported and tested for direct `SPEED 0`,
zero speed patterns and numeric status settings, including EPS stop/start
transitions. The Anytown and zero-speed regressions are part of
`epanet_physical_models`.

### Reproduce and inspect results

```sh
python3 tests/setup_epanet_reference.py
cmake -S . -B build -DSTACI_EPANET_LIBRARY=/absolute/path/libepanet2.so -DSTACI_EPANET_EXECUTABLE=/absolute/path/runepanet
cmake --build build
ctest --test-dir build --output-on-failure
python3 tests/test_public_networks.py --binary build/staci --reference-library /absolute/path/libepanet2.so --require-reference --output-dir build/public-network-reference-tests
python3 tests/test_epanet_pressure_elements.py --binary build/staci --reference-library /absolute/path/libepanet2.so
```

Use `libepanet2.dylib` on macOS or `epanet2.dll` on Windows, with the actual
paths printed by setup. Multi-configuration builds need the appropriate
configuration and executable subdirectory. The standard suite includes
`application_diagnostics`, `solver_options`, `public_epanet_networks`,
`public_epanet_reference`, `epanet_reference_units`, `epanet_reference_valves`,
`epanet_pressure_elements` and `epanet_physical_models`.

Saved results from this development are in
`tests/test-results/pressure-elements/`: `summary.json`, `ctest.log`,
`targeted-tests.log`, `public-execution/report.json` and
`public-reference/report.json`, with per-network diagnostic logs. Reported
errors include the affected asset and measured differences. Reports contain
local execution paths; reproduce the commands above on another machine.

Detailed implementation notes remain available in
[pressure elements](doc/epanet_pressure_elements.md),
[EPANET valves](doc/epanet_valves.md),
[reference comparisons](doc/epanet_reference.md),
[error classification](doc/error_classification.md),
[Anytown adaptation](doc/anytown_physical_model.md) and
[the test reference](doc/testing.md).
