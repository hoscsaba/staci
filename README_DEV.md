# STACI developer integration guide

This guide describes the interfaces actually present in this repository. Use the
command-line programs from a web backend, Python, or any GUI. MATLAB also has an
in-memory MEX interface. The repository does not expose an HTTP server, a stable
C ABI, a Python extension module or a separate `staci_eps` executable.

## Build and repository layout

CMake 3.16+, a C++17 compiler, nlohmann-json and a supported sparse linear solver
are required. See [README](README.md#installation) for platform dependencies.
The optional optimizers additionally need pagmo2, Eigen3 and igraph. HDF5 enables
binary EPS output; CSV and JSON remain available without it.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTACI_BUILD_OPTIMIZERS=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Enable `STACI_BUILD_OPTIMIZERS=ON` for all four applications. MATLAB is optional:
`STACI_BUILD_MATLAB_MEX=ON` requires an installed MATLAB and compatible compiler.
For a multi-configuration generator use `cmake --build build --config Release`
and binaries in `build/Release/`. Installed binaries retain external shared
library dependencies; package these or provision them on the target machine.

| Path | Purpose |
|---|---|
| `src/Staci.cpp`, `include/staci/Staci.h` | Hydraulic assembly, nonlinear solve, CLI operations |
| `src/epanet_reader.cpp`, `epanet_writer.cpp` | EPANET import/export and SI conversion |
| `src/epanet_extended_simulation.cpp` | EPS events, controls, tank levels, quality configuration |
| `src/epanet_water_age.cpp` | Transient age and chemical transport |
| `src/diagnostics.cpp` | Shared CLI error boundary and JSONL logging |
| `include/staci/StaciSession.h` | In-memory C++ facade |
| `matlab/` | MEX dispatcher and `StaciModel` handle class |
| `tests/public_networks/manifest.json` | Pinned corpus, hashes and expected outcomes |

## Process contract

Give every job a unique working directory. Copy inputs there before running:
INP files are preserved, but native SPR operations can update the input, and
some programs use fixed result filenames. Create output parent directories first.
Use an absolute executable path and an argument array, not a shell command
assembled from an uploaded filename. Parallel jobs use separate processes and
separate input/configuration copies.

| Application | Invocation and inputs | Main result |
|---|---|---|
| `staci` | `-s network.inp` | `network.inp.hydraulics.json`, legacy `.ros/.rps/.rrs` files |
| `staci` EPS | `--epanet-eps network.inp -o result` | `result.meta.json`, `result-{nodes,links,tanks,summary}.csv`, optional `result.h5` |
| `staci` steady quality | `--steady-quality network.inp -o result` | `result-steady-{nodes,links,summary}.csv` |
| `staci_split` | `--settings FILE.xml` or `FILE.json`; default `staci_split_settings.xml` then `.json`; optional `--seed N` | `membership.txt` and configured logs |
| `staci_calibrate` | `--settings FILE.xml` or `FILE.json`; default `staci_calibrate_settings.xml` then `.json`; optional `--seed N` | Configured logs, fitted network and history files |
| `staci_flush` | See `staci_flush --help` and [flushing API](doc/flushing.md) for plan/options | Plan/scenario results and diagnostics |

All four applications accept `--diagnostics-file /absolute/path/errors.jsonl`.
Otherwise `STACI_DIAGNOSTICS_FILE` supplies the path, falling back to
`staci-diagnostics.jsonl` in the working directory. The log is appended to;
parallel writes are locked. Records can interleave, so group by `run_id`.
A fresh per-job log is simplest for services. Keep console output in a separate
file: stdout is human-oriented and is not a stable result protocol.

```sh
/opt/staci/bin/staci --diagnostics-file /srv/jobs/123/diagnostics.jsonl \
  --epanet-eps /srv/jobs/123/network.inp -o /srv/jobs/123/result \
  --head-tolerance-m 0.0001 --mass-tolerance-kg-s 1e-8 --max-iterations 1000
```

The INP default head-equation tolerance is 0.1 mm **RMS residual**, not a bound
on every head value or on the difference from EPANET. Explicit SPR settings are
retained unless overridden. Reference tests use a stricter profile.

Auxiliary files also have JSON variants: optimizer settings, calibration measurements,
`staci -i` initial values, flushing configuration and hydrant lists. See
[JSON input formats](doc/json_inputs.md) for schemas, units and examples. Network
definitions remain SPR/XML or EPANET INP.

Calibration CSV measurements are checked as finite numbers in the selected
period window. Invalid cells return exit code 2 and `INPUT.CONFIG` with the
file, measurement row, element ID and zero-based period; clients must not treat
these runs as fitted results. See the [measurement contract](doc/json_inputs.md#calibration-measurements)
for the multi-period pool-state update convention.

`staci_flush` accepts `"mode": "single"` (default) or `"mode": "multi"` in its
JSON config (`MODE` is an alias). Multi opens all configured nodes in one solve;
the aggregate report has one scenario row, while `scenario_hydrants.csv` holds
individual outlet values. See the [flushing output contract](doc/flushing.md)
for aggregate pressure/flow semantics and route diagnostics.

## Completion and errors

| Process exit code | Meaning | Consumer behavior |
|---|---|---|
| 0 | Success, possibly warnings | Read results and display warnings |
| 1 | Calculation/execution failure | Display diagnostics; do not use output as a successful result |
| 2 | Invalid input/configuration or unsupported feature | Display the offending input/feature |
| 3 | Partial EPS/scenario result | Mark partial; inspect convergence flags for individual frames |

A timeout or signal termination is a process-manager outcome, not STACI exit
code 3. A killed process may lack `run_end`; treat such a job as incomplete.
Do not infer success from the existence of a result file.

JSONL records have `schema_version: 1`, `timestamp` (UTC), `program`, `run_id`,
`sequence`, `event`, `severity`, `code`, and human-readable `message`. Optional
fields include `network`, `section`, `element` and `line`. Consumers must tolerate
new fields/codes. `run_start` begins a run and `run_end` includes `exit_code`,
`error_count`, and `warning_count`. Require a matching completion record and
process exit code. If the log cannot be opened, stderr reports `DIAGNOSTICS_OPEN`
and the process exits 1 without a usable log.

Useful codes include:

- `CLI_ARGUMENT`, `CLI.INPUT`: invalid arguments or missing input.
- `INPUT.XML`, `EPANET.COMPATIBILITY`: input/compatibility errors.
- `HYDRAULICS.NONCONVERGENCE`: residuals, limits, worst node and worst link.
- `HYDRAULICS.VALVE_CONSTRAINT`: a closed valve carrying significant flow or an
  active FCV whose computed flow violates its setpoint; a numerical regularizer
  cannot substitute for a physical supply path.
- `EPANET.POWER_PUMP_DEAD_END`: a running constant-power pump has no finite
  working point because its terminal outlet requires zero flow.
- `EPANET.FCV_UNATTAINABLE`: fully open FCV cannot attain the requested flow.
- `EPANET.QUALITY_MIXING`, `EPANET.QUALITY_REACTION`: unsupported chemical EPS
  mixing/reaction law; no falsely successful chemical result is produced.
- `EPANET.EPS_PARTIAL_FAILURE`: some hydraulic states failed.
- `FLUSH.BASELINE`, `FLUSH.SCENARIO`: flushing baseline/scenario failure.

The message is suitable for the GUI; use `code` for categorization. A warning
alone does not imply failure. For the full schema see [diagnostics](doc/diagnostics.md).

## Python

[examples/integration/run_staci.py](examples/integration/run_staci.py) is a
standard-library adapter with isolated job directories, preserved input copies,
argument-array invocation, timeout handling and completion-record checks.

```sh
python3 examples/integration/run_staci.py \
  --binary /absolute/path/build/staci --network /absolute/path/network.inp \
  --job-root /absolute/path/jobs --eps --timeout 300
```

From Python, put the example directory on `sys.path` or copy the adapter into your
application:

```python
from run_staci import run_network
import json

job = run_network('/opt/staci/bin/staci', 'network.inp', 'jobs', eps=False)
if job['state'] == 'success':
    with open(job['result_file'], encoding='utf-8') as stream:
        hydraulics = json.load(stream)
else:
    for item in job['diagnostics']:
        if item.get('severity') in ('warning', 'error'):
            print(item.get('code'), item.get('message'))
```

The adapter may raise normal Python `OSError` exceptions for missing executables,
unreadable inputs or unwritable job directories. Handle these as service errors.
Its own command-line exit code is 0 only for success; the JSON `exit_code` retains
the solver's actual code, including 3 for partial output.

Read EPS CSV with `csv.DictReader` or pandas; HDF5 consumers can use h5py. Output
fields carry SI units: head/level m, flow m³/s, velocity m/s, age s, chemical
concentration kg/m³. **1 mg/L = 0.001 kg/m³.** Missing quality quantities may be
`nan` in CSV; consult metadata before displaying them.

## Web application / remote calls

Run STACI on the server through a worker process. A browser uploads an INP or
selects a server-side model; the backend creates a job, runs the adapter, then
returns job state and parsed results. The browser does not execute the native
binary. A useful API is `POST /jobs`, `GET /jobs/{id}`, and
`GET /jobs/{id}/results`, with logs retained per job.

Framework-neutral worker example:

```python
def calculate_job(validated_input_path, configured_binary, configured_job_root):
    return run_network(configured_binary, validated_input_path,
                       configured_job_root, eps=True, timeout=300)
```

Run this worker outside the HTTP event loop. Use your framework's background
queue or thread/process executor to wait for the subprocess. Store job status in
server-side storage; translate the returned `job_dir`/`result_file` into a job ID
and controlled downloads rather than exposing server paths. Limit simultaneous
workers according to available memory/CPU, set upload size/time limits, and keep
executable paths and CLI modes server-configured. Optimizer jobs need their own
XML/JSON settings and referenced model copies in the same job directory. They can
run substantially longer than a single hydraulic solve.

Cancellation belongs to the worker: terminate the owned solver process, retain
the log, mark the job cancelled, and do not invent a STACI `run_end` event.
No HTTP authentication, queue, database or web deployment is provided here.

## MATLAB

For the same process interface as the Python example, use
[examples/integration/run_staci_cli.m](examples/integration/run_staci_cli.m).
It requires MATLAB with Java enabled, but no Python or MEX build. It passes an
argument list directly to the executable, creates a separate directory for each
job, captures console output, enforces a timeout and checks the matching JSONL
completion record. Its return fields mirror the Python adapter; `diagnostics`
is a cell array of structs, and an unavailable `exit_code`/`result_file` is
`[]`/`''`. Filesystem and process-start failures raise MATLAB exceptions.

```matlab
addpath('/absolute/path/staci/examples/integration');
job = run_staci_cli('/absolute/path/build/staci', ...
    '/absolute/path/network.inp', '/absolute/path/jobs', ...
    'EPS', true, 'Timeout', 300);
if strcmp(job.state, 'success')
    metadata = jsondecode(fileread(job.result_file));
    nodes = readtable(fullfile(job.job_dir, 'result-nodes.csv'));
else
    fprintf('STACI: %s\n', job.state);
    disp(job.diagnostics);
end
```

Omit `'EPS', true` for a steady hydraulic solve; `result_file` then points to
the hydraulic JSON. [example_staci_cli.m](examples/integration/example_staci_cli.m)
demonstrates both modes, reading results and displaying errors/warnings. Partial
EPS output has state `partial`, not `success`; inspect frame convergence flags
before using it. Job files remain available after the function returns.

The CLI examples were run in MATLAB R2026a on macOS. The standalone regression
also checks invalid input, timeout handling and paths containing spaces:

```matlab
addpath('/absolute/path/staci/tests/matlab');
test_staci_cli_adapter('/absolute/path/build/staci');
```

The supported in-memory interface is documented in [matlab/README.md](matlab/README.md).
It requires MATLAB and a compiled MEX module; it does not spawn the CLI per solve.

```matlab
addpath('/absolute/path/staci/matlab');
build_staci_mex();
model = StaciModel('/absolute/path/network.inp');
cleanup = onCleanup(@() model.release());
try
    status = model.solveHydraulics();
    assert(status.converged, 'Hydraulic solve did not converge');
    nodes = model.nodeTable();
    links = model.linkTable();
    quality = model.solveSteadyQuality("both");
catch exception
    fprintf(2, '%s: %s\n', exception.identifier, exception.message);
    rethrow(exception);
end
```

MEX errors use `STACI:MatlabAPI`; they are MATLAB exceptions, not CLI process
exit codes. Check `status.converged` as well. The CLI's `run_start/run_end` JSONL
contract is a process interface and is not promised by direct MEX/C++ calls.
Each parallel MATLAB worker must construct its own model; do not share MEX
handles. `StaciModel` exposes steady hydraulics, property access, steady quality
and sensitivities. For **EPS from MATLAB**, use `run_staci_cli` above and read its
CSV/HDF5 output. Do not confuse `solveSteadyQuality` with EPS.

## C++ and extending the code

Link the CMake `staci_core` target and include `StaciSession.h` for in-memory
steady calculations. See `examples/cpp/` and the C++ example tests. `StaciSession`
owns its model; call `solve_hydraulics()` before quality/sensitivity operations.
IDs and public facade values use SI. This is a source-level C++ interface, not a
stable binary ABI; rebuild dependent applications with the same headers/compiler.
The legacy core and its diagnostics have shared state: isolate concurrent service
calculations in processes instead of assuming thread-safe sessions.

## Validation and known limits

The [package audit](doc/package_validation.md) records the validated build options,
CTest/portable/MATLAB results and installation checks. The optional MEX interface
was built and tested, in addition to both CLI adapters.

Provision the pinned independent EPANET reference with
`python3 tests/setup_epanet_reference.py`, configure `STACI_EPANET_LIBRARY` and
`STACI_EPANET_EXECUTABLE`, rebuild, then run CTest. Missing reference tools are
reported as skips, not numerical agreement. The public steady corpus and the
full-period EPS corpus establish different kinds of coverage.

`epanet_eps_chemical_reference` compares Net1, Net2 and Net2-CL2 over their configured
periods without changing inputs or loosening tolerances. It covers bulk/wall
reaction with mass transfer, complete-mix tank storage, an external concentration
source and segment transport. The chemical EPS implementation supports MIXED
tanks and first-order reactions. Other tank mixing laws and reaction orders are
explicit compatibility errors. The separate water-age implementation still has
limited tank-mixing support. Residual circulation in stagnant branches can still
affect reported chemical node concentrations; see the reference status for the
known PRV example and the unavailable stagnant Batch reference. Full-corpus EPS hydraulic/chemical equivalence is
not established; see [current reference status](doc/epanet_reference.md).

When extending an element, update import, unit conversion, hydraulic equation and
Jacobian, status transitions, export and diagnostics together. Add independent
EPANET comparisons for applicable regimes and malformed-input tests. Keep public
INP hashes intact; record intentional model adaptations in the manifest. Never
turn an unavailable/unbalanced reference or a missing comparison into a pass.
