# STACI MATLAB interface

`StaciModel` is a handle class backed by the general `staci_mex` module. The
constructor reads an SPR or INP network once. The C++ model then remains in
memory until `release`, object deletion, or MATLAB exit; repeated property
changes and calculations do not use intermediate files.

```matlab
addpath('matlab')
model = StaciModel('tests/anytown_1med.spr');
pipes = model.linkIds("pipe");
model.setLinkProperty(pipes, "diameter_m", ...
    model.getLinkProperty(pipes, "diameter_m") * 1.05);
status = model.solveHydraulics();
nodes = model.nodeTable();
links = model.linkTable();
```

All public numerical values use SI units. Set/get operations accept one id or
an id vector. Values may be scalar (broadcast to every selected element) or a
vector matching the ids.

Supported node properties:

| Property | Unit | Read | Write |
|---|---:|:---:|:---:|
| `elevation_m` | m | yes | no |
| `pressure_head_m` | m | yes | yes |
| `total_head_m` | m | yes | no |
| `demand_m3s` | m3/s | yes | yes |
| `water_age_s` | s | yes | yes |
| `concentration_kgm3` | kg/m3 | yes | yes |
| `source_concentration_kgm3` | kg/m3 | yes | yes |

Common link results are `flow_rate_m3s`, `mass_flow_rate_kgs`, `velocity_ms`,
`status`, `water_age_start_s`, and `water_age_end_s`. Pipe properties are
`diameter_m`, `length_m`, `friction_coeff`, `minor_loss`,
`bulk_reaction_per_s`, and `wall_reaction_m_s`. Type-specific writable
properties include pump `speed`/`power_w`, valve `position`/`tcv_setting`, pool
`water_level_m`/`bottom_level_m`, and fixed-pressure `boundary_head_m`.

Quality and sensitivity calls require a converged hydraulic state:

```matlab
quality = model.solveSteadyQuality("both");
sensitivity = model.hydraulicSensitivity("PIPE74", "diameter_m");
```

Sensitivity parameters currently supported are pipe `diameter_m`, pipe
`friction_coeff`, and nodal `demand_m3s`. Result derivatives are returned for
every link flow and node pressure head in consistent SI units.

Build explicitly with `build_staci_mex()` or allow the first `StaciModel`
constructor to build automatically. MATLAB, CMake, a C++17 compiler, Eigen3
and nlohmann/json are required. For a preconfigured build (including an explicit
macOS SDK when needed), use CMake with `STACI_BUILD_MATLAB_MEX=ON`, build target
`staci_mex`, and add the build's `matlab` directory to MATLAB's path. On Windows the CMake generator must be compatible with the
installed MATLAB MEX compiler.

Each MATLAB worker is a separate process. Parallel optimization workers must
therefore construct their own `StaciModel`; a MEX handle cannot be shared
between workers.

## Validation and external-process alternative

The in-memory interface is covered by `tests/matlab/run_matlab_tests.m` (six
tests). The current MEX module was built and these tests passed in MATLAB
R2026a on macOS arm64. Each build must match the MATLAB architecture.

For a process-based call, including EPS, use
[run_staci_cli.m](../examples/integration/run_staci_cli.m) and
[example_staci_cli.m](../examples/integration/example_staci_cli.m). This requires
MATLAB's JVM but no MEX or Python. The separate
`test_staci_cli_adapter('/absolute/path/staci')` checks successful steady/EPS
runs, invalid input, timeouts and paths containing spaces. See
[README_DEV](../README_DEV.md#matlab) for usage and diagnostics.
