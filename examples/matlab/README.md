# STACI MATLAB examples

These examples use the general `StaciModel` interface. A network is read once;
all property changes, hydraulic calculations and result transfers then happen
in memory without temporary SPR, INP or result files. All public numerical
values use SI units.

Start MATLAB in the repository root and run:

```matlab
addpath('matlab')
addpath('examples/matlab')

getting_started_hydraulics();
modify_network_in_memory();
steady_water_quality_demo();
hydraulic_sensitivity_demo();
```

The first call automatically builds `staci_mex` when needed. The examples are:

- `getting_started_hydraulics.m`: load a model, solve hydraulics, retrieve
  MATLAB tables and plot pressures and flows;
- `modify_network_in_memory.m`: read and modify all pipe diameters in one
  vectorized call, then re-solve without writing a network file;
- `steady_water_quality_demo.m`: calculate asymptotic steady water age from a
  converged hydraulic state;
- `hydraulic_sensitivity_demo.m`: obtain implicit link-flow and nodal-pressure
  sensitivities with respect to a selected pipe diameter;
- `optimize_anytown_residence_time.m`: Global Optimization Toolbox `ga`
  example that minimizes demand-weighted water age subject to pressure limits.

Run the optimization separately because its default configuration evaluates a
larger, 50-generation population:

```matlab
result = optimize_anytown_residence_time();
```

It writes its final reports and plot to `matlab/matlab-results/`. Objective
evaluations themselves remain file-free.
