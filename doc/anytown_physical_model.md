# Physical Anytown test model

The actual corpus file `tests/public_networks/wntr/wntr/tests/networks_for_testing/Anytown.inp` has been changed at the user's explicit request. It is now an adapted physical model rather than the upstream all-pumps-stopped design template. The 96 other retained files remain byte-identical to their upstream sources. The manifest records the adapted SHA256 and the original upstream SHA256 separately, with an explicit list of modifications. The original input is archived in `doc/model_sources/Anytown_upstream.inp` (also copied into `tests/test-results/anytown-numerics/original-input.inp`) and remains available at the pinned upstream URL.

## Model choices

| Part | Original | Revised |
| --- | --- | --- |
| Pumps 78, 79, 80 | Zero multipliers throughout speed patterns 2, 3, 4 | Explicit SPEED 1.0 and all 24 multipliers of each pattern set to 1.0 |
| Candidate pipes 110, 113, 114, 115, 116, 125 | 0.0001-inch diameters | 16 inches (406.4 mm), keeping all six connections |
| Inlet trunks 1, 2, 3 | 12, 12, 16 inches | 24 inches (609.6 mm) |
| Other distribution pipes | 8, 10, 12, 16 inches | 12, 16, 18, 24 inches respectively |
| Existing pipe 4 | 30 inches | 30 inches |

Pump speed is a dimensionless ratio of nominal rotational speed, not a velocity in m/s. A setting of 1.0 means 100% nominal speed; it does not change the supplied pump curve.

The network topology, source elevations, demands and demand pattern, pipe lengths/roughness, pump curves, tank dimensions/levels/limits and all timestep settings are retained. Changing just the tiny pipes and starting the pumps gave a numerically converged network with negative pressures. The inlet/distribution upgrades provide a usable pressure envelope at the specified demands, including the original daily peak. This is a hydraulic test model, not a claim that a real system has been engineered or cost-optimized.

## Validation

* Default 0.1 mm and strict 1e-12 m head residual profiles both converge and numerically match the independent pinned EPANET 2.2 initial-state reference. EPANET reports no warnings for that state.
* Initial junction pressure heads are approximately 31.31–76.71 m; maximum pipe velocity is approximately 1.83 m/s.
* A full 24-hour STACI run at 1e-9 m residual tolerance completes all 1441 minute-by-minute hydraulic states and 25 hourly frames, with no failed states. Recorded junction pressures stay approximately 31.31–86.00 m, recorded pipe velocities remain below 3 m/s, and tank levels remain within their original minimum/maximum limits. Initial tanks start at minimum level and fill from the operating pumps. Maximum-level clamp warnings are retained.
* The new standard `epanet_physical_models` CTest verifies these checks and the zero-speed cases below. Independent numerical agreement for this model is asserted at t=0; the 24-hour checks establish STACI execution and physical bounds, not full EPANET/STACI EPS equivalence.

## Zero pump speed

Both HEAD and POWER pump equations support zero speed: the pump enforces zero flow, its Jacobian avoids division by zero, and its state is reported as stopped/closed. It does not act as an open bypass pipe. A usable alternative supply lets the rest of the network converge. If every supply is unavailable, a network-level failure is expected even though the stopped pump itself is handled correctly; `EPANET.NO_AVAILABLE_SUPPLY` explains that situation.

Twelve independent EPANET comparisons cover both pump types: explicit SPEED 0, zero speed-pattern multiplier, numeric STATUS setting 0, and EPS shutdown/restart/shutdown at 0/1/2 hours. The rest of the fixture is supplied by a second reservoir, so these are valid hydraulic cases. All pass, with zero STACI pump flow when stopped and positive flow when restarted.

The POWER-pump steady JSON export now uses its actual operating status, fixing a previously misleading enabled flag at zero speed. EPS RULES status premises likewise use the pump operating status. EPANET's input SPEED 0 can retain a raw OPEN flag despite no pumping (a stopped HEAD pump can also retain a tiny penalty leakage and an insufficient-head warning). Reference extraction preserves `raw_status` and `pump_speed`, and normalizes `enabled` to require a positive speed. This normalization is common to all pumps and does not change their reference flows or pressures.

## Current corpus

The corpus has 97 retained models: 96 upstream originals plus this adapted
Anytown. Current strict initial-state validation has 88 numerical matches,
three EPANET input rejections and six unreliable references. `ky10` now matches.
See [current reference status](epanet_reference.md) for execution coverage,
full-period limitations and reproduction commands. Historical counts from the
original adaptation are superseded by that report.

The original-Anytown investigation and its 170 unsuccessful numerical-setting
experiments are retained in [anytown_convergence.md](anytown_convergence.md).
