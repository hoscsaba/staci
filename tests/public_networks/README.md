# Public EPANET hydraulic regression corpus

This directory contains 97 INP models (95 distinct byte sequences): 96 unmodified upstream files and one intentionally adapted physical Anytown:
42 from USEPA/WNTR and 55 from OpenWaterAnalytics/epanet-example-networks.
They range from a single pipe to 57,460 pipes and 28,899 junctions. The corpus
includes Net1/2/3/6, Kentucky 4/10, Anytown, Hanoi, Exeter, large NW models,
and targeted zero-flow, check-valve, pump, tank, valve, control, emitter and
pressure-dependent-demand cases. Upstream intentionally invalid input fixtures
are included as diagnostic cases.

Sources:

- https://github.com/USEPA/WNTR — revision
  `11eec7d5bf0a9367ffe04d75707918f0c0a890a1`
- https://github.com/OpenWaterAnalytics/epanet-example-networks — revision
  `4f0ceb2641c35adece59145881f10232c40c26b7`

`manifest.json` records the exact upstream path, revision, raw download URL,
SHA256, element counts, valve types, expected outcome and diagnostic checks. Anytown records both its upstream and adapted SHA256, and its local modifications.
Upstream WNTR license and OWA README are preserved beside this file. Network
files retain their original comments, encodings and line endings except the explicitly documented Anytown adaptation. No download
or Python package installation is needed to run tests.

## Scope and expected outcomes

Each case calls **only the ordinary `staci -s` executable** on an isolated copy.
It checks the exit code, `.rrs` convergence marker and that the INP remains
byte-identical. It does not run `staci_split`, `staci_flush` or `staci_calibrate`.

This is an initial hydraulic snapshot: tanks use their initial levels, demand
and source patterns use their first multipliers, dynamic controls/rules are
not executed, and chemical/MSX transport is not tested. The importer prints
these limitations; the per-network report preserves all compatibility
warnings. A converged snapshot does not certify EPS or quality compatibility
or agreement with a numerical EPANET reference. Existing targeted EPANET
reference-comparison tests remain separate.

Expected outcomes are explicit rather than accepting any failed run:

- `solved`: exit 0 and an `OK` convergence marker.
- `compatibility`: exit 2, no convergence marker, and a diagnostic naming the
  network, section/element/line and specific unsupported physics or input
  problem. Malformed inputs, unknown
  sections and invalid hydraulic records are rejected before solving a
  silently altered network. Import/list/export modes still preserve metadata.
- `nonconvergence`: a known physical/numerical failure, exit 1, an `ERROR!` marker,
  RMS residuals with units and limits, and worst node/link IDs. A future
  converged result is accepted as an improvement. These cases are **not**
  counted as successful hydraulic solves.

Timeouts, crashes, unexpected errors, changed inputs and degraded diagnostics
always fail the regression suite. The intentionally malformed time record in
WNTR `bad_times.inp` is unused by a snapshot; its successful snapshot does not
assert that its EPS time configuration is valid.

## Run

Included in the normal CTest suite and in `tests/run_tests.py`:

```sh
ctest --test-dir build --output-on-failure
ctest --test-dir build -R 'public_epanet_networks|epanet_hydraulic_diagnostics' --output-on-failure
```

For an independent corpus run:

```sh
python3 tests/test_public_networks.py --binary build/staci \
  --output-dir build/public-network-tests
```

On Windows use `build/Release/staci.exe` and add `-C Release` to CTest.
`--network wntr/examples/networks/Net3.inp` selects one case. `--jobs` defaults
to 4 and `--timeout` to 120 seconds per network. CTest's overall timeout is
1800 seconds. Results include `report.txt`, `report.json`, per-case console
logs, solver logs and convergence markers in the output directory.

`tests/test_epanet_diagnostics.py` additionally checks malformed numbers,
invalid dimensions, duplicate IDs, missing endpoints/curves, source-free
components, unknown units/sections, unsupported head-loss/demand models,
malformed emitters, PDA settings and valve curves. A positive case verifies that
EPANET node and link IDs may share a name in their separate namespaces.

## Current numerical verification (2026-10-02)

Strict initial-state validation has **88 numerical matches**, three EPANET
input rejections and six unreliable references across the 97 retained files.
Standalone strict execution has 89 solved cases, four input/compatibility
rejections and four physical/numerical failures. `ky10` and `JEP5-13` now match.
All EPANET valve types, emitters and PDA are implemented; this does not establish
full-period or water-quality equivalence.

```sh
python3 tests/setup_epanet_reference.py
python3 tests/test_public_networks.py --binary build/staci \
  --reference-library /absolute/path/libepanet2.so --require-reference \
  --output-dir build/public-network-reference-tests
```

Use `libepanet2.dylib` on macOS or `epanet2.dll` on Windows. Configure
`STACI_EPANET_LIBRARY` in CMake for the standard numerical corpus test.
The independent unit tests cover all ten flow units and legacy formats, and
verify that perturbed outputs are rejected. See [the reference protocol](../../doc/epanet_reference.md)
for tolerances, scope, convergence handling and full-period limitations.
Expected diagnostic cases and unavailable references are never numerical matches.
Use `--head-tolerance-m 0.0001` to compare the default profile; the reference
harness default is 1e-12 m.

`wolf-3` was removed after its independent reference failed; provenance remains
in `excluded_networks.json`. Reference acceptance includes independent junction
conservation and physical pump checks. The actual `Anytown.inp` uses nominal
pump speed and physically sized pipes; see [physical-model validation](../../doc/anytown_physical_model.md).
The original input investigation is retained explicitly as historical analysis.
