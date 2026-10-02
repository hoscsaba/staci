# Package validation — 2026-10-02

Validation platform: macOS arm64, Release build, Xcode toolchain, MATLAB R2026a.
Optimizers, HDF5, C++ examples and the pinned independent EPANET reference were
enabled. The optional MEX module was built separately with the Eigen backend.
These results do not certify Windows/Linux builds or full-corpus EPS equivalence.

| Check | Result |
|---|---|
| Complete CTest suite, including optional C++ EPANET EPS example | **84 passed, 0 failed, 0 skipped** |
| Portable `tests/run_tests.py` with required independent reference | **46 passed, 0 failed** |
| Large native SPR solves omitted by the portable runner's default size limit | **2 passed**, exit 0 and `OK` markers |
| MATLAB MEX interface | **6 passed, 0 failed, 0 incomplete** |
| MATLAB CLI adapter | Steady/EPS success, invalid input, timeout and paths with spaces passed |
| Temporary installation | All four installed applications returned exit 0 for `--help` |
| Public input integrity | All **97** SHA256 values match the manifest |
| Maintained Markdown links and JSON examples | Local file links resolve; example JSON parses |
| Generated HTML/API documentation | Regenerated with Doxygen 1.17.0 from current sources and guides |

The extra C++ toolkit example raises the configured CTest count from 83 to 84.
Python integration, JSON auxiliary inputs, common application diagnostics,
optimization, hydraulic/quality models and independent comparisons are included
in CTest. Expected diagnostic failures are regression passes, not solved models.

The two additional large native solves were
`LOV-LOVOTV-2-input_mod2.spr` and `VIZ-SOPTVR-2-660-input.spr`. Tests ran on copies;
the corpus and original models were preserved.

## Fixes made during this audit

- Linked the optional `staci_matlab_core` target to
  `nlohmann_json::nlohmann_json`. A clean MEX build previously failed to find the
  JSON header even though ordinary application builds succeeded.
- Updated the detailed usage guide, Hungarian split manual, corpus README,
  Anytown status, historical error classification, MATLAB guide and testing
  documentation. Removed outdated current claims that ky10 failed, valves and
  emitters were missing, tank chemistry was absent, or settings required XML.
- Corrected the developer guide's table formatting and documented explicit SDK
  selection for macOS installations with mismatched Command Line Tools/Xcode.

Historical investigation results remain explicitly marked as historical.
The maintained usage/validation documents are Markdown; `doc/html` was regenerated
from the current sources and guides using `doxygen doxygen.config`. The Doxygen
input list now includes the developer guide, current documentation and examples,
and uses GitHub-compatible Markdown anchors.

## Reference scope and remaining limitations

Both default and strict snapshot profiles give 89 solved networks, four input
rejections and four physical/numerical failures. The strict independent reference
has **88 numerical matches**, three EPANET input rejections and six unreliable
references. No input or expected-error case is counted as numerical agreement.

Full-period discrepancies remain in Net3/Net6 variants and the stagnant PRV
chemical branch. Batch's stagnant chemical reference remains unavailable;
non-MIXED chemical tanks and non-first-order reaction laws are explicitly
unsupported. See [reference status](epanet_reference.md) for details. Existing
legacy `nrutil_nr.h` uninitialized-pointer compiler warnings remain; the build
is successful, but is not warning-free.

## Reproduction

Provision the pinned reference as described in [testing](testing.md), then:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DSTACI_EPANET_LIBRARY=/absolute/path/epanet/build/lib/libepanet2.dylib \
  -DSTACI_EPANET_EXECUTABLE=/absolute/path/epanet/build/bin/runepanet \
  -DSTACI_EPANET_TOOLKIT_ROOT=/absolute/path/epanet
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 tests/run_tests.py --binary build/staci \
  --calibrate-binary build/staci_calibrate --split-binary build/staci_split \
  --epanet-binary /absolute/path/epanet/build/bin/runepanet \
  --epanet-library /absolute/path/epanet/build/lib/libepanet2.dylib \
  --require-epanet-reference --full-spr-hydraulics
```

Use the platform's library suffix. The portable runner recreates its
`test-results` directory; use `--tests-dir` with an isolated copy to preserve
previous reports. The audit used isolated copies and ran the two large SPR
solves separately. For MATLAB, build `staci_mex`, add its output directory and
`tests/matlab` to the MATLAB path, then run `run_matlab_tests()` and
`test_staci_cli_adapter('/absolute/path/staci')`.
