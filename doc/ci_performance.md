# GitHub build performance

## Measured baseline

Successful run [39258308194](https://github.com/hoscsaba/staci/actions/runs/39258308194)
(commit `923943e`, 2026-10-19), before the cache/shared-core changes:

| Platform | Install dependencies | Build | CTest |
| --- | ---: | ---: | ---: |
| Linux x64 | 31 s | 2 min 41 s | 42 s |
| macOS | 1 min 5 s | 53 min 53 s | 1 min 2 s |
| Windows x64 | 52 min 30 s | 8 min 6 s | 1 min 6 s |

Windows rebuilt the vcpkg dependency tree (including Boost, BLAS/LAPACK,
SuiteSparse and pagmo) on each fresh runner. macOS spent its time compiling,
not installing Homebrew dependencies. The previous CMake configuration compiled
the hydraulic sources separately for six normal targets; the workflow also used
`--parallel` with no limit. The logs show many simultaneous compiler starts.
Memory contention is a plausible contributor to the highly variable macOS times;
the old logs do not contain memory/swap measurements to establish it directly.

## Changes

- The CLI, optimizers, flushing tool and C++ examples reuse `staci_core`.
  MATLAB keeps a separate Eigen/PIC solver variant because its compilation
  requirements differ. Solver functionality and the test suite remain enabled.
- Build parallelism is bounded: 2 jobs on macOS, 4 on Linux/Windows.
- Linux/macOS use CMake's `ccache` compiler launcher. The cache is limited to
  500 MB per platform and compiler content is checked. Cache statistics are
  printed after the build. Homebrew dependencies are still installed normally:
  their measured install time is small compared with compilation.
- Windows uses vcpkg's `files` binary cache, restored/saved with
  `actions/cache`. Installation always runs so vcpkg checks ABI hashes and
  rebuilds incompatible/missing packages. The runner's vcpkg revision is part
  of the preferred cache prefix; fallback archives can supply ABI-compatible
  packages after a revision change. No installed tree or CMake build directory
  is reused blindly.
- Binary caches are saved directly after dependency installation; compiler caches
  after a successful build. A later application-test failure does not discard
  an expensive successful dependency build.
- Each run/attempt uses a new cache key with platform-specific restore prefixes,
  allowing updated entries to be saved despite immutable Actions cache keys.
  Old caches remain subject to GitHub's normal quota/eviction policy.
- A newer push cancels an obsolete run on the same ref. Different refs remain
  independent. CTest remains serial to preserve existing file-isolation behavior.

The first Windows run after enabling caching still needs a cold dependency
build. A warm run is needed to measure the benefit; a restored archive alone is
not proof of a vcpkg/ccache hit. Runner image/compiler/package updates can require
rebuilding some packages. No fixed speedup is promised before measured runs.

## How to verify or troubleshoot

1. Compare step durations for the same commit using an initial run and a
   `workflow_dispatch` run on `master` after the first run has saved its caches.
2. Check `Restore Windows binary cache`, then vcpkg's restore/build counts.
   For Linux/macOS check `Report compiler cache statistics` (hits and misses).
3. Confirm **all three** `Test hydraulics and flushing` steps, executable smoke
   tests and bundle uploads succeed. Timing improvements must not reduce coverage.
4. Hosted CI now builds the pinned EPANET reference and requires all reference
   comparisons to execute. `check_ci_test_results.py` rejects missing required
   tests, skips or failures in the JUnit report. Local builds can still omit the
   optional toolkit; see [package validation](package_validation.md).
5. To force a cold cache, change the `v1` cache namespace in the workflow or
   delete the relevant cache through GitHub Actions. Do not disable vcpkg ABI
   checks or ccache compiler/header checks to manufacture hits.

Configuration: [build-binaries.yml](../.github/workflows/build-binaries.yml).
Upstream references: [Actions cache](https://github.com/actions/cache),
[vcpkg binary caching](https://learn.microsoft.com/en-us/vcpkg/reference/binarycaching),
[ccache manual](https://ccache.dev/manual/latest.html),
[CMake parallel builds](https://cmake.org/cmake/help/latest/manual/cmake.1.html#build-a-project).

## EPANET reference in hosted CI

All three platforms build official OWA EPANET 2.2 at commit
`4d8d82ddc260fce216af9321fc3d9a4646ac6827` before configuring STACI.
`tests/setup_epanet_reference.py --cmake-cache FILE` loads the shared library,
checks `ENgetversion == 20200`, and writes explicit CLI/shared-library paths
for CMake. Both CLI report comparisons and Python C API comparisons are enabled.

The source and Release build are cached under `.epanet-reference`, keyed by
OS, architecture, hosted-runner image version and setup-script hash (which
contains the pinned revision). No cross-image fallback is used. On a hit,
CMake still reconfigures/builds and the library version is checked. On a miss,
the archive is downloaded and built. The cache is saved before STACI tests.

EPANET's reference build uses CMake and the runner's native compiler; its own
Boost-dependent unit tests and Python wrapper are disabled. It adds no external
package dependencies. STACI's existing vcpkg and compiler caches remain active.

The 18 previously skipped tests now run, and supplying the shared library also
registers two EPS reference tests: 85 CTest cases in the hosted configuration.
The optional C++ EPS demonstration remains outside this CI configuration (a
local build with `STACI_EPANET_TOOLKIT_ROOT` can enable that 86th test). The
reference is a test tool, not a new runtime dependency of the shipped programs.

CTest writes `build/ctest-results.xml`. The post-test gate requires the 20
reference cases by name and rejects any skipped/failed test. This prevents an
unavailable reference or accidentally disabled registration from turning green.

### Pinned EPANET parser safety correction

The setup script applies one bounded-loop correction to EPANET 2.2 `hour()` in
`src/input2.c`: `n <= 3` becomes `n < 3` for the three-element `double y[3]`.
The unchanged upstream `bad_times.inp` fixture contains a four-component time;
AddressSanitizer confirms a stack-buffer overflow in the unpatched source, and
Linux's stack protector aborts the public-corpus test. The correction preserves
the existing first-three-field interpretation and changes no hydraulic equations
or valid time values. This is a pinned EPANET 2.2 reference with this explicit
parser safety patch, not a byte-identical upstream build. The setup-script hash
invalidates the CI cache when this correction changes. The public-corpus test
continues to check the original fixture and its hash without weakening tolerances.
