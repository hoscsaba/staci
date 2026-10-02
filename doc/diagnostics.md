# Common application diagnostics

`staci`, `staci_split`, `staci_calibrate` and `staci_flush` append errors and warnings to the same JSON Lines file. The default is `staci-diagnostics.jsonl` in the working directory. A GUI should supply an absolute path using `--diagnostics-file PATH` for every invocation. Alternatively set `STACI_DIAGNOSTICS_FILE`; the command-line option takes precedence.

Each line is one JSON object. Common fields are `schema_version` (1), `timestamp` (UTC), `program`, `run_id`, `sequence`, `event`, `severity`, `code` and `message`. A run starts with `run_start` and finishes with `run_end`; the latter includes `exit_code`, `error_count` and `warning_count`. Diagnostic records contain human-readable messages and, where available, network, section, element and line fields extracted from the importer diagnostic. Consumers should tolerate additional fields.

| Exit code | Meaning |
| --- | --- |
| 0 | Successful operation, possibly with warnings |
| 1 | Calculation or execution failure |
| 2 | Invalid input, configuration or unsupported feature |
| 3 | Partial result: some frames or scenarios failed |

The GUI should group records by `run_id`, display warning/error `message` values, and use the process exit code and matching `run_end` to determine completion. An interrupted or killed process may have no `run_end`. Append writes are locked across processes; records from simultaneous invocations may interleave, but each record stays intact. A failure to open the diagnostics file is reported on stderr with `DIAGNOSTICS_OPEN` and exits with code 1.

Existing console output and network-specific result files remain available. Errors are also visible on stderr. Legacy fatal paths now pass through the common exception boundary, so they cannot silently exit successfully. Legacy messages use fallback codes such as `INPUT_OR_CONFIGURATION`; newer paths have specific codes including `INPUT.XML`, `CLI.INPUT`, `FLUSH.BASELINE`, `FLUSH.SCENARIO` and `EPANET.EPS_PARTIAL_FAILURE`. Numeric failures of discarded calibration candidates are warnings; a failed initial calibration baseline is an error.

The integration test `application_diagnostics` checks all four programs, invalid inputs, partial results, environment overrides, Unicode paths, log-open failures and simultaneous writers. POSIX locking is exercised on macOS; the Windows implementation has not been run in this validation.

Successful `staci -s network.inp` runs also write `network.inp.hydraulics.json` with convergence and SI node/link results. A new steady solve removes any previous hydraulic export before validating the input, preventing GUI reuse of stale results after rejection. See [epanet_reference.md](epanet_reference.md) for numerical validation.

## Hydraulic solver options

All four applications accept `--head-tolerance-m VALUE`, `--mass-tolerance-kg-s VALUE` and `--max-iterations INTEGER`, including `--option=value` syntax. Values must be positive and finite; iteration counts must be integers. Invalid values produce `CLI_ARGUMENT` and exit code 2. CLI values override input settings and application solver settings.

For INP files the default head residual limit is **0.0001 m (0.1 mm)** and the continuity residual limit is **1e-8 kg/s**. The head limit applies to the RMS equation residual, not the largest individual residual or the difference from an EPANET result. Existing explicit SPR settings remain in effect unless overridden. The iteration limit remains the input setting.

```sh
staci --epanet-eps ky16.inp -o ky16 --head-tolerance-m 0.0001 --mass-tolerance-kg-s 1e-8 --max-iterations 100
```

The `solver_options` integration test checks help and invalid values in all four applications, and exercises default, strict and iteration-limited ky16 simulations.
