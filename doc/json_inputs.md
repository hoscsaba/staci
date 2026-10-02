# JSON auxiliary inputs

All runtime file inputs other than network definitions have JSON variants.
Network definitions remain native SPR/XML or EPANET INP. Settings embedded in
the network stay part of that network format. Existing XML, measurement CSV and
hydrant text files remain supported. JSON comments are not allowed.

| Input | Legacy format | JSON format / selection |
|---|---|---|
| `staci_split` settings | XML `<settings>` | Flat JSON object; `--settings PATH` |
| `staci_calibrate` settings | XML `<settings>` | Flat JSON object; `--settings PATH` |
| Calibration measurements | Semicolon-separated CSV | `measurements` array; path in `sollwert_dfile` |
| `staci -i` initial values | Native XML result fields | `nodes` / `edges` arrays |
| Flushing configuration | Already JSON | `--config PATH.json`; see [flushing](flushing.md) |
| Flushing hydrants | One node ID per text line | `hydrants` array; `--hydrants PATH.json` |

## Optimizer settings

```sh
staci_split --settings partition.json --seed 12345
staci_calibrate --settings calibration.xml --seed 12345
```

The extension selects XML or JSON, case-insensitively; there is no format flag
or content sniffing. Other settings extensions are rejected. Without `--settings`,
the historical `staci_split_settings.xml` / `staci_calibrate_settings.xml` is
used when present, otherwise its `.json` sibling. If both exist, XML wins;
explicit `--settings` always wins. An invalid selected XML never silently falls
back to JSON.

JSON uses the same case-sensitive keys as the XML child names, without a
`settings` wrapper. Numeric fields accept JSON numbers or numeric strings;
nonfinite numbers and fractional integer settings are rejected. Paths and
enumerations are strings. `Spoil_Active_Pipes` accepts `true` / `false` as well
as the legacy `"yes"` / `"no"`. Required fields and conditional pipe-selection
fields are the same as in XML.

Paths inside either format remain relative to the **process working directory**,
including output logs. Selecting settings in another directory does not change
this base. Calibration retains its `dir_name` + `fname_prefix` + period + `.spr`
network naming convention; `dir_name` needs its trailing separator. Measurement
paths are also formed as `dir_name` + `sollwert_dfile`.

Complete examples: [split](../examples/config/staci_split_settings.json),
[calibrate](../examples/config/staci_calibrate_settings.json). Replace their
network paths and IDs before running. These are configurations, not bundled
network definitions.

## Calibration measurements

```json
{"measurements": [
  {"id": "NODE64", "type": "node", "values": [12.3, 12.5]},
  {"id": "POOL1", "type": "pool", "values": [3.1, 3.0]}
]}
```

Values are numerical pressure heads (`node`) or water levels (`pool`), in metres,
one value per period. `Start_of_Periods` is a zero-based offset into `values`;
at least `Start_of_Periods + Num_of_Periods` values are required. JSON has no
trailing empty CSV column. Select it with `"sollwert_dfile": "targets.json"`.
XML settings may reference JSON measurements and JSON settings may reference CSV.

Selected CSV values must be finite numbers. Empty cells, nonnumeric text,
trailing junk, `nan`, `inf` and out-of-range values are rejected with exit code
2 and an `INPUT.CONFIG` diagnostic identifying the file, measurement row,
element ID and zero-based period. Surrounding whitespace and scientific notation
are accepted; invalid measurements are never silently replaced with zero.

For multi-period calibration, each candidate starts with the measured initial
pool levels. Each following period's pool level is calculated from that
candidate's preceding solved flow and `dt`, before solving the new period.
`Start_of_Periods` selects the input window; storage updates still begin with
the second period within that window.

## Hydraulic initial values

```sh
staci -s network.spr -i initial_values.json
```

```json
{
  "nodes": [{"id": "NODE64", "pressure_pa": 98100,
             "concentration_kg_m3": 0.0005}],
  "edges": [{"id": "PIPE1", "mass_flow_rate_kg_s": 1.2}]
}
```

Both arrays may be partial; omitted elements retain their loaded initial values.
At least one array must be present. Concentration is optional and nonnegative.
Node pressure is in Pa, converted using the legacy XML convention of
1000 kg/m³ and 9.81 m/s²; edge mass flow is in kg/s. Unknown IDs, duplicate IDs,
unknown fields, missing values and invalid numeric types are input errors. This
JSON is an initialization file, not a JSON network or the hydraulic output schema.
The `-i` option now retains its requested input through solver setup, for both
legacy XML and JSON.

## Flushing

The existing [configuration example](../examples/flushing/flushing_config.json)
is already JSON. A text hydrant list has this equivalent JSON form:

```json
{"hydrants": [{"node_id": "J1", "status": "matched"}]}
```

An optional `dxf_handle` identifies the hydrant asset; otherwise its node ID is
used. The configuration's `hydrant_node_ids` array is another supported way to
specify hydrants. See [flushing documentation](flushing.md) for precedence and
network requirements.

## Errors and verification

New auxiliary input validation uses the shared diagnostic log, exit code 2 and
`INPUT.CONFIG`, with the file and offending field/ID. Malformed legacy XML keeps
its existing XML diagnostics. Flushing retains its existing diagnostic contract.
JSON support does not alter the hydraulic algorithms or optimizer objectives.

`json_auxiliary_inputs` compares seeded XML/JSON optimizer results, tests JSON-only
default discovery, mixed-case extensions, measurement loading, initial values
and malformed-input diagnostics. Existing tests continue to cover legacy XML/CSV
and flushing JSON inputs.
