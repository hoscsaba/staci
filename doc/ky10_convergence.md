# ky10: reference convergence investigation

The original `tests/public_networks/wntr/examples/networks/ky10.inp` is unchanged. STACI converges on its initial steady snapshot; the failing public-network test is a numerical disagreement with the default EPANET 2.2 reference, not a STACI convergence failure.

## Cause

The constant-power pump `~@Pump-11` (20 HP) feeds the pressure-reducing valve `~@RV-4` (139.99 PSI). STACI finds pump flow 0.0115681413253 m³/s, pump head gain 131.5222733 m and an active PRV with head loss 85.1305793 m.

The default EPANET reference (relative accuracy 1e-6, 1000 trials) returns effectively zero pump flow, a closed PRV and pump head gain 7.6095788 m, without an instability warning. Its maximum hydraulic head-equation error, queried through `EN_getstatistic(EN_MAXHEADERROR)`, is **7.6095788228 m**. This exceeds a 0.1 mm energy criterion by about 76,096 times.

EPANET caps the constant-power pump gradient at very low flow. Its default convergence test accepts relative flow change without requiring an absolute head-error limit. The closed PRV and near-zero pump flow therefore leave this initialization at an unreliable hydraulic state. Junction continuity alone does not detect this error.

## Independent checks

1. Temporarily solve the same EPANET project with RV-4 open, then restore its original 139.99 PSI setting and solve again. The final model retains the original valve setting. All **7021** comparisons pass against STACI. Maximum node-head difference is **8.827737872e-6 m (0.00883 mm)**; maximum flow difference is **9.903334147e-8 m³/s**.
2. Without valve preconditioning, require absolute `HEADERROR = 0.0001 m` (converted to feet for this US-unit input) and use `DAMPLIMIT = 0.01`. EPANET reaches the active-PRV solution, with pump flow 0.011568143210 m³/s and head gain 131.522251257 m. Requiring head error alone does not converge within 1000 trials; damping is necessary in this experiment.
3. Query the default reference head-error statistic across the previously comparable corpus. ky10 is the clear outlier; the next largest successfully measured head error is below 6e-7 m.

These are initial-snapshot checks, not extended-period equivalence checks. The 0.00883 mm comparison is from the warm-start experiment; the numerical-settings experiment separately confirms the same pump/valve operating state.

## Consequence for the test harness

No new hydraulic element or ky10 input modification is indicated. The reference harness should validate absolute energy convergence and retry unreliable references with suitable numerical settings. Applying damping globally was also tested: it introduces small-flow comparison failures in six otherwise matching cases. A future harness change should preserve the successful initial method and apply a checked fallback, rather than globally changing settings or relaxing comparison tolerances.

This investigation does not change the production solver or reference harness. Saved experiment summaries are in `tests/test-results/ky10-investigation/`.
