"""Audit saved measurements without fitting parameters or rerunning circuits."""
import csv
from pathlib import Path
import numpy as np
from .logical_model import branches, multiplier
from .metrics import wilson


def audit_samples(output):
    from .experiment import write_csv
    output = Path(output)
    with (output/"short_state_results.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    hist_rows, conditional = [], []
    worst = 0.0
    for row in rows:
        theta, t = float(row["theta"]), int(row["cycles"])
        if row["kind"] != "encoded" or row["decoder"] != "standard" or row["state"] != "0":
            continue
        if float(row["readout_error"]) or t == 0 or theta == 0:
            continue
        data = np.load(output/"records"/f"{row['record_id']}.npz")
        syn, failure = data["syndromes"], data["outcome"]
        # Complete joint two-cycle law; single-cycle marginal for T=1/4.
        labels = syn[:, 0] + 8*syn[:, 1] if t == 2 else syn.ravel()
        q = branches(theta).q
        p = np.outer(q, q).ravel() if t == 2 else q
        counts = np.bincount(labels, minlength=len(p))
        for label, (n, expected) in enumerate(zip(counts, p)):
            ci = wilson(int(n), len(labels))
            hist_rows.append(dict(theta=theta, cycles=t, model=row["model"],
                distribution="joint_s0_plus_8s1" if t == 2 else "pooled_syndromes",
                label=label, counts=int(n), samples=len(labels),
                observed=float(n/len(labels)), expected=float(expected),
                ci_low=ci[0], ci_high=ci[1]))
        k = np.count_nonzero(syn, axis=1)
        for kk in np.unique(k):
            selected = k == kk
            n = int(selected.sum())
            success = int(np.count_nonzero(failure[selected] == 0))
            expected = float((1+multiplier(theta, t, kk, row["model"]))/2)
            ci = wilson(success, n)
            conditional.append(dict(theta=theta, cycles=t, model=row["model"],
                nonzero_count=int(kk), shots=n, successes=success, observed=success/n,
                ideal_probability=expected, ci_low=ci[0], ci_high=ci[1]))
            # Do not use Gaussian diagnostics in sparse or deterministic strata.
            if n*expected*(1-expected) >= 5:
                z = abs(success/n-expected)/np.sqrt(expected*(1-expected)/n)
                worst = max(worst, float(z))
    write_csv(output/"syndrome_histograms.csv", hist_rows)
    write_csv(output/"conditional_channel_validation.csv", conditional)
    with (output/"calibration_batches.csv").open(encoding="utf-8") as f:
        cal = list(csv.DictReader(f))
    summaries = []
    for theta in (0.1, 0.05):
        for model in ("F", "R"):
            selected = [r for r in cal if float(r["theta"]) == theta and r["true_model"] == model]
            n = sum(int(r["budget"]) for r in selected)
            count = sum(int(r["nonzero"]) for r in selected)
            expected = float(selected[0]["exact_nonzero_rate"])
            ci = wilson(count, n)
            z = abs(count/n-expected)/np.sqrt(expected*(1-expected)/n)
            summaries.append(dict(theta=theta, model=model, shots=n, nonzero=count,
                observed=count/n, exact_probability=expected, ci_low=ci[0], ci_high=ci[1],
                standardized_residual=float(z)))
    write_csv(output/"calibration_rates.csv", summaries)
    if worst > 6 or max(r["standardized_residual"] for r in summaries) > 6:
        raise AssertionError("Measured conditional channels/calibration violate prespecified 6-SE check")
    return dict(conditional_max_standardized_residual=worst,
                calibration_max_standardized_residual=max(r["standardized_residual"] for r in summaries),
                full_two_cycle_syndrome_law="saved for F and R; ideal predictions identical",
                calibration_rates="pooled independent nonempty batches; known angle")
