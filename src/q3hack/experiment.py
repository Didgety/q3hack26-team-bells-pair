"""One end-to-end experiment shared by the CLI and notebook."""
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from scipy.stats import binom
from .calibration import fit_model, nonzero_probability
from .circuits import (run_storage, run_calibration, run_bare, storage_circuit,
                       calibration_circuit, resource_counts)
from .code import STATES
from .decoders import Decoder
from .logical_model import (branches, physical_branches, exact_fidelity,
                            bare_fidelity)
from .metrics import (mean_fidelity, wilson, sample_logical_records, score,
                      paired_hierarchical_bootstrap)

NAMES = ("standard", "independence", "calibrated", "oracle")


@dataclass(frozen=True)
class Config:
    """Run sizes, not physical device settings.

    short_shots is per short-circuit condition; budgets are calibration shots
    per independent batch. Each batch has evaluation_shots separate long-model
    records, reused across decoders/budgets for paired comparisons.
    bootstrap_repeats sets Monte Carlo precision of the reported interval;
    seed is the root of disjoint stage-specific reproducibility streams.
    Angles and storage lengths are prespecified below, in radians and intervals.
    """
    seed: int = 20261014
    short_shots: int = 2048
    calibration_repeats: int = 24
    evaluation_shots: int = 4096
    bootstrap_repeats: int = 1500
    budgets: tuple = (0, 32, 128, 512, 2048)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, obj):
    def default(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj, indent=2, default=default, allow_nan=False), encoding="utf-8")


def _decoders(theta, model, fit):
    return (Decoder(theta), Decoder(theta, "R"), Decoder.calibrated(fit), Decoder(theta, model))


def _six_state_summary(rows):
    grouped = {}
    for row in rows:
        key = tuple(row[k] for k in ("theta", "cycles", "model", "readout_error", "kind", "decoder"))
        grouped.setdefault(key, []).append(row)
    out = []
    for key, group in grouped.items():
        assert len(group) == 6
        ci = np.array([wilson(r["successes"], r["shots"], confidence=1-0.05/6) for r in group]).mean(axis=0)
        value = np.mean([r["state_fidelity"] for r in group])
        exact = np.mean([r["ideal_state_fidelity"] for r in group])
        out.append(dict(zip(("theta", "cycles", "model", "readout_error", "kind", "decoder"), key),
                        f_mean=float(value), ci_low=float(ci[0]), ci_high=float(ci[1]),
                        fe_inferred=float((3*value-1)/2), ideal_f_mean=float(exact),
                        shots=sum(r["shots"] for r in group),
                        interval="conservative six-state Bonferroni-Wilson 95%"))
    return out


def exact_selection_probability(theta, model, shots):
    ks = np.arange(shots+1)
    choose_f = np.array([fit_model(int(k), shots, theta).model == "F" for k in ks])
    return float(np.dot(binom.pmf(ks, shots, nonzero_probability(theta, model)), choose_f))


def run(output, config=Config(), progress=print):
    """Execute short quantum circuits and calibration, then long logical sampling.

    Write CSV/JSON summaries, compressed raw records, figures, and REPORT.md
    under output; return the manifest dict. Main and calibration Aer runs
    have engineered Rx noise but no NoiseModel. Only short report-error
    sensitivity adds ancilla ReadoutError. Long records are classical samples
    of validated conditional channels, never long Qiskit/QPU executions.
    """
    if sys.prefix == sys.base_prefix:
        raise RuntimeError("Execute using the project .venv")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "records"
    raw.mkdir(exist_ok=True)
    start = perf_counter()
    # Disjoint namespaces: short evaluation=11, short calibration=12,
    # long evaluation=21, temporal calibration=22, bootstrap=31.
    used_seeds = {}

    def seed(domain, *indices):
        value = int(np.random.SeedSequence([config.seed, domain, *indices]).generate_state(1)[0])
        key = "/".join(map(str, (domain, *indices)))
        if value in used_seeds.values() and key not in used_seeds:
            raise RuntimeError("Unexpected seed collision")
        used_seeds[key] = value
        return value

    accounting = dict(qiskit_storage_shots=0, qiskit_calibration_shots=0,
                      qiskit_bare_shots=0, qiskit_circuits=0,
                      qiskit_execution_and_transpile_seconds=0.0,
                      logical_evaluation_records=0, logical_sampling_seconds=0.0,
                      decoder_scoring_seconds=0.0, bootstrap_seconds=0.0)
    checks = []
    for theta in (0, 0.05, 0.1, 0.2):
        for xi in (-1, 1):
            actual, leakage = physical_branches(theta, xi)
            residual = float(np.max(np.abs(actual-branches(theta, xi).kraus)))
            if residual > 1e-12 or leakage > 1e-12:
                raise AssertionError("Physical branches disagree; refusing long-model run")
            checks.append(dict(theta=theta, sign=xi, max_operator_residual=residual,
                               max_leakage_norm=float(leakage)))
    write_json(output / "branch_validation.json", checks)
    progress("Qiskit branch validation passed; executing short repeated-QEC circuits.", flush=True)

    short_rows, short_cal = [], []
    case = 0
    # Ordinary short validation plus a limited ancilla-readout sensitivity.
    settings = [(theta, t, model, 0.0)
                for theta in (0.1, 0.2) for t in (0, 1, 2, 4) for model in ("F", "R")]
    settings += [(0.0, 4, "F", 0.0)]
    settings += [(0.1, 4, model, p) for model in ("F", "R") for p in (0.01, 0.03)]
    fits = {}
    for ti, theta in enumerate((0, 0.1, 0.2)):
        for mi, model in enumerate(("F", "R")):
            syn, seconds, count = run_calibration(theta, model, 512, seed(12, ti, mi))
            fit = fit_model(int(np.count_nonzero(syn)), len(syn), theta)
            fits[theta, model] = fit
            short_cal.append(dict(true_model=model, **asdict(fit)))
            np.savez_compressed(raw / f"short_cal_theta{theta}_{model}.npz", syndromes=syn)
            accounting["qiskit_calibration_shots"] += len(syn)
            accounting["qiskit_circuits"] += count
            accounting["qiskit_execution_and_transpile_seconds"] += seconds
    write_json(output / "short_calibration.json", short_cal)
    for theta, t, model, readout in settings:
        fit = fits[theta, model]
        for si, state in enumerate(STATES):
            rec = run_storage(theta, t, model, state, config.short_shots, seed(11, case, si), readout)
            rec_id = f"short_{case}_{si}"
            np.savez_compressed(raw / f"{rec_id}.npz", syndromes=rec.syndromes, outcome=rec.outcome)
            accounting["qiskit_storage_shots"] += config.short_shots
            accounting["qiskit_circuits"] += rec.circuit_count
            accounting["qiskit_execution_and_transpile_seconds"] += rec.seconds
            for name, decoder in zip(NAMES, _decoders(theta, model, fit)):
                scores = rec.score(decoder.actions(rec.decoder_input()))
                successes = int(scores.sum())
                fe = exact_fidelity(theta, t, model, decoder.assumed_model)
                exact = 1.0 if state in ("+", "-") else fe
                ci = wilson(successes, len(scores))
                short_rows.append(dict(theta=theta, cycles=t, model=model,
                    readout_error=readout, kind="encoded", state=state, decoder=name,
                    state_fidelity=float(scores.mean()), successes=successes,
                    shots=len(scores), ci_low=ci[0], ci_high=ci[1],
                    ideal_state_fidelity=exact, record_id=rec_id))
            if readout == 0:
                bare = run_bare(theta, t, model, state, config.short_shots, seed(13, case, si))
                np.savez_compressed(raw / f"bare_{case}_{si}.npz", outcome=bare.outcome)
                accounting["qiskit_bare_shots"] += config.short_shots
                accounting["qiskit_circuits"] += bare.circuit_count
                accounting["qiskit_execution_and_transpile_seconds"] += bare.seconds
                successes = int(np.count_nonzero(bare.outcome == 0))
                ci = wilson(successes, config.short_shots)
                short_rows.append(dict(theta=theta, cycles=t, model=model,
                    readout_error=readout, kind="bare", state=state, decoder="none",
                    state_fidelity=successes/config.short_shots, successes=successes,
                    shots=config.short_shots, ci_low=ci[0], ci_high=ci[1],
                    ideal_state_fidelity=1.0 if state in ("+", "-") else bare_fidelity(theta, t, model),
                    record_id=f"bare_{case}_{si}"))
        case += 1
        if case % 8 == 0:
            progress(f"Short-circuit settings completed: {case}/{len(settings)}", flush=True)
    write_csv(output / "short_state_results.csv", short_rows)
    six = _six_state_summary(short_rows)
    write_csv(output / "short_six_state_results.csv", six)
    # A predetermined loose 6-SE gate for sampled validation, never a fit.
    validation = [r for r in short_rows if r["readout_error"] == 0]
    max_standardized = max(abs(r["state_fidelity"]-r["ideal_state_fidelity"]) /
                          np.sqrt(max(r["ideal_state_fidelity"]*(1-r["ideal_state_fidelity"]), 1/r["shots"])/r["shots"])
                          for r in validation)
    if max_standardized > 6:
        raise AssertionError(f"Short circuit disagrees with ideal model: {max_standardized:.2f} SE")
    progress("Short dynamic circuits passed. Running independent Qiskit calibration batches.", flush=True)

    long_start = perf_counter()
    comparisons, calibration_rows, batch_rows, exact_rows = [], [], [], []
    for pi, (theta, t) in enumerate(((0.1, 128), (0.05, 929))):
        for mi, model in enumerate(("F", "R")):
            evaluation = []
            for b in range(config.calibration_repeats):
                logical_start = perf_counter()
                rec, failure = sample_logical_records(theta, t, model, config.evaluation_shots, seed(21, pi, mi, b))
                accounting["logical_sampling_seconds"] += perf_counter()-logical_start
                evaluation.append((rec, failure))
                np.savez_compressed(raw / f"logical_{pi}_{mi}_{b}.npz",
                                    nonzero_count=rec.nonzero_count, failure=failure)
                accounting["logical_evaluation_records"] += config.evaluation_shots
            for budget in config.budgets:
                all_scores = []
                for b, (rec, failure) in enumerate(evaluation):
                    syn, seconds, count = run_calibration(theta, model, budget, seed(22, pi, mi, budget, b))
                    fit = fit_model(int(np.count_nonzero(syn)), budget, theta)
                    np.savez_compressed(raw / f"cal_{pi}_{mi}_{budget}_{b}.npz", syndromes=syn)
                    accounting["qiskit_calibration_shots"] += budget
                    accounting["qiskit_circuits"] += count
                    accounting["qiskit_execution_and_transpile_seconds"] += seconds
                    calibration_rows.append(dict(theta=theta, true_model=model, budget=budget, batch=b,
                        nonzero=fit.nonzero, selected_model=fit.model, correct=int(fit.model == model),
                        log_odds=fit.log_odds, exact_nonzero_rate=nonzero_probability(theta, model),
                        intervals=2*budget, seed=seed(22, pi, mi, budget, b), seconds=seconds))
                    scoring_start = perf_counter()
                    scores = np.column_stack([score(failure, d.actions(rec)) for d in _decoders(theta, model, fit)])
                    accounting["decoder_scoring_seconds"] += perf_counter()-scoring_start
                    all_scores.append(scores)
                    batch_rows.append(dict(theta=theta, cycles=t, true_model=model, budget=budget, batch=b,
                        selected_model=fit.model, **{f"fe_{n}":float(x) for n, x in zip(NAMES, scores.mean(axis=0))}))
                bootstrap_start = perf_counter()
                stats = paired_hierarchical_bootstrap(all_scores, config.bootstrap_repeats,
                                                     seed(31, pi, mi, budget))
                accounting["bootstrap_seconds"] += perf_counter()-bootstrap_start
                prob_f = exact_selection_probability(theta, model, budget)
                expected_cal = prob_f * exact_fidelity(theta, t, model, "F") + (1-prob_f)*exact_fidelity(theta, t, model, "R")
                expected = (exact_fidelity(theta, t, model), exact_fidelity(theta, t, model, "R"),
                            expected_cal, exact_fidelity(theta, t, model, model))
                selection = [r["correct"] for r in calibration_rows if r["theta"] == theta
                             and r["true_model"] == model and r["budget"] == budget]
                for j, name in enumerate(NAMES):
                    comparisons.append(dict(theta=theta, cycles=t, true_model=model, budget=budget,
                        calibration_intervals_per_batch=2*budget, decoder=name,
                        fe=float(stats["mean"][j]), ci_low=float(stats["ci"][j, 0]), ci_high=float(stats["ci"][j, 1]),
                        f_mean=float(mean_fidelity(stats["mean"][j])), exact_expected_fe=expected[j],
                        gain_vs_standard=float(stats["difference_mean"][j]),
                        gain_ci_low=float(stats["difference_ci"][j, 0]),
                        gain_ci_high=float(stats["difference_ci"][j, 1]),
                        model_selection_accuracy=float(np.mean(selection)),
                        exact_selection_accuracy=prob_f if model == "F" else 1-prob_f,
                        batches=config.calibration_repeats,
                        evaluation_shots_per_batch=config.evaluation_shots,
                        engine="sampled validated logical model; calibration Qiskit"))
                progress(f"Calibration/evaluation: theta={theta}, model={model}, N={budget}", flush=True)
        for tt in sorted(set([0, 1, 2, 4, t] + list(range(8, 2*t+1, max(1, t//48))))):
            for model in ("F", "R"):
                for label, assumed in (("standard", None), ("independence", "R"), ("oracle", model)):
                    fe = exact_fidelity(theta, tt, model, assumed)
                    exact_rows.append(dict(theta=theta, cycles=tt, model=model, decoder=label,
                                           fe=fe, f_mean=float(mean_fidelity(fe))))
        # Every integer time resolves the bare memory's fast coherent revivals.
        # The encoded curves use the smaller sweep above.
        for tt in range(2*t+1):
            for model in ("F", "R"):
                bare = bare_fidelity(theta, tt, model)
                exact_rows.append(dict(theta=theta, cycles=tt, model=model, decoder="bare",
                                       fe=bare, f_mean=float(mean_fidelity(bare))))
    write_csv(output / "decoder_comparison.csv", comparisons)
    write_csv(output / "calibration_batches.csv", calibration_rows)
    write_csv(output / "decoder_batches.csv", batch_rows)
    write_csv(output / "exact_storage_sweep.csv", exact_rows)
    long_seconds = perf_counter()-long_start
    from .validation import audit_samples
    sample_audit = audit_samples(output)

    resources = []
    for t in (0, 1, 2, 4, 128, 929):
        # Compile resource estimates only; long physical circuits are NOT run.
        resources.append(dict(kind="storage", cycles=t,
            executed=(t <= 4), **resource_counts(storage_circuit(0.1, (1,)*t))))
    resources.append(dict(kind="calibration", cycles=2, executed=True,
                          **resource_counts(calibration_circuit(0.1, (1, 1)))))
    write_json(output / "resources.json", resources)
    from .plots import make_plots
    make_plots(output, exact_rows, comparisons, six, resources)
    packages = {name: importlib.metadata.version(name) for name in
                ("qiskit", "qiskit-aer", "numpy", "scipy", "matplotlib",
                 "pytest", "ipykernel", "nbformat", "nbclient")}
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), config=asdict(config),
                    python=sys.version, executable=sys.executable, prefix=sys.prefix,
                    base_prefix=sys.base_prefix, platform=platform.platform(), packages=packages,
                    seeds=used_seeds, accounting=accounting,
                    seconds_total=perf_counter()-start,
                    seconds_long_stage_including_qiskit_calibration=long_seconds,
                    sample_audit=sample_audit,
                    short_validation_max_standardized_residual=float(max_standardized),
                    bootstrap="95% percentile; resample independent batches, then paired evaluation score vectors",
                    main_regime="known angle; ideal repeated X correction; no postselection",
                    noise_regimes=dict(
                        main="engineered common-sign Rx on F and R; no Aer NoiseModel",
                        calibration="two engineered intervals; no Aer NoiseModel",
                        sensitivity="engineered Rx plus q7-only symmetric ReadoutError p=.01/.03; T=4",
                        long="classically sampled validated logical channel; no Aer or device noise"),
                    scope="Qiskit dynamic circuits T=0,1,2,4; T=128,929 sampled logical model",
                    hardware="none; all resource depths are simulator compilation estimates")
    write_json(output / "manifest.json", manifest)
    _report(output, comparisons, six, resources, manifest)
    progress(f"Complete: {output} ({manifest['seconds_total']:.1f} s)", flush=True)
    return manifest


def _report(output, comparisons, six, resources, manifest):
    lines = ["# Executed MVP results", "",
             "Known-angle coherent X-noise experiment with otherwise ideal operations. Long storage uses the Qiskit-validated logical model; "
             "only 0, 1, 2, and 4 repeated-correction cycles were executed as full quantum circuits.", "",
             "**F and R both contain injected rotations.** Main storage, bare controls and calibration use "
             "Aer without an added NoiseModel. Only the short sensitivity test adds ancilla report flips. "
             "See [noise settings and I/O](../docs/noise-and-io.md).", "",
             "| theta | cycles | true model | calibration shots/batch | decoder | observed Fe [95% CI] | exact expected Fe | gain over standard [95% CI] |",
             "|---:|---:|:---:|---:|:---|:---|---:|:---|"]
    for r in comparisons:
        if r["budget"] in (0, 128, 512) and r["decoder"] in ("standard", "calibrated", "oracle"):
            if r["decoder"] != "calibrated" and r["budget"] != 512:
                continue
            lines.append(f"| {r['theta']} | {r['cycles']} | {r['true_model']} | {r['budget']} | {r['decoder']} | "
                         f"{r['fe']:.4f} [{r['ci_low']:.4f}, {r['ci_high']:.4f}] | {r['exact_expected_fe']:.6f} | "
                         f"{r['gain_vs_standard']:.4f} [{r['gain_ci_low']:.4f}, {r['gain_ci_high']:.4f}] |")
    lines += ["", "Fe is entanglement fidelity. For the ideal I/X channel it also equals Z-input survival. "
              "The six-state mean is (2 Fe + 1)/3. The standard and independence-assuming decisions "
              "coincide at these parameters.", "",
              "The same held-out records score all four decoders, and the same evaluation batches are reused "
              "across calibration budgets. Confidence intervals include empirical calibration-batch variation "
              "and paired evaluation shot resampling. They are approximate and can be conservative; if all "
              "batches select the same model they do not reveal unobserved rare misclassification. "
              "The exact finite-binomial expectations shown beside the observations provide that check. "
              "In the full 24-batch release run for F at theta=.1, T=128, N=128, all batches selected correctly: observed Fe=.8238 "
              "but expected Fe=.807603. At N=512, expected Fe=.822941. "
              "Calibration budgets are per batch, amortized over separate held-out records.", "",
              "## Short-circuit sensitivity", "",
              "| true model | ancilla report flip | mean six-state fidelity [conservative 95% CI] | ideal prediction |",
              "|:---:|---:|:---|---:|"]
    for r in six:
        if r["theta"] == 0.1 and r["cycles"] == 4 and r["kind"] == "encoded" and r["decoder"] == "standard":
            lines.append(f"| {r['model']} | {r['readout_error']:.2%} | {r['f_mean']:.4f} "
                         f"[{r['ci_low']:.4f}, {r['ci_high']:.4f}] | {r['ideal_f_mean']:.4f} |")
    lines += ["", "Only syndrome reports have the added flip probability; other gates, reset, and final readout "
              "remain ideal. The reduced model is not used to predict this disturbed circuit. "
              "The oracle label in sensitivity rows means an ideal-model oracle, not an optimal noisy-circuit decoder. "
              "Calibration probes remain ideal, so this does not test calibration robustness.", "",
              "## Cost and scope", "",
              f"- Qiskit storage shots: {manifest['accounting']['qiskit_storage_shots']:,}.",
              f"- Qiskit calibration shots: {manifest['accounting']['qiskit_calibration_shots']:,} "
              f"({2*manifest['accounting']['qiskit_calibration_shots']:,} physical noise intervals).",
              f"- Qiskit bare-memory shots: {manifest['accounting']['qiskit_bare_shots']:,}.",
              f"- Unique long logical-model evaluation records: {manifest['accounting']['logical_evaluation_records']:,}.",
              f"- Total pipeline time: {manifest['seconds_total']:.1f} seconds on the recorded host.",
              f"- Long-model sampling: {manifest['accounting']['logical_sampling_seconds']:.3f} s; "
              f"decoder scoring: {manifest['accounting']['decoder_scoring_seconds']:.3f} s; "
              f"bootstrap: {manifest['accounting']['bootstrap_seconds']:.2f} s. "
              "These are host computation times, not physical storage times.",
              "- Storage uses seven data qubits plus one reusable ancilla. Each cycle has 12 CX, 3 resets, "
              "3 syndrome measurements, and at most one executed conditional X. The final decoder is a "
              "destructive-readout Pauli frame with no added quantum gates.",
              "- Calibration at N shots costs 2N noise intervals, N encoder preparations, 12N extraction CX, "
              "3N resets and 3N measurements. Full compiled totals, including encoder/decoder, are in resources.json.",
              "- No QPU, routing, gate timing, general fault tolerance, unknown-angle estimation, bath witness, "
              "or long physical-circuit performance was established.", "",
              "## Track 3 interpretation", "",
              "The known answer is preservation of an arbitrary encoded qubit. Zero-noise controls, noisy storage, "
              "a fixed correction circuit, and calibrated final recovery are supplied with measured and modeled "
              "costs. Syndrome agreement alone does not establish recovery reliability in this constructed example. "
              "Calibration is worth its cost when it changes a consequential decision; it gives little benefit "
              "at short times or when the independence assumption already holds. This is a simulator demonstration "
              "of evidence needed for trust, not a claim of device readiness or a new error-correction theorem.", "",
              "Reproduction commands and mathematical assumptions are in ../README.md and ../docs/math-note.md. "
              "The proposed path to devices and application metrics is in "
              "[hardware and workloads](../docs/hardware-workloads.md)."]
    (output/"REPORT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
