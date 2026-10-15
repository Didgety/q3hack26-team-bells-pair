"""Export figures from the actual run tables, with explicit engine labels."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_plots(output, sweep, comparison, six, resources):
    output = Path(output)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    colors = {"standard": "#d36a35", "independence": "#d36a35",
              "oracle": "#147d92", "bare": "#8b8d98", "calibrated": "#4b60b3"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.6), sharey=True, constrained_layout=True)
    for row, theta in enumerate((0.1, 0.05)):
        for col, model in enumerate(("F", "R")):
            ax = axes[row, col]
            for decoder in ("bare", "standard", "oracle"):
                values = [r for r in sweep if r["theta"] == theta and r["model"] == model and r["decoder"] == decoder]
                ax.plot([r["cycles"] for r in values], [r["fe"] for r in values],
                        label="standard = independence" if decoder == "standard" else decoder,
                        color=colors[decoder], alpha=0.65 if decoder == "bare" else 1,
                        linewidth=1.0 if decoder == "bare" else 2)
            ax.set(title=f"theta = {theta:.2f} | {'fixed' if model == 'F' else 'refreshed'} sign",
                   xlabel="Storage intervals", ylabel="Entanglement fidelity", ylim=(0, 1.04))
            ax.grid(alpha=0.18)
    axes[0, 0].legend(loc="lower left", fontsize=9)
    fig.suptitle("Exact validated logical model — long physical circuits were not executed", fontsize=13)
    fig.savefig(output/"storage_fidelity.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.6), sharey=True, constrained_layout=True)
    for row, theta in enumerate((0.1, 0.05)):
        for col, model in enumerate(("F", "R")):
            ax = axes[row, col]
            for decoder in ("standard", "calibrated", "oracle"):
                values = [r for r in comparison if r["theta"] == theta and r["true_model"] == model and r["decoder"] == decoder]
                x = np.arange(len(values))
                ys = np.array([r["fe"] for r in values])
                low, high = np.array([[r["ci_low"], r["ci_high"]] for r in values]).T
                ax.errorbar(x, ys, yerr=np.maximum(0, [ys-low, high-ys]),
                            marker="o", markersize=4, capsize=3, color=colors[decoder], label=decoder)
                if decoder == "calibrated":
                    ax.plot(x, [r["exact_expected_fe"] for r in values], "k:", alpha=0.55, label="exact calibrated expectation")
                ax.set_xticks(x, [r["budget"] for r in values])
            ax.set(title=f"theta = {theta:.2f} | {'fixed' if model == 'F' else 'refreshed'} sign",
                   xlabel="Qiskit calibration shots (2 intervals each)",
                   ylabel="Entanglement fidelity", ylim=(0.18, 0.92))
            ax.grid(alpha=0.18)
    axes[0, 0].legend(loc="lower right", fontsize=8)
    fig.suptitle("Separate temporal calibration + held-out logical-model evaluation | 95% paired bootstrap", fontsize=12)
    fig.savefig(output/"calibration_budget.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for model, color in (("F", "#147d92"), ("R", "#b86530")):
        values = [r for r in six if r["theta"] == 0.1 and r["cycles"] == 4
                  and r["kind"] == "encoded" and r["decoder"] == "standard" and r["model"] == model]
        ys = np.array([r["f_mean"] for r in values])
        low, high = np.array([[r["ci_low"], r["ci_high"]] for r in values]).T
        axes[0].errorbar([100*r["readout_error"] for r in values], ys,
                        yerr=[ys-low, high-ys], marker="o", capsize=3, color=color, label=model)
    axes[0].set(title="Executed 4-cycle Qiskit circuits", xlabel="Ancilla report flip probability (%)",
                ylabel="Six-state mean fidelity", ylim=(0.85, 1.005))
    axes[0].legend(title="Temporal model")
    axes[0].grid(alpha=0.18)
    cal = next(r for r in resources if r["kind"] == "calibration")
    budgets = sorted(set(r["budget"] for r in comparison))
    positions = np.arange(len(budgets))
    axes[1].bar(positions, [cal["cx"]*n for n in budgets], color="#4b60b3")
    axes[1].set_xticks(positions, budgets)
    axes[1].set(title="Compiled calibration cost per batch",
                xlabel="Calibration shots", ylabel="CX operations (including encoder)")
    axes[1].text(0.04, 0.94, f"{cal['qubits']} qubits | {cal['cx']} CX/shot | 2 intervals/shot\nAll-to-all; no device timing",
                 transform=axes[1].transAxes, va="top", fontsize=9)
    fig.suptitle("Sensitivity and overhead — ideal reset, gates, and final readout", fontsize=13)
    fig.savefig(output/"cost_sensitivity.png")
    plt.close(fig)
