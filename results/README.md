# Results and raw records

These files come from the full experiment run. The notebook uses smaller budgets
by default and writes its separate outputs to notebook-run/.

| File | Contents |
|---|---|
| REPORT.md | Results, fidelity definitions, uncertainty, costs, and limitations |
| manifest.json | Run parameters, seed namespaces, package versions, timings, and shot counts |
| branch_validation.json | Independent Qiskit versus algebraic branch-operator checks |
| short_state_results.csv | Per-input survival, Wilson intervals, and raw record IDs |
| short_six_state_results.csv | Six-state means and conservative intervals |
| syndrome_histograms.csv | Measured syndrome frequencies and ideal probabilities |
| conditional_channel_validation.csv | Survival conditioned on nonzero syndrome count |
| calibration_rates.csv | Measured delayed-probe rates and exact predictions |
| short_calibration.json | Separate calibration used for short-circuit decoder comparisons |
| calibration_batches.csv | Independent calibration counts, fitted models, seeds, and costs |
| decoder_batches.csv | Per-batch held-out performance for all four decoders |
| decoder_comparison.csv | Observed fidelities, paired bootstrap intervals, and exact finite-calibration expectations |
| exact_storage_sweep.csv | Exact encoded and bare-memory curves |
| resources.json | All-to-all compiled gate counts and depths; long circuits are compiled only |
| *.png | Plots of storage fidelity, calibration budgets, and sensitivity/cost |
| records/*.npz | Raw arrays readable with numpy.load, without pickle |
| pytest.txt, pip_check.txt | Output of the test suite and dependency check |

For short-state rows, record_id identifies records/<record_id>.npz. Encoded
records contain syndromes shaped (shots, cycles) and outcome shaped (shots,).
Outcome 0 means survival in the prepared input's measurement basis. Bare-memory
records contain outcomes only. Calibration records contain three-bit syndrome
integers. Long-model records contain nonzero_count and failure arrays.

Decoders receive storage length and syndrome count, together with their allowed
known or learned model. Final measurement outcomes are used only for scoring.
No raw record exposes the injected hidden signs to a decoder.

Fe denotes entanglement fidelity; f_mean is the six-state average fidelity.
For a trace-preserving qubit channel, f_mean = (2 Fe + 1)/3. Exact expectations
are model calculations, not extra measured shots. Long-storage records are
samples of the validated logical model, not long Aer or hardware executions.
