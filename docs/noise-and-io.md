# Noise settings, inputs, and outputs

The short answer: **both F and R use deliberately injected coherent noise.**
The main Aer runs have no added simulator NoiseModel. A separate short-circuit
sensitivity test adds a NoiseModel that flips the reported ancilla bit.
An ideal simulator can execute a noisy *experiment* when error rotations are
explicit gates in the circuit.

## Which engine and which noise?

| Experiment | Execution engine | Programmed error | Added Aer noise | What it establishes |
|---|---|---|---|---|
| Ideal controls | Qiskit Aer shots | theta = 0 | None | Known-answer baseline, including encoder and correction |
| Main short encoded memory, T = 0, 1, 2, 4 | Qiskit Aer dynamic circuits | Seven Rx(2 xi theta) gates per interval, F or R | None | Actual measurement, reset, and conditional correction match the ideal branch model |
| Bare memory | Qiskit Aer shots | One Rx(2 xi theta) per interval, F or R | None | Unencoded control under the same temporal noise rule |
| Calibration probes | Separate Qiskit Aer shots | Two uninterrupted intervals, then extraction | None | Learn F versus R at a known angle |
| Report-error sensitivity, T = 4, theta = 0.1 | Qiskit Aer dynamic circuits | Same F/R rotations | Symmetric 1% or 3% ReadoutError on ancilla q7 only | Show fragility to imperfect syndrome reports |
| Long storage, T = 128 or 929 | Classical sampling of the Qiskit-validated logical model | F/R encoded in conditional logical channels | No Aer execution or device NoiseModel | Compare four decoders on identical held-out records |
| Exact curves and branch checks | Algebraic enumeration; independent Qiskit Statevector check | Ideal coherent model | None | Reference predictions, without shot uncertainty |

In the sensitivity experiment, a flipped ancilla report can trigger the wrong
physical X correction. Final data readout, gates, and resets remain ideal.
The calibration probes themselves are not given report errors. This test therefore
does **not** establish robustness of calibration to device noise.

The project contains no depolarizing gate model, T1/T2 relaxation, leakage,
crosstalk, hardware-derived noise model, or QPU data.

## Inputs

| Input | Meaning and units | Who can see it? |
|---|---|---|
| theta | Known coherent amplitude angle, in radians; Qiskit's Rx gate angle is 2 xi theta | Circuit builder and model-based decoders |
| xi = +1 or -1 | Fair hidden sign, common to all seven data qubits in an interval | Injection code only; never a decoder input |
| F / R | F retains a single sign for an entire shot; R redraws it each interval | Evaluator knows the true label; calibrated decoder learns it from separate probes |
| T / cycles | Number of rotation-plus-extraction-and-correction intervals | Every decoder; not seconds or circuit depth |
| state | Logical input: 0, 1, +, -, +i, or -i | Preparation and scoring only |
| N / budget | Number of calibration shots in **one** batch | Calibrator; N is not added to each held-out storage shot |
| short_shots | Shot budget per short-circuit condition, default 2048 | Experiment controller |
| calibration_repeats | Independent calibration/evaluation batches, default 24 | Statistics; captures finite calibration variability |
| evaluation_shots | Held-out logical records per batch, default 4096 | Statistics; reused across decoders and budgets |
| readout_error | Probability of a symmetric reported-bit flip on q7, 0 / .01 / .03 | Short-circuit sensitivity engine only |
| seed | Reproducibility control; disjoint namespaces for each stage | Sampling code; not a physical parameter |

F and R are *both* spatially correlated: every data qubit shares the sign in an
interval. Only sign persistence across time differs. R does not mean noise-free,
and F does not mean a fixed deterministic positive angle.

## Outputs and interpretation

| Output | Meaning |
|---|---|
| syndrome | Integer 0–7 from three measured Z checks, s0 + 2 s1 + 4 s2; identifies a correction class, not a revealed hidden error |
| syndromes | Array shaped (shots, T) for short encoded runs |
| nonzero_count / k | Number of nonzero syndromes in a record; sufficient here because all seven nonzero branches share a logical form |
| outcome | Final measurement after undoing preparation; 0 means input-basis survival |
| failure | Sampled Z-basis failure bit for long-model scoring; evaluator only |
| fitted model, likelihood, log odds | Evidence from separate two-interval probes for F versus R under the known-angle, two-candidate model |
| decoder action | Final I or logical X, implemented as a destructive-readout Pauli frame; no reusable quantum output is produced |
| state_fidelity | Fraction of an input's measurements that survive in its own basis |
| f_mean | Mean fidelity across all six Pauli inputs |
| Fe | Entanglement fidelity of the logical channel; equals Z-input survival only in the ideal I/X channel here |
| exact finite-calibration expectation | Model average over every possible binomial calibration outcome; exposes rare mistakes absent from sampled batches |
| confidence interval | Finite-shot uncertainty conditional on this model and experiment; not a bound on device mismatch |
| resources / depth | Compiled all-to-all counts and instruction layers; not scheduled duration or QPU cost |

For a trace-preserving qubit channel, Favg = (2 Fe + 1)/3; the six Pauli inputs
give this average. These are different quantities and are labeled separately.
The four decoders are standard (final I), independence (assume R), calibrated
(select F/R using probes), and oracle (knows the distribution but not hidden signs).

At theta = .1, T = 128, N = 128, observed calibrated Fe is .8238 with bootstrap
CI [.8202, .8275], but exact finite-calibration expectation is .807603. All 24
observed batches chose correctly, so that CI misses rare wrong-model batches.
At N = 512, the corresponding exact expectation is .822941, close to observed
.8238. The result files retain both expectations and observations.

See [the code walkthrough](code-walkthrough.md) for equations mapped to functions,
[the data dictionary](../results/README.md) for file columns and record IDs, and
[hardware and workloads](hardware-workloads.md) for the proposed device path.
