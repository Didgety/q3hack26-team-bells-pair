# Executed MVP results

Known-angle coherent X-noise experiment with otherwise ideal operations. Long storage uses the Qiskit-validated logical model; only 0, 1, 2, and 4 repeated-correction cycles were executed as full quantum circuits.

**F and R both contain injected rotations.** Main storage, bare controls and calibration use Aer without an added NoiseModel. Only the short sensitivity test adds ancilla report flips. See [noise settings and I/O](../docs/noise-and-io.md).

| theta | cycles | true model | calibration shots/batch | decoder | observed Fe [95% CI] | exact expected Fe | gain over standard [95% CI] |
|---:|---:|:---:|---:|:---|:---|---:|:---|
| 0.1 | 128 | F | 0 | calibrated | 0.3832 [0.3794, 0.3871] | 0.384512 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | F | 128 | calibrated | 0.8238 [0.8202, 0.8275] | 0.807603 | 0.4407 [0.4347, 0.4471] |
| 0.1 | 128 | F | 512 | standard | 0.3832 [0.3793, 0.3872] | 0.384512 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | F | 512 | calibrated | 0.8238 [0.8203, 0.8275] | 0.822941 | 0.4407 [0.4347, 0.4466] |
| 0.1 | 128 | F | 512 | oracle | 0.8238 [0.8203, 0.8275] | 0.823052 | 0.4407 [0.4347, 0.4466] |
| 0.1 | 128 | R | 0 | calibrated | 0.6101 [0.6057, 0.6143] | 0.611117 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | R | 128 | calibrated | 0.6101 [0.6057, 0.6140] | 0.604249 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | R | 512 | standard | 0.6101 [0.6056, 0.6146] | 0.611117 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | R | 512 | calibrated | 0.6101 [0.6056, 0.6146] | 0.611076 | 0.0000 [0.0000, 0.0000] |
| 0.1 | 128 | R | 512 | oracle | 0.6101 [0.6056, 0.6146] | 0.611117 | 0.0000 [0.0000, 0.0000] |
| 0.05 | 929 | F | 0 | calibrated | 0.2551 [0.2515, 0.2585] | 0.254879 | 0.0000 [0.0000, 0.0000] |
| 0.05 | 929 | F | 128 | calibrated | 0.7594 [0.6845, 0.8311] | 0.687190 | 0.5043 [0.4281, 0.5752] |
| 0.05 | 929 | F | 512 | standard | 0.2551 [0.2517, 0.2585] | 0.254879 | 0.0000 [0.0000, 0.0000] |
| 0.05 | 929 | F | 512 | calibrated | 0.7598 [0.6672, 0.8307] | 0.805976 | 0.5047 [0.4128, 0.5753] |
| 0.05 | 929 | F | 512 | oracle | 0.8314 [0.8283, 0.8343] | 0.829665 | 0.5763 [0.5706, 0.5819] |
| 0.05 | 929 | R | 0 | calibrated | 0.7428 [0.7388, 0.7465] | 0.743822 | 0.0000 [0.0000, 0.0000] |
| 0.05 | 929 | R | 128 | calibrated | 0.6756 [0.6093, 0.7277] | 0.688986 | -0.0673 [-0.1335, -0.0163] |
| 0.05 | 929 | R | 512 | standard | 0.7428 [0.7388, 0.7467] | 0.743822 | 0.0000 [0.0000, 0.0000] |
| 0.05 | 929 | R | 512 | calibrated | 0.7101 [0.6590, 0.7446] | 0.728451 | -0.0328 [-0.0830, 0.0000] |
| 0.05 | 929 | R | 512 | oracle | 0.7428 [0.7388, 0.7467] | 0.743822 | 0.0000 [0.0000, 0.0000] |

Fe is entanglement fidelity. For the ideal I/X channel it also equals Z-input survival. The six-state mean is (2 Fe + 1)/3. The standard and independence-assuming decisions coincide at these parameters.

The same held-out records score all four decoders, and the same evaluation batches are reused across calibration budgets. Confidence intervals include empirical calibration-batch variation and paired evaluation shot resampling. They are approximate and can be conservative; if all batches select the same model they do not reveal unobserved rare misclassification. The exact finite-binomial expectations shown beside the observations provide that check. In the full 24-batch release run for F at theta=.1, T=128, N=128, all batches selected correctly: observed Fe=.8238 but expected Fe=.807603. At N=512, expected Fe=.822941. Calibration budgets are per batch, amortized over separate held-out records.

## Short-circuit sensitivity

| true model | ancilla report flip | mean six-state fidelity [conservative 95% CI] | ideal prediction |
|:---:|---:|:---|---:|
| F | 0.00% | 0.9848 [0.9767, 0.9896] | 0.9835 |
| R | 0.00% | 0.9849 [0.9768, 0.9897] | 0.9847 |
| F | 1.00% | 0.9794 [0.9704, 0.9851] | 0.9835 |
| F | 3.00% | 0.9574 [0.9457, 0.9659] | 0.9835 |
| R | 1.00% | 0.9784 [0.9692, 0.9842] | 0.9847 |
| R | 3.00% | 0.9565 [0.9447, 0.9651] | 0.9847 |

Only syndrome reports have the added flip probability; other gates, reset, and final readout remain ideal. The reduced model is not used to predict this disturbed circuit. The oracle label in sensitivity rows means an ideal-model oracle, not an optimal noisy-circuit decoder. Calibration probes remain ideal, so this does not test calibration robustness.

## Cost and scope

- Qiskit storage shots: 258,048.
- Qiskit calibration shots: 264,192 (528,384 physical noise intervals).
- Qiskit bare-memory shots: 208,896.
- Unique long logical-model evaluation records: 393,216.
- Total pipeline time: 84.6 seconds on the recorded host.
- Long-model sampling: 0.052 s; decoder scoring: 0.242 s; bootstrap: 2.54 s. These are host computation times, not physical storage times.
- Storage uses seven data qubits plus one reusable ancilla. Each cycle has 12 CX, 3 resets, 3 syndrome measurements, and at most one executed conditional X. The final decoder is a destructive-readout Pauli frame with no added quantum gates.
- Calibration at N shots costs 2N noise intervals, N encoder preparations, 12N extraction CX, 3N resets and 3N measurements. Full compiled totals, including encoder/decoder, are in resources.json.
- No QPU, routing, gate timing, general fault tolerance, unknown-angle estimation, bath witness, or long physical-circuit performance was established.

## Track 3 interpretation

The known answer is preservation of an arbitrary encoded qubit. Zero-noise controls, noisy storage, a fixed correction circuit, and calibrated final recovery are supplied with measured and modeled costs. Syndrome agreement alone does not establish recovery reliability in this constructed example. Calibration is worth its cost when it changes a consequential decision; it gives little benefit at short times or when the independence assumption already holds. This is a simulator demonstration of evidence needed for trust, not a claim of device readiness or a new error-correction theorem.

Reproduction commands and mathematical assumptions are in ../README.md and ../docs/math-note.md. The proposed path to devices and application metrics is in [hardware and workloads](../docs/hardware-workloads.md).
