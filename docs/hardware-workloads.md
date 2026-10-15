# Connection to hardware and useful workloads

This experiment asks whether we can preserve a logical state during **storage**.
Its result is a controlled example of a diagnostic changing a recovery decision.
It does not yet establish an improvement in a chemistry, optimization, or
simulation algorithm. The following device and workload paths are proposals;
the factual examples are cited separately.

## Why the hardware connection is concrete

The circuit uses the same basic ingredients needed for feedback: ancilla
measurement, reset, and a conditional operation before the quantum state is
measured finally. IBM's current dynamic-circuit interface supports mid-circuit
measurement and classical feedforward through Qiskit's if_test. Its documentation
also explains that measurement and classical operations can take longer than
two-qubit gates. Backend capabilities and measurement instructions must therefore
be checked when adapting this circuit.
[IBM dynamic-circuit guide](https://quantum.cloud.ibm.com/docs/en/guides/execute-dynamic-circuits).

Repeated QEC and classical decoding already matter in physical quantum memories.
Google Quantum AI demonstrated below-threshold surface-code memory, including
a distance-five experiment with real-time decoding. That is hardware evidence
for a different code and noise regime; our ideal Steane fidelities are not a
comparison with its error-per-cycle results.
[Nature paper](https://www.nature.com/articles/s41586-024-08449-y).

Our inference from these examples: a decoder belongs in a measurement-and-control
system whose accuracy, data access, throughput, and response time all matter.
In this MVP the final decoder runs offline, and its choice only changes the
interpretation of destructive readout. It has not been put into a real-time
controller or used to preserve a qubit for subsequent logical gates.

## What maps, and what still has to be measured

| Current MVP | Device interpretation | Evidence required before a device claim |
|---|---|---|
| Eight virtual qubits, all-to-all compilation | Seven data sites and one ancilla with repeated interactions | Inspect native gate set, final layout, routing and measurement map; include routing gates |
| Programmed Rx(2 xi theta) | Controlled coherent error injection for a diagnostic benchmark | Separate intended rotations from native errors; estimate/drift-track theta |
| Ideal Z-check extraction | Repeated parity readout and classical correction | Characterize gate, measurement and reset faults; add full X/Z protection and fault-tolerant gadgets for general noise |
| Report flips on q7 | One limited part of measurement uncertainty | Add duration-dependent decoherence and measured correlations; validate calibration with those imperfections too |
| F/R model selection | A small temporal model-selection problem | Test whether these two candidate distributions explain device records; expand or reject the model if they do not |
| 70 + 12T compiled CX, 3T + 1 measurements | Resource accounting before routing | Schedule on a target; report measurement/reset/feedforward time, not only gate depth |
| N calibration shots | A reusable diagnostic budget | Measure drift and specify a validity window; charge 47N CX and 2N intervals plus actual device time |
| Final I/X readout frame | Classical interpretation of measured output | For a continuing computation, track the frame through later logical operations or implement recovery |

The 47N CX calibration cost includes a fresh 35-CX encoder plus 12-CX extraction
per probe. Its benefit can be amortized over multiple storage shots only while
the relevant noise behavior is stable. The full study uses 24 calibration batches
and 4096 held-out records per batch; N is not a per-evaluation-shot charge.

There is no physical interval duration in the present theta parameter. Assigning
seconds, T1/T2 values, or a hardware error probability to it would require a new
model and measurements. Similarly, compiled instruction depth is not latency.

## Example workload: chemistry and many-body energies

Sample-based quantum diagonalization (SQD) prepares and samples quantum states,
then uses classical configuration recovery and diagonalization in the sampled
subspace to estimate Hamiltonian eigenvalues. IBM documents chemistry and
fermionic lattice applications of this workflow.
[IBM SQD documentation](https://quantum.cloud.ibm.com/docs/en/addons/qiskit-addon-sqd).

Our proposed transfer is the **validation method**: diagnose a noise assumption
independently, freeze the chosen processing rule, then score fresh data against
a meaningful target and count the extra cost. SQD's configuration recovery and
this project's logical recovery are different algorithms. This decoder has not
been integrated into SQD, and improved memory Fe does not establish improved
energy estimates or quantum advantage.

A concrete extension would use a small, classically solvable Hamiltonian and a
fixed state-preparation circuit. Compare energy error and uncertainty at matched
total shot/time budgets, with independent calibration and held-out samples.
Measure whether any diagnostic is still useful once algorithm gates, native
noise, and sampling overhead are present. A positive memory result alone cannot
answer that question.

## Proposed next device experiment

1. Select an accessible backend and inspect its dynamic-circuit capabilities.
   Map only the short T = 0, 1, 2, 4 circuits first.
2. Save native circuit layout, measurement map, schedule, and calibration timestamp.
   Recount actual gates and time; make suppression/twirling settings explicit.
3. Run ideal and device-derived simulation controls before requesting QPU access.
   Add calibration-probe faults as well as storage faults, then revalidate the
   reduced model or report its failure.
4. Interleave fresh probes and held-out storage under a stated drift window.
   Evaluate all decoders on the same data without giving them hidden labels.
5. Only after this works, introduce a workload with a known classical answer and
   its own accuracy metric. Keep a matched-cost unprotected baseline.

Current platform facts above were checked on 2026-10-14.
