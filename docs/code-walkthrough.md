# How the code, mathematics, and evidence fit together

The experiment stores one arbitrary logical qubit in seven data qubits and asks whether a separately calibrated temporal noise model improves the **final recovery decision**. Both candidate models apply deliberately programmed coherent rotations. Their ordinary syndrome records have the same distribution, yet the appropriate final recovery can differ.

Start here for implementation meaning; use the [mathematical note](math-note.md) for derivations and citations, the [noise and input/output guide](noise-and-io.md) for result interpretation, and the [run report](../results/REPORT.md) for saved results. The canonical implementation is [`src/q3hack`](../src/q3hack), shared by the command-line experiment and notebook.

## 1. What runs, and what “noise” means

“Aer simulator” names the execution engine. It does **not** by itself say which errors are present. The coherent noise is a deliberate circuit operation; the optional Aer `NoiseModel` adds a separate imperfection.

| Evidence path | Engine and entry point | Deliberate coherent noise | Additional simulator noise | What the output establishes |
|---|---|---|---|---|
| No-noise control | Aer shots, `run_storage(theta=0, readout_error=0)` | None | None | Encoding, checks, correction, and decoding preserve the tested inputs |
| Main short memory experiment | Aer dynamic-circuit shots, `run_storage`; released runs have 0, 1, 2, or 4 cycles | Explicit data-qubit `rx(2*sign*theta)` gates in both F and R | None; no Aer `NoiseModel` attached | Measured short-circuit behavior with ideal operations around the injected rotations |
| Readout sensitivity | Aer dynamic-circuit shots, `run_storage(..., readout_error=0.01 or 0.03)` at 4 cycles | Same explicit coherent gates | Symmetric flips of **ancilla q7 reports only** | Short-circuit degradation when reported syndromes can cause wrong corrections |
| Temporal calibration | Aer shots, `run_calibration` | Two uninterrupted rotation intervals before extraction | None; no Aer `NoiseModel` | Observed probe counts used to select F/R |
| Bare short control | Aer shots, `run_bare` | Same marginal sign process on one unencoded qubit | None; no Aer `NoiseModel` | Unencoded state survival for comparison |
| Branch validation | Exact Qiskit `Statevector`, `physical_branches` | Same coherent operator | None; no Aer shot sampling | Independent matrix check of rotation → projection → correction → decoding |
| Exact curves | NumPy/SciPy algebra, `exact_fidelity` / `bare_fidelity` | Included in the derived channel | No additional device noise | Exact predictions **within the specified model**, without sampling error bars |
| Long evaluation | Classical sampling, `sample_logical_records`, at 128/929 modeled cycles | Included in the validated logical channel | No additional device noise | Sampled ideal-model decoder performance; **no long physical-circuit execution** |
| Resource estimates | Qiskit transpilation, `resource_counts` | Circuit instructions counted | No backend noise or timing model | All-to-all compiled gate/depth counts; long circuits are compiled only |

For the readout sensitivity, the final data-qubit readout, encoding/decoding gates, rotations, and resets remain ideal. There is no comprehensive gate-error, relaxation, dephasing, crosstalk, or device calibration model. Its ideal-model “oracle” is a comparison rule, not a proven optimum for the disturbed circuit.

## 2. The vocabulary of inputs and outputs

| Name | Meaning, units, and shape | Who may use it |
|---|---|---|
| `theta` / $\theta$ | Known programmed **half-rotation angle**, radians per interval | Generator and every model-based decoder |
| `cycles` / $T$ | Nonnegative integer number of storage/correction intervals; not seconds | Generator, decoder, and scorer |
| `model`, `true_model` | Distribution label `F` or `R`; not a register or prepared quantum state | Generator/evaluator; the distribution oracle also knows it |
| `assumed_model` | Decoder's hypothesis `F`/`R`; `None` means standard final identity | Decoder; calibration supplies its selected hypothesis |
| `signs`, `xi` / $\xi$ | Hidden ±1 rotation signs; one common sign across all seven data qubits per interval | Circuit generator only; erased from decoder records |
| `state` | Prepared logical input: `0`, `1`, `+`, `-`, `+i`, `-i`, the six Pauli eigenstates | Preparation and scoring; excluded from final decisions |
| `shots` / $N$ | Number of independent executions/records, distinct from cycle count | Sampling and uncertainty calculations |
| `seed` | Reproducible random-number initialization; not a physical input | Sampling only; calibration and evaluation use distinct seeds |
| `readout_error` / $p$ | Dimensionless probability of an ancilla-report flip; default zero | Aer backend configuration only |
| `syndromes` | Integer array `(shots, cycles)` with entries 0–7, in chronological order | Measurement record; condensed before final decoding |
| `nonzero_count` / $k$ | Integer array `(shots,)`; count of cycles with syndrome != 0 | Final decoder; not a count of physical faults or set syndrome bits |
| `outcome` | Array `(shots,)`; 0 means return to the prepared state in its measurement basis, before final action | Scorer only |
| `actions` | Boolean array `(shots,)`; `False` = final logical I, `True` = final logical X | Decoder output, then scorer input |
| `failure` | Boolean array `(shots,)` from **long-model sampling**; Z-input failure before final action | Scorer only; never decoder input |
| `score` | Float 0/1 per shot after final action; 1 means survival | Evaluation; mean is the relevant survival estimate |
| `seconds` | Host wall-clock execution/compilation or classical compute time, depending on the field | Cost reporting; never device storage duration |

`F` means **fixed sign within a run**: draw one fair sign, then keep it for all $T$ intervals. `R` means **refreshed sign**: redraw a fair sign independently each interval. In both models, the sign is common across the seven data qubits in each interval. Thus both have spatial correlation; F also has temporal memory. Neither uses an explicit bath register, and classical hidden-sign memory is sufficient for this example.

The six test-state labels refer to *states*, not noise models. In particular, `+i` means $(|0\rangle+i|1\rangle)/\sqrt2$. The experiment tests arbitrary-state preservation through the six-state benchmark and an arbitrary-input encoder; it is not a useful quantum algorithm executed on application data.

## 3. Follow one storage shot through the code

The path is:

**Prepare → encode → [rotate → extract syndrome → correct] repeated $T$ times → decode and measure → choose a final Pauli using only syndromes → score its readout frame.**

The final Pauli is computed after data collection in this offline implementation. For the measured observables, changing the classical readout frame is equivalent to applying that Pauli before measurement. This procedure does not return a live corrected quantum state for a later gate.

### Encoding and bit order: `code.py`

[`codewords()`](../src/q3hack/code.py) constructs a complex `(128, 2)` matrix $V$ whose columns are $|0_L\rangle$ and $|1_L\rangle$:

$$
|0_L\rangle=\frac1{\sqrt8}\sum_{v\in C_X}|v\rangle,
\qquad |1_L\rangle=X^{\otimes7}|0_L\rangle.
$$

$C_X$ is the binary span of the three rows of the parity-check matrix `H`. The constant named `CX` stores this set of binary words; it is **not** a CNOT gate. `encoder()` synthesizes the same map as a Clifford circuit from an arbitrary input on q0 and zeros on q1–q6. It uses gates rather than initializing a desired encoded state directly. `preparation()` selects the test state; its inverse turns final survival into measurement outcome 0.

An integer support mask uses bit $j$ for q$j$. The displayed basis ket is $|q_6\ldots q_0\rangle$. A syndrome uses $s=s_0+2s_1+4s_2$, where `H[a]` is check $a$. Qiskit displays a three-bit syndrome as `s2 s1 s0`.

### Coherent injection and repeated correction: `noise.py`, `circuits.py`

Each interval applies

$$
U_\xi(\theta)=\bigotimes_{j=0}^{6}
(\cos\theta\,I-i\xi\sin\theta\,X_j).
$$

Qiskit's $R_X(\phi)=e^{-i\phi X/2}$, so the implementation must use `rx(2*sign*theta)`. For example, `theta=0.10` means a ±0.20-radian Bloch rotation per data qubit per interval. It is not a 10% stochastic error rate.

`sign_mixture()` enumerates the allowed hidden schedules and probabilities. `_run_mixture()` draws multinomial shot counts, executes each selected assignment in Aer, combines/shuffles the records, and discards assignment labels. Random allocations preserve independent mixture membership; each assignment still contains coherent gates, rather than sampled Pauli faults.

`extract()` resets q7 before each Z check, accumulates parity through data-to-ancilla CX gates, measures q7, and applies the fixed single-X correction associated with the three-bit report. That **per-cycle correction is identical for every final decoder**.

`run_storage()` parses Qiskit's displayed memory `out s(T-1) ... s0`. The earliest register occupies the lowest three bits. It returns `ShotRecords` containing chronological syndrome integers, final outcomes, state label, elapsed seconds, and the number of sign-assignment circuits actually executed.

### Worked syndrome example: the same report can hide different logical effects

Consider a single X on q0. Its mask is `1` (`0000001`), and column 0 of `H` is $(0,0,1)^T$. Therefore:

| Step | Value | Interpretation |
|---|---:|---|
| Injected support $e$ | `1` | X on q0 |
| `syndrome(1)` | `4` (`100`) | $s_0=0,s_1=0,s_2=1$ |
| `RECOVERY[4]` | `1` | Apply X on q0 |
| Corrected support $e\oplus r_4$ | `0` | Identity logical action |

Now consider support `6` (`0000110`), X on q1 and q2. Its syndrome is also 4, so the same correction is applied. The corrected support is `6 XOR 1 = 7`: X on q0, q1, and q2. It has zero syndrome but is in the **logical-X coset**, since `7 XOR 127 = 120` is a stabilizer row. Thus it has a different logical effect despite the same measured syndrome and correction.

These supports illustrate terms in the coherent expansion; the main simulation does **not** choose one such support as a classical error label. Both and many other supports contribute amplitudes to the same branch. Measuring the syndrome distinguishes branches, but it does not distinguish coherent terms *within* a branch. This is why `branches()` adds amplitudes before squaring.

After a history such as `[4, 0, 4, 0]`, the final decoder receives `cycles=4, nonzero_count=2`. It does not receive the physical supports, hidden signs, prepared state, or final outcome. At the main short-regime angle both model rules can still choose I; the long decision-sensitive points reveal the consequential difference.

## 4. From physical amplitudes to long-memory mathematics

[`branches(theta, xi)`](../src/q3hack/logical_model.py) loops over all 128 X-support masks. A weight-$w$ support contributes

$$
c_e=(\cos\theta)^{7-w}(-i\xi\sin\theta)^w.
$$

The example above explains how correction classifies each support into either the identity or logical-X coset. Summing within each class produces

$$
K_s=\alpha_s I+\beta_s X,
\qquad \widetilde K_s=a_sI-i\xi b_sX,
\qquad q_s=a_s^2+b_s^2.
$$

`kraus` stores the eight raw $(2\times2)$ matrices. The rephased coefficients `a`, `b` are real; the phase changes no physical channel. The scalar identity $K_s^\dagger K_s=q_sI$ makes each syndrome probability independent of logical input and sign. Their total is one. `physical_branches()` independently checks the matrices with Qiskit statevector rotation/project/correct/decode operations.

The normalized branch's effect on $|+\rangle\langle-|$ coherence is

$$
g_s(\xi)=\frac{(a_s-i\xi b_s)^2}{q_s}.
$$

All seven nonzero branches have the same probability and multiplier. `symmetric_parameters()` asserts this equality before compressing a history to count $k$. Under either temporal model,

$$
p_T(k)={T\choose k}q_0^{T-k}(1-q_0)^k.
$$

The temporal distinction enters through **when the hidden signs are averaged**:

$$
\lambda_F(k)=\operatorname{Re}\left(g_0^{T-k}g_n^k\right),
\qquad
\lambda_R(k)=(\operatorname{Re}g_0)^{T-k}(\operatorname{Re}g_n)^k.
$$

`multiplier()` implements these expressions. The conditional logical channel is

$$
\mathcal L_{\lambda}(\rho)=
\frac{1+\lambda}{2}\rho+\frac{1-\lambda}{2}X\rho X.
$$

`count_distribution()` and `exact_fidelity()` sum all $T+1$ counts with no tail truncation. `sample_logical_records()` instead draws $k$ and a logical Z-readout failure from this channel. These are classical long-model records; they contain no sequence of executed physical syndrome measurements. Changes to coupling uniformity or noise that break the branch symmetry require revisiting the reduction, not just changing a label.

## 5. Calibration and final decisions without information leakage

Ordinary ideal syndrome histories cannot distinguish F/R because both have the same $p_T(k)$ and full history distribution. Calibration therefore uses a **different probe schedule**: prepare $|0_L\rangle$, apply two uninterrupted intervals, then extract once.

Two equal signs add, while opposite signs cancel. With $Q=1-q_0(2\theta)$, the nonzero probe rate is $Q$ for F and $Q/2$ for R. [`fit_model(K, N, theta)`](../src/q3hack/calibration.py) receives only the observed nonzero count $K$, probe budget $N$, and known programmed angle. At equal priors its log odds are

$$
\log\frac{L_F}{L_R}=K\log2+(N-K)\log\frac{1-Q}{1-Q/2}.
$$

Positive odds select F; ties select R, including the default zero-shot case. No unknown-angle estimation is performed. The short readout-noise sensitivity still uses ideal calibration probes; it is not a joint fit of readout noise and temporal dynamics.

[`Decoder.actions(record)`](../src/q3hack/decoders.py) returns final X iff its assumed $\lambda(k)<0$. The conditional entanglement fidelity is

$$
f_e(k,d)=\frac{1+\eta_d(k)\lambda_{\mathrm{true}}(k)}2,
\qquad \eta_d=+1\text{ for I},\; -1\text{ for X}.
$$

| Decoder role | Constructor meaning | Allowed information |
|---|---|---|
| Standard | `Decoder(theta, None)` | Always final I after the fixed per-cycle correction |
| Independence | `Decoder(theta, "R")` | Known angle, syndrome count, assumed R |
| Calibrated | `Decoder.calibrated(fit)` | Known angle, syndrome count, separately selected F/R |
| Distribution oracle | `Decoder(theta, true_model)` | Known angle, syndrome count, true distribution; **no realized signs** |

The simulator/evaluator necessarily knows which distribution generated data. That label does not enter the calibrated fit or its decisions. All decoders score the same held-out records; evaluation batches are also reused across calibration budgets for paired comparison.

The pipeline derives seeds from the run seed and separate namespaces: 11 for short storage, 12 for short calibration, 13 for bare controls, 21 for long evaluation, 22 for long-stage Qiskit calibration, and 31 for bootstrapping. Seeds and environment versions are saved in the manifest. Reusing one seed for calibration and evaluation would violate the intended separation.

An illustrative API path, run after the notebook's setup cell has selected the project virtualenv and added src to the import path, is:

```python
import numpy as np
from q3hack.circuits import run_calibration
from q3hack.calibration import fit_model
from q3hack.decoders import Decoder
from q3hack.metrics import sample_logical_records, score

theta = 0.10
# Generator-only truth. Both calls deliberately inject the F distribution.
probe_syndromes, seconds, circuit_count = run_calibration(
    theta, model="F", shots=512, seed=1001)  # actual Aer shots, no NoiseModel
fit = fit_model(int(np.count_nonzero(probe_syndromes)), 512, theta)

record, failure = sample_logical_records(
    theta, cycles=128, model="F", shots=4096, seed=2001)  # classical model
actions = Decoder.calibrated(fit).actions(record)  # no failure/state/sign input
estimated_fe = score(failure, actions).mean()
```

This snippet illustrates the interfaces, not the full release's uncertainty analysis. The complete experiment repeats independent calibration/evaluation batches and uses the manifest's recorded seed scheme.

## 6. What a fidelity and confidence interval mean

For short circuits, `ShotRecords.score()` gives 0/1 survival of the named input after the final readout frame. Logical X flips Z/Y-basis outcomes and leaves X-basis survival unchanged. A mean of these scores is **state fidelity for that input**. Wilson intervals describe its sampled uncertainty.

Entanglement fidelity $F_e$ measures preservation when the stored qubit is entangled with an untouched reference. For a trace-preserving qubit channel with an input-independent recovery,

$$
\overline F=\frac{2F_e+1}{3}.
$$

The equally weighted six Pauli-input survivals yield $\overline F$. `mean_fidelity()` performs this conversion; the notebook/report distinguish it from $F_e$. For the derived ideal I/X channel specifically, Z-input survival equals $F_e$, which justifies scoring long-model Z records as entanglement fidelity. This shortcut is not a universal rule for noisy qubit channels.

`paired_hierarchical_bootstrap()` takes a binary array `(independent batches, evaluation shots, decoders)`. It resamples batches and then paired score vectors, returning decoder means, 95% percentile intervals, and paired differences against decoder column 0. Calibration errors not present in the observed batches cannot appear in an empirical bootstrap; exact finite-binomial expected performance supplies a separate check. The intervals do not cover hardware drift, unknown angle, or an incorrect noise model.

## 7. Connection to hardware and workloads

The workload being modeled is **logical memory**: preserve a qubit, potentially entangled with the rest of a computation, while it waits. This is a component of a larger computation rather than an application algorithm. The implemented measurement/reset/conditional-correction loop and calibrated final recovery are the connection to hardware operation. No application-level benefit or live-device performance has been measured here.

The [hardware and workloads guide](hardware-workloads.md) connects these operations to current platforms and explains the remaining requirements: device layout, native gates, broader error characterization, timing, and integration of a Pauli frame with continued computation. `resource_counts()` provides compiled protocol counts, not those device measurements. Hundreds of modeled cycles cannot be presented as hundreds of demonstrated hardware cycles.

## 8. Function-to-equation reading map

| Code location | Mathematical object | Mathematical note section |
|---|---|---|
| `code.codewords`, `encoder`, `syndrome` | Isometry $V$, binary parity $s=He$, correction $R_s$ | §1 |
| `noise.sign_mixture`, `circuits.storage_circuit` | $U_\xi(\theta)$ and F/R sign distributions | §2 |
| `logical_model.branches`, `physical_branches` | $K_s=V^\dagger R_s\Pi_sU_\xi V$; $q_s,g_s$ | §3 |
| `logical_model.multiplier`, `count_distribution` | $\lambda_F(k),\lambda_R(k),p_T(k)$ | §4 |
| `decoders.Decoder.actions`, `logical_model.exact_fidelity` | Final I/X rule and weighted $F_e$ | §4 |
| `circuits.calibration_circuit`, `calibration.fit_model` | Delayed-probe rates $Q,Q/2$ and binomial log odds | §5 |
| `metrics.mean_fidelity`, `score`, `paired_hierarchical_bootstrap` | $\overline F=(2F_e+1)/3$, survival, paired uncertainty | §6 |
| `circuits.bare_circuit`, `logical_model.bare_fidelity` | Bare-memory $\cos(2T\theta)$ or $\cos(2\theta)^T$ coherence | §7 |
| `circuits.resource_counts` | Per-shot protocol resources; no hardware time conversion | §7 |

Module and function docstrings supply argument contracts beside the implementation. The [mathematical note](math-note.md) supplies the proofs, explicit branch polynomials, assumptions, and primary-source citations.
