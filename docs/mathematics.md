# Mathematical construction for the MVP

This compact construction develops the Steane instance and its temporal-decoding
rules. [Steane's code construction](https://arxiv.org/abs/quant-ph/9601029) and
[Morales-Hashemi, Smith, and Cao](https://arxiv.org/abs/2610.03545v1) provide the
published foundations; full citations are in [references.bib](../references.bib).
The [mathematical note](math-note.md) gives the detailed implemented conventions,
branch polynomials, uncertainty analysis, and executed checks. Section 7 below
adds an explicit-environment construction and extension considerations that
were not implemented in this experiment.

## 1. Explicit code and ordering

Number data qubits \(0,\ldots,6\). Let each column below correspond to that order:

\[
H=\begin{pmatrix}
0&0&0&1&1&1&1\\
0&1&1&0&0&1&1\\
1&0&1&0&1&0&1
\end{pmatrix}.
\]

The six stabilizer generators are \(X^{h_a}\) and \(Z^{h_a}\), for rows \(h_a\), with positive signs. All row overlaps have even parity, so the generators commute. Let \(C_X=\operatorname{rowspan}_{\mathbb F_2}(H)\), containing eight binary strings, and \(\ell=(1,1,1,1,1,1,1)\). Fix

\[
|0_L\rangle=\frac1{\sqrt8}\sum_{v\in C_X}|v\rangle,
\qquad |1_L\rangle=X^{\ell}|0_L\rangle,
\qquad \bar X=X^{\ell},\quad \bar Z=Z^{\ell}.
\]

Implement a Clifford encoder for an arbitrary input qubit, not a different state-preparation recipe for each test input. Validate all six logical Pauli eigenstates. Under the core X-only disturbance, the three Z checks detect the reachable syndrome sectors. The other three checks remain trivial in the ideal model. For generic noise this simplification need not hold.

For a binary X support \(e\), \(s(e)=He\pmod2\). Choose \(r_0=0\); for nonzero syndrome \(s\), choose the unique single-qubit support whose column equals \(s\). The per-cycle correction is \(R_s=X^{r_s}\). The syndrome value in the reference check is \(s_0+2s_1+4s_2\). Qiskit displays Pauli labels and classical strings in a different visual order; document the mapping explicitly.

## 2. Engineered coherent noise and memory

For one physical interval define

\[
U_\xi(\theta)=\bigotimes_{j=0}^6e^{-i\xi\theta X_j},
\qquad \xi\in\{+1,-1\}.
\]

Qiskit implements each factor as `rx(2 * xi * theta, j)`. Thus \(\theta\) is half the Bloch rotation angle. A physical Pauli-twirled factor has X probability \(\sin^2\theta\), but twirling every interval changes a coherent temporal process. It is not an interchangeable representation for this experiment.

Define two run-level models:

- \(\mathsf F\): draw one fair sign and keep it for all intervals of one run.
- \(\mathsf R\): draw a fair sign independently in every interval, shared by all seven data qubits in that interval.

Both have the same isolated one-interval channel

\[
\mathcal N(\rho)=\tfrac12U_+\rho U_+^\dagger+
\tfrac12U_-\rho U_-^\dagger.
\]

The refreshed model still has spatial correlations. The experiment varies temporal persistence while keeping the spatial rule fixed. Runs are independent conditional on the model. Keeping a random sign in a Python variable is a valid implementation of classical temporal memory; it is not proof that a quantum bath is required.

Each data interval is followed by syndrome extraction and the fixed \(R_s\). Record every syndrome. The main comparison changes only the final logical recovery, not the extraction schedule.

## 3. Conditional logical Kraus operators

Let \(V\) denote the encoder isometry and \(\Pi_s\) the syndrome projector. Define

\[
K_s(\xi)=V^\dagger R_s\Pi_sU_\xi V.
\]

Expand \(U_\xi\) as a coherent sum over all 128 X supports:

\[
U_\xi=\sum_{e\in\{0,1\}^7}c_e(\xi)X^e,
\qquad c_e(\xi)=(\cos\theta)^{7-|e|}(-i\xi\sin\theta)^{|e|}.
\]

Group amplitudes before squaring. After applying the correction, \(e\oplus r_s\) lies in either \(C_X\) or \(\ell\oplus C_X\). Therefore

\[
K_s(\xi)=\alpha_s(\xi)I+\beta_s(\xi)\bar X,
\]

\[
\alpha_s=\sum_{s(e)=s,\ e\oplus r_s\in C_X}c_e,
\qquad
\beta_s=\sum_{s(e)=s,\ e\oplus r_s\in\ell\oplus C_X}c_e.
\]

Multiplying \(K_s\) by the irrelevant branch phase \((i\xi)^{|r_s|}\) makes the identity coefficient real and the X coefficient imaginary. Write that representative as \(a_sI-i\xi b_s\bar X\), with \(a_s,b_s\) real and independent of the sign. Then

\[
K_s^\dagger K_s=q_sI,\qquad q_s=a_s^2+b_s^2,
\qquad \sum_s q_s=1.
\]

This directly proves input-independent, sign-independent syndrome probabilities for this instance. Corrected conditional evolution is a logical X rotation, with multiplier on X-basis coherence

\[
g_s(\xi)=\frac{(a_s-i\xi b_s)^2}{q_s},
\qquad |g_s|=1,\quad g_s(-)=g_s(+)^*.
\]

The sign of an intermediate coefficient depends on our basis and phase convention. Use the explicit amplitude sums and compare operators, not just signed constants printed in another convention. Zero-probability branches are omitted.

## 4. History-dependent final recovery

For \(y=(s_1,\ldots,s_T)\), both models give

\[
p(y)=\prod_{t=1}^Tq_{s_t}.
\]

Their conditional logical coherence multipliers differ:

\[
\lambda_{\mathsf F}(y)=\operatorname{Re}\prod_{t=1}^Tg_{s_t}(+),
\qquad
\lambda_{\mathsf R}(y)=\prod_{t=1}^T\operatorname{Re}g_{s_t}(+).
\]

The first averages the hidden sign only after its coherent accumulation. The second averages each fresh sign separately. Here both multipliers are real; the resulting conditional channel is dephasing about logical X.

For final logical recovery \(I\) or \(\bar X\), the conditional entanglement fidelities are

\[
f_h(y,I)=\tfrac12(1+\lambda_h(y)),\qquad
f_h(y,\bar X)=\tfrac12(1-\lambda_h(y)).
\]

Thus the model-aware Pauli rule chooses I for nonnegative \(\lambda_h\) and logical X otherwise. The oracle knows h and \(\theta\), not the realized signs. A calibrated rule plugs in allowed estimates. The mean fidelity is \(F_e(h,d)=\sum_yp(y)f_h(y,d(y))\). For this real dephasing channel the same maximum is attainable with general deterministic CPTP recovery; this special-case observation is not a universal statement about optimal Pauli recovery.

Keep the distinction between the fixed per-cycle table \(R_s\) and the final history-dependent logical action \(d(y)\). A classifier based only on syndrome records cannot identify h here because their complete distributions coincide.

For long runs exploit this instance's symmetry: all seven nonzero syndromes have the same q and g. Group histories by nonzero count k with probability \(\binom Tkq_0^{T-k}(1-q_0)^k\). Use log-probabilities or stable recurrences. This avoids enumerating \(8^T\) histories. Record why this shortcut applies; disable it if a perturbation breaks symmetry.

## 5. Separate calibration with delayed extraction

Prepare a known encoded state, apply TWO noise intervals without an intervening extraction, then measure the three Z checks. Let

\[
Q(\theta)=1-q_0(2\theta),
\]

where \(q_0(x)\) is the one-extraction zero-syndrome probability at rotation parameter x. Equal signs double the angle; opposite signs cancel. Thus the nonzero-syndrome calibration event has probability

\[
p_{\rm nz}(\mathsf F)=Q(\theta),\qquad
p_{\rm nz}(\mathsf R)=\tfrac12Q(\theta).
\]

With N independent calibration shots and k events, the known-angle log-likelihood ratio is

\[
\log\frac{L_{\mathsf F}}{L_{\mathsf R}}=
k\log2+(N-k)\log\frac{1-Q}{1-Q/2}.
\]

Use an explicit prior and deterministic tie rule. Main evaluation runs retain ordinary per-interval correction. Count two noise intervals per delayed calibration shot. If \(\theta\) is unknown, estimate it using separate one-interval probes and jointly fit the temporal model; do not secretly pass true \(\theta\) to the learned rule. Known-angle and unknown-angle experiments must be labeled separately. Do not assert a universal shot-cost law from this finite experiment.

## 6. Fidelity, uncertainty, and independent checks

Use six logical Pauli eigenstates as the Qiskit end-to-end test ensemble:

\[
\overline F=\tfrac16\sum_{\psi\in\{0,1,+,-,+i,-i\}}
\langle\psi|\rho_{\rm out}^{(\psi)}|\psi\rangle.
\]

For a trace-preserving logical-qubit channel and a state-independent recovery, \(\overline F=(2F_e+1)/3\). State-specific recovery or postselection needs its own treatment. A reference qubit can measure entanglement fidelity directly, at the cost of another physical qubit.

The model gives these **analytical check values, not measured Qiskit or hardware results**:

| theta | cycles | Oracle F: Fe | Oracle R: Fe | R-assuming decoder evaluated on F: Fe |
|---:|---:|---:|---:|---:|
| 0.10 | 128 | 0.8230515273 | 0.6111174360 | 0.3845121731 |
| 0.05 | 929 | 0.8296646494 | 0.7438220766 | 0.2548793741 |

These points establish a decision-sensitive ideal regime. They do not promise feasible hardware storage times or that finite calibration achieves the oracle value. The small-angle expansion gives \(1-\operatorname{Re}\sum_sq_sg_s\sim126\theta^4\) and \(|\operatorname{Im}\sum_sq_sg_s|\sim28\theta^3\). The mathematical note provides the branch polynomials and normalization/sign-conjugacy checks.

For sampled fidelity and calibration comparisons use paired evaluation and repeated independent calibration batches. Bootstrap both sources of uncertainty. Separate shot noise, parameter-estimation variability, and exact-model error. Compute exact quantities before adding statistical estimation, so the uncertainty protocol cannot hide a model mismatch.

## 7. Optional explicit environment

The commuting interaction \(e^{-i\theta Z_B\sum_jX_j}\), with balanced bath Z populations, realizes the two common rotation signs. A retained bath versus a reset/reprepared bath maps to F versus R under the stated ideal controls. The bath must not be measured or supplied to the decoder in evaluation.

For genuinely quantum-memory questions add noncommuting bath dynamics, specify the full joint unitary, and recompute the conditional process through actual syndrome interventions. A product of rotation gates with a bath pulse is a split-step circuit model. It is not exactly a simultaneous Hamiltonian exponential when its terms do not commute. Use that circuit as ground truth or quantify Trotter error. Check that first-interval channels match for the chosen comparison.

Echo and randomization modify the data process. Benchmark their overhead separately. Built-in ZNE targets implementation error around the programmed circuit, including its intentional interaction; it does not automatically erase the engineered disturbance. Formal trace-distance backflow needs two unconditioned initial states and a specified system boundary ([Laine, Piilo, and Breuer](https://arxiv.org/abs/1002.2583)). Neither that witness nor master-equation fitting is implemented in the core decoder experiment.
