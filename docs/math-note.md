# Temporal calibration for a Steane quantum memory

This note specifies the implemented experiment and derives its finite-angle predictions.
The code, amplitude enumeration, Clifford encoder synthesis, Qiskit cross-checks, sampling,
and resource tables are this project's calculations. The underlying code and temporal-decoding
phenomenon are published results, not a novelty claim. The experiment uses Qiskit
to test a Steane memory under restricted coherent X noise.

## Sources and claim boundaries

Steane's simplex-code/coset construction supplies the encoding foundation
(steane1996interference, §3.3, manuscript pp. 17–19). Our physical ordering and CNOT network
are explicitly synthesized below, rather than copied from that circuit diagram.
[Full text](https://arxiv.org/pdf/quant-ph/9601029).

The F/R experiment follows Morales-Hashemi, Smith, and Cao, **v1**, §§IV–V.
Their branch structure and recovery reduction are in Eq. (20) and Appendix E.
We independently evaluate the specified seven-qubit instance, and do not claim to test
their asymptotic cost theorem. The finite-angle regime matters; an asymptotic near-perfect
oracle is not a promise at our angles (morales2026memory).
[Full text](https://arxiv.org/html/2610.03545v1).

Calibration has broader precedent: Hockings, Doherty, and Harper study calibrated decoder
priors for circuit-level Pauli noise in surface-code memories, using ACES. That setting
does not supply our coherent temporal likelihood, and we import no hardware timing or
calibration-cost figure from it (hockings2025noiseaware, main text, “Quantum codes and
decoders” and “Numerical results”).
[Full text](https://arxiv.org/html/2502.21044v2).

Kam et al. compare temporal correlation structures at matched marginal error rates.
Their surface-code study motivates matched controls, not a numerical prediction for
our Steane circuit (kam2025detrimental, v4, §III.1 and §V).
[Full text](https://arxiv.org/html/2410.23779v4).

Nielsen gives the trace-preserving-channel relation between entanglement and average
fidelity (nielsen2002average, v2, Eq. (3), proof pp. 1–2).
[Full text](https://arxiv.org/pdf/quant-ph/0205035v2).
The two-state trace-distance criterion discussed only as an extension comes from
Laine, Piilo, and Breuer (laine2010measure, v2, §II.A, Eqs. (1), (6)–(8)).
[Full text](https://arxiv.org/pdf/1002.2583v2).

The complete bibliography is in [references.bib](../references.bib).
Process-tensor and Petz decoders are not implemented in this experiment.
The [compact mathematical construction](mathematics.md) also discusses an
optional explicit environment and the conditions for extending this model.

## 1. Code, tensor order, and encoder

Data qubits are \(q_0,\ldots,q_6\). Integer bit \(j\) denotes \(q_j\), so statevector index
\(\sum_j2^j v_j\) represents the displayed ket \(|v_6\cdots v_0\rangle\).
Define

\[
H=\begin{pmatrix}
0&0&0&1&1&1&1\\
0&1&1&0&0&1&1\\
1&0&1&0&1&0&1
\end{pmatrix},\quad C_X=\mathrm{rowspan}_{\mathbb F_2}H,\quad \ell=(1,\ldots,1).
\]

The row masks are 120, 102, 85. \(HH^\mathsf T=0\) over \(\mathbb F_2\);
the three \(X^{h_a}\) and three \(Z^{h_a}\) stabilizers therefore commute. They are
independent positive-sign generators, leaving a two-dimensional code space. Set

\[
V|0\rangle=|0_L\rangle=8^{-1/2}\sum_{v\in C_X}|v\rangle,\qquad
V|1\rangle=|1_L\rangle=X^\ell|0_L\rangle.
\]

The disjoint cosets establish orthonormality. Row weights are four and their span has
seven weight-four nonzero words. Thus \(Z^\ell\) is +1 on the first coset and −1 on
the second; \(X^\ell\) exchanges them. These are \(\bar Z,\bar X\).
The seven distinct nonzero columns of \(H\) make the minimum nonzero weight in
\(\ker H\) equal to three; its nontrivial logical coset has seven weight-three words.
The symmetric X/Z construction gives distance three. We test only the X-correction
part of the code; the full six-check general-noise protocol is not being executed.

The encoder begins with arbitrary input on \(q_0\), zeros elsewhere, and H gates on
\(q_1,q_2,q_3\). An invertible binary map whose first four columns are
\((\ell,h_0,h_1,h_2)\) maps the basis variables to
\(x\ell+u_0h_0+u_1h_1+u_2h_2\). Completing these columns to a basis and reversing
Gaussian row elimination synthesizes a CNOT/SWAP circuit. This proves the same
Clifford encoder accepts every input \(a|0\rangle+b|1\rangle\).
Qiskit tests check all six Pauli eigenstates, the six stabilizers, and both logical
operators. The resulting simple, unoptimized network costs 35 CX after SWAP
decomposition; a decode costs another 35.

For physical X support \(e\), syndrome \(s=He\) is stored as
\(s_0+2s_1+4s_2\). The table is:

| Syndrome integer | Corrected data qubit |
|---:|---:|
| 0 | none |
| 1 | 3 |
| 2 | 1 |
| 3 | 5 |
| 4 | 0 |
| 5 | 4 |
| 6 | 2 |
| 7 | 6 |

Let \(r_s\) be that weight-zero/one mask. Then \(R_s=X^{r_s}\).
Each physical cycle measures Z checks using data-to-ancilla CNOTs and ancilla
Z measurement, resetting the same eighth qubit before every check. A distinct
three-bit register stores each cycle. Qiskit's displayed memory string is
“out s(T-1) ... s0”; each displayed three-bit chunk is “s2 s1 s0”.
The parser removes spaces, reads a binary integer, extracts low three-bit chunks
in chronological order, and reads the final output at bit \(3T\).
The all-single-X tests exercise this mapping, not just the table in isolation.

## 2. Specified noise and memory

One interval is

\[
U_\xi(\theta)=\bigotimes_{j=0}^6(\cos\theta\,I-i\xi\sin\theta\,X_j),
\qquad \xi\in\{-1,+1\}.
\]

The Qiskit instruction is rx(2*xi*theta, j). The angle \(\theta\) is in radians,
half the Bloch rotation. Model F draws one fair sign per run; model R independently
draws a fair common sign each interval. In both cases the sign is shared across
the seven data qubits in an interval. This is spatial correlation in both models,
and classical temporal memory only in F. Both have the same single-interval channel
\(\mathcal N=\tfrac12\operatorname{Ad}_{U_+}+\tfrac12\operatorname{Ad}_{U_-}\).

The data circuit performs \(U_\xi\), three checks, and \(R_s\) every interval.
It is identical for all final decoders. Sign schedules are internal to circuit
generation. For short circuits, F has two assignments and R has \(2^T\).
We draw multinomial assignment counts for each shot batch, execute each selected
assignment in Aer, combine and shuffle the measured records, and discard assignment
labels. This samples the specified iid run mixture; it does not supply hidden signs
to a decoder or substitute Pauli errors for coherent rotations.

A per-qubit twirl would yield X probability \(\sin^2\theta\), but changes the temporal
experiment. No twirl, echo, bath measurement, or postselection is performed.

## 3. Independently enumerated coherent branches

Expand

\[
U_\xi=\sum_{e=0}^{127}c_eX^e,\qquad
c_e=c^{7-|e|}(-i\xi u)^{|e|},\quad c=\cos\theta,\ u=\sin\theta.
\]
After projection and correction \(e\oplus r_s\) belongs either to \(C_X\) or
\(\ell+C_X\). Hence

\[
K_s(\xi)=V^\dagger R_s\Pi_sU_\xi V
=\alpha_s I+\beta_s X,
\]
where \(\alpha_s,\beta_s\) sum the corresponding amplitudes **before** squaring.
An irrelevant phase \((i\xi)^{|r_s|}\) makes this
\(\widetilde K_s=a_sI-i\xi b_sX\). Our enumeration gives the following explicit
polynomials:

\[
\begin{aligned}
a_0&=c^7+7c^3u^4,&b_0&=-(7c^4u^3+u^7),\\
a_{s\ne0}&=c^6u-4c^4u^3+3c^2u^5,&
b_{s\ne0}&=3c^5u^2-4c^3u^4+cu^6.
\end{aligned}
\]

For example, the identity-class weights in a nonzero syndrome are 1 (one support),
3 (four supports), and 5 (three supports); the logical-X class has complementary
weights 6, 4, and 2. This establishes the factors and signs without fitting.

Because \(a_s,b_s\) are real, the cross term cancels:

\[
K_s^\dagger K_s=(a_s^2+b_s^2)I=q_sI,\qquad \sum_s q_s=1.
\]
The syndrome probability is independent of input and sign. For \(q_s>0\),
the normalized channel is an X rotation, whose multiplier on the
\(|+\rangle\langle-|\) coherence is

\[
g_s(\xi)=\frac{(a_s-i\xi b_s)^2}{q_s},
\quad |g_s|=1,\quad g_s(-)=g_s(+)^*.
\]
Zero-probability branches can have arbitrary unit-modulus placeholders because they
carry zero weight. The code uses 1 there. Changing a logical-axis convention can
conjugate the multiplier; probabilities and fidelities must stay invariant.

As a small-angle check, \(G=\sum_s q_sg_s(+)\) satisfies

\[
1-\operatorname{Re}G=2\sum_s b_s^2=126\theta^4+O(\theta^6),\qquad
\operatorname{Im}G=-2\sum_sa_sb_s=-28\theta^3+O(\theta^5).
\]
The leading logical drift is not removed by syndrome averaging.
The existing standard-library reference and new Qiskit check are separate
implementations. The branch_validation.json result records operator residuals.

## 4. Histories, recovery, and information limits

Every syndrome branch has scalar probability, including after prior rotations.
Therefore induction gives \(p(y)=\prod_t q_{s_t}\), under either model and for any
input. Even the complete ideal syndrome history cannot identify F versus R.
Multiplying the normalized branch rotations and averaging over unobserved signs
gives, for \(y=(s_1,\ldots,s_T)\),

\[
\lambda_F(y)=\operatorname{Re}\prod_tg_{s_t}(+),\qquad
\lambda_R(y)=\prod_t\operatorname{Re}g_{s_t}(+).
\]
The conditional channel is

\[
\mathcal L_{\lambda}(\rho)=\tfrac{1+\lambda}{2}\rho
+\tfrac{1-\lambda}{2}X\rho X.
\]
It is a mixture only **after** averaging hidden signs in the specified temporal
order. Treating each physical error path as a classical event earlier loses the
interference that defines this experiment.

Let \(\eta_d(y)=+1\) denote final I and −1 denote final logical X.
Entanglement fidelity conditioned on the history is

\[
f_e(y,d)=\tfrac12[1+\eta_d(y)\lambda_h(y)],\qquad
F_e(h,d)=\sum_y p(y)f_e(y,d).
\]
Thus choose I if the assumed multiplier is nonnegative, X otherwise; ties choose I.
The distribution oracle uses true \(h,\theta\), but never a realized sign or input.
The standard final decoder is I, and the independence decoder always assumes R.
The calibrated decoder plugs in a model selected from separate probe counts.
Its action interface accepts only cycles and nonzero_count, with known angle
and selected model stored in the decoder object. Outcomes and input labels exist
only in scoring. At our two main points the standard and R-assuming rules coincide.

For completeness, any CPTP recovery with Kraus \(A_j\) has
\(f_I=\frac14\sum_j|\operatorname{Tr}A_j|^2\) and
\(f_X=\frac14\sum_j|\operatorname{Tr}A_jX|^2\).
Expansion in the four orthogonal Paulis and trace preservation give
\(f_I+f_X\le1\). Therefore
\(p_I f_I+p_X f_X\le\max(p_I,p_X)\), attained by I or X.
This proves optimality for this real I/X channel only, consistent with the
published reduction above. It is not an assertion about arbitrary coherent noise.

All seven nonzero syndrome branches have the same \(q_n,g_n\).
Only the count \(k\) of nonzero syndromes matters:

\[
p_T(k)={T\choose k}q_0^{T-k}(1-q_0)^k,\quad
\lambda_F(k)=\cos[(T-k)\arg g_0+k\arg g_n],\quad
\lambda_R(k)=(\operatorname{Re}g_0)^{T-k}(\operatorname{Re}g_n)^k.
\]
We sum all \(T+1\) counts with SciPy's stable binomial PMF; no tail is dropped.
Symmetry assertions guard the shortcut. Long records store this sufficient count
and a sampled Z-readout outcome, and are explicitly labeled logical-model records.
No long physical-circuit execution is implied.

## 5. Separate known-angle temporal calibration

Prepare \(|0_L\rangle\); apply two uninterrupted intervals; then measure three checks.
Equal signs double the angle and opposite signs cancel. With
\(Q(\theta)=1-q_0(2\theta)\), the nonzero rate is \(Q\) for F and \(Q/2\) for R.
For \(N\) fresh shots, \(K\) is binomial and the equal-prior log odds are

\[
\log\frac{L_F}{L_R}
=K\log2+(N-K)\log\frac{1-Q}{1-Q/2}.
\]
Positive odds choose F; ties, including \(N=0\), choose R. At \(\theta=0\) neither
model is identifiable and both channels are identity, so the default has no loss.
The known programmed angle is given fairly to all model-based decoders. No
unknown-angle inference is claimed.

Main budgets are 0, 32, 128, 512, 2048 shots, with 24 independent calibration
batches for each angle/model/budget. Every nonempty batch comes from actual Qiskit
two-interval circuits. Evaluation seeds are disjoint and independent.
The simulator knows the injected model to generate data; fit_model receives
only \(K,N,\theta\). Evaluation records are reused across all four decoders and
budgets, so gains are paired. Each evaluation batch has 4096 logical-model shots.

Let \(P_F(N|h)\) be the exact binomial probability that the likelihood rule selects F.
The expected calibrated fidelity is

\[
P_F(N|h)F_e(h,d_F)+[1-P_F(N|h)]F_e(h,d_R).
\]
We include this exact finite-shot expectation alongside empirical estimates.
It exposes rare calibration mistakes that a finite set of batches may miss.
It also explains why finite calibration can hurt under R.
For small angle \(Q=28\theta^2+O(\theta^4)\), but this project reports finite
shot/interval budgets, not a universal optimal cost law.

## 6. Fidelity and uncertainty

Entanglement fidelity uses a maximally entangled input and an untouched reference:

\[
F_e=\langle\Phi|(\mathrm{id}\otimes\mathcal E)(|\Phi\rangle\langle\Phi|)|\Phi\rangle.
\]
For a trace-preserving qubit channel with state-independent recovery,
\(\bar F=(2F_e+1)/3\), by the cited Nielsen relation. A qubit's six Pauli eigenstates
have zero average Bloch vector and isotropic second moment
\(\frac16\sum_{\psi}r_\psi r_\psi^\mathsf T=I/3\), hence their mean survival equals
this average fidelity, even for a nonunital channel. In our ideal X channel the
two X eigenstates survive exactly and each Y/Z eigenstate survives with probability
\(F_e\), independently verifying the relation.

Short circuits prepare and inverse-read all six inputs, including ±i. All
per-state intervals are 95% Wilson intervals; the six-state summary averages
six Bonferroni-adjusted Wilson bounds for a conservative approximate 95% interval.
The noisy-report circuit still defines a trace-preserving decoded qubit channel
because other data qubits are implicitly traced out, with no postselection.
Here we foreground the measured six-state mean, not a claim of continued ideal
code-space confinement. The final readout itself is ideal in this sensitivity test.

For long experiments, logical-Z survival equals \(F_e\) specifically because of
the derived I/X channel. We resample the 24 independent batches and, within each
selected batch, resample common evaluation score vectors. This paired two-level
bootstrap propagates observed calibration variability and shot noise.
The 2.5th/97.5th percentiles of 1500 replicates form approximate 95% intervals.
Resampling both levels can be conservative; empirical intervals cannot quantify
unobserved rare model-selection failures. Exact expectations are the accompanying
model check. Exact sums have no artificial error bars.

These uncertainties are conditional on the specified known-angle simulator.
They do not quantify drift, model misspecification, hardware error, or systematic
calibration bias. No state-specific recovery, postselection, or hidden-sign oracle
is used.

## 7. Bare control, overhead, and verification

For one unencoded qubit with the same marginal sign process,

\[
\lambda_F^{\rm bare}=\cos(2T\theta),\qquad
\lambda_R^{\rm bare}=[\cos(2\theta)]^T,\qquad
F_e^{\rm bare}=(1+\lambda^{\rm bare})/2.
\]
The six-state mean again follows the fidelity relation. Short bare circuits are
executed, and long curves are exact model calculations. A bare coherent memory can
revive; encoding is not uniformly better at every storage time.

The executed storage circuit uses eight qubits, \(3T+1\) classical bits, \(7T\)
coherent rotations, \(12T\) extraction CX, \(3T\) resets and \(3T+1\) measurements.
Encoding/decoding add 70 CX in the unoptimized chosen synthesis. Seven conditional
X branches per cycle are compiled, but at most one executes. All-to-all compiled
instruction depth is recorded, not a device duration or routed resource bound.
Long resource circuits are compiled only. The calibration circuit costs 47 CX
(35 encoder + 12 extraction), 14 rotations, 3 resets, and 3 measurements per shot;
its two noise intervals are counted separately. The final readout frame adds
classical arithmetic and no quantum gates. A live output would require a physical
logical action or a maintained frame; this project ends in destructive measurement.

Verification layers are intentionally distinct:

1. Algebra: 128-support branch enumeration, normalization, conjugacy, symmetry,
   known finite-angle checkpoints, and the existing reference script.
2. Independent physical calculation: Qiskit rotation/project/correct/decode
   matrices for both signs at several angles, and selected multi-cycle histories.
3. Dynamic circuits: measurement, reset, and register-conditioned correction,
   every single-X fault, all six inputs, and 0/1/2/4-cycle tests.
4. End-to-end samples: short fidelities against predictions, actual delayed probes,
   repeated independent calibration, paired held-out decisions, and saved raw records.

The prespecified points reproduce these exact values:

| theta | T | F: standard Fe | F: oracle Fe | R: oracle Fe |
|---:|---:|---:|---:|---:|
| 0.10 | 128 | 0.3845121731 | 0.8230515273 | 0.6111174360 |
| 0.05 | 929 | 0.2548793741 | 0.8296646494 | 0.7438220766 |

Measured calibration performance, confidence intervals, runtimes, and resource
details are generated in [the run report](../results/REPORT.md), with CSV/JSON and
compressed shot records beside it.

## 8. Limitations and proposals

Only the X-error-correction sector (measuring Z checks) with ideal operations supports the analytic model. We add
independent syndrome-report flips of 1% and 3% in short circuits; corrupted reports
can cause wrong physical corrections and residual sectors. The ideal-model
distribution oracle is not a noisy-circuit optimum. We do not extrapolate this
sensitivity to hundreds of cycles or arbitrary Pauli noise.

Known angle, uniform coupling, two candidate models, stationary independent runs,
and ideal state preparation are restrictive. Unknown angle requires fresh probes,
an identifiable parameter region, and uncertainty propagation. Changes that break
nonzero-syndrome symmetry invalidate count compression. Longer noncommuting
environment dynamics require a new physical process calculation.

A retained/reset spectator, physical randomization, echo, process-tensor decoding,
and master equations remain proposals. They change access or the storage process
and need separate comparisons. Classical hidden-sign memory here is not evidence
that a quantum bath is necessary. Trace-distance backflow would require a specified
unconditioned map and two inputs; a single population revival is insufficient.
We measured no such witness. Hardware performance and fault-tolerance thresholds
remain unverified.
