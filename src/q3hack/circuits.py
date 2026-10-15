"""Qiskit dynamic circuits and measured short-experiment shot records.

There are two independent noise controls. ``theta`` programs coherent RX
gates on the data (the intended experiment); ``readout_error`` optionally
adds an Aer NoiseModel for ancilla measurement reports. With the default
readout_error=0, Aer adds no device noise, but theta != 0 still injects
coherent noise. Only theta=0 AND readout_error=0 is the no-noise control.
Calibration and bare controls always use the simulator without a NoiseModel.

Seven data qubits plus reusable ancilla q7 implement the X-correction
sector of Steane memory. This is a memory workload, not a general algorithm
or a hardware run. See math-note.md §§1–2, 5–7 and code-walkthrough.md.
"""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, ReadoutError
from .code import H, RECOVERY_QUBIT, encoder, preparation
from .decoders import SyndromeRecord
from .noise import sign_mixture


def extract(qc, register, correct=True):
    """Append three Z-check measurements and optional physical X correction.

    Mutates ``qc`` in place and returns None. ``qc`` must contain data
    q0..q6, ancilla q7, and the three-bit ClassicalRegister ``register``.
    For check a, resetting q7 and using data-to-ancilla CX gates measures
    row a of H into register[a]. Thus its integer value is s0+2*s1+4*s2.

    ``correct=True`` appends the table-based R_s = X**RECOVERY[s], using
    the reported register value. At most one conditional X executes.
    If an Aer readout error corrupts that value, the physical correction
    can be wrong; the ideal reduced logical model no longer applies.
    """
    for a, support in enumerate(H):
        qc.reset(7)
        for q, bit in enumerate(support):
            if bit:
                qc.cx(q, 7)
        qc.measure(7, register[a])
    if correct:
        for syndrome, q in RECOVERY_QUBIT.items():
            with qc.if_test((register, syndrome)):
                qc.x(q)


def storage_circuit(theta, signs, state="0", injected_x=None, final_x=False):
    """Build, but do not execute, one encoded storage sign assignment.

    Args:
        theta: Programmed half-rotation angle in radians, per interval.
            RX(phi)=exp(-i*phi*X/2), so phi=2*sign*theta.
        signs: Length-T sequence of ±1, common to all data qubits in each
            interval. This is hidden generator information, not decoder input.
        state: Input label from code.STATES; used for preparation and scoring.
        injected_x: Optional data-qubit index 0..6 for a deterministic X
            fault in the first interval; a validation hook, not an Aer error.
        final_x: If True, physically apply logical X on all seven data qubits
            before decoding. Normal decoder comparisons use a readout frame
            in ShotRecords.score instead, on identical measured records.

    Returns:
        An eight-qubit QuantumCircuit with 3*T+1 classical bits. Each cycle
        does rotation -> extraction -> fixed correction; inverse encoding
        and inverse preparation make final output 0 signify input survival.
        Classical registers are s0, ..., s(T-1), out. Qiskit's display order
        reverses them: "out s(T-1) ... s0", with each chunk "s2 s1 s0".
        The circuit itself specifies no backend NoiseModel or hardware layout.
    """
    q = QuantumRegister(8, "q")
    registers = [ClassicalRegister(3, f"s{t}") for t in range(len(signs))]
    out = ClassicalRegister(1, "out")
    qc = QuantumCircuit(q, *registers, out)
    qc.compose(preparation(state), [0], inplace=True)
    qc.compose(encoder(), range(7), inplace=True)
    qc.barrier()
    for t, sign in enumerate(signs):
        # The factor of two translates the paper's exp(-i*xi*theta*X)
        # convention into Qiskit's Bloch-sphere RX rotation convention.
        qc.rx(2 * sign * theta, range(7))
        if injected_x is not None and t == 0:
            qc.x(injected_x)
        extract(qc, registers[t])
        qc.barrier()
    if final_x:
        qc.x(range(7))
    qc.compose(encoder().inverse(), range(7), inplace=True)
    qc.compose(preparation(state).inverse(), [0], inplace=True)
    qc.measure(0, out[0])
    return qc


def calibration_circuit(theta, signs):
    """Build a |0_L> temporal probe with two uninterrupted noise intervals.

    ``theta`` is the half-rotation angle in radians; ``signs`` contains
    exactly two ±1 values. There is no extraction between the intervals:
    equal signs add to ±2*theta and opposite signs cancel. The returned
    eight-qubit circuit measures only one three-bit syndrome, without
    correction or final logical readout. It supplies evidence about F/R
    that ordinary per-interval syndrome records cannot supply (§5).
    This builder does not choose a simulator or add a NoiseModel.
    """
    if len(signs) != 2:
        raise ValueError("Calibration uses exactly two uninterrupted intervals")
    q, s = QuantumRegister(8, "q"), ClassicalRegister(3, "s")
    qc = QuantumCircuit(q, s)
    qc.compose(encoder(), range(7), inplace=True)
    qc.barrier()
    for sign in signs:
        qc.rx(2 * sign * theta, range(7))
        qc.barrier()
    extract(qc, s, correct=False)
    return qc


def bare_circuit(theta, signs, state="0"):
    """Build the one-qubit unencoded control for a named Pauli input.

    ``theta`` is in radians and ``signs`` lists the hidden ±1 intervals.
    There is no encoding, syndrome extraction, or correction. The returned
    circuit measures one inverse-prepared output bit: 0 means survival.
    The noise process has the same single-qubit marginal as encoded storage.
    """
    qc = QuantumCircuit(1, 1)
    qc.compose(preparation(state), inplace=True)
    for sign in signs:
        qc.rx(2 * sign * theta, 0)
        qc.barrier()
    qc.compose(preparation(state).inverse(), inplace=True)
    qc.measure(0, 0)
    return qc


def simulator(readout_error=0.0):
    """Return Aer statevector simulation, optionally with noisy syndrome reports.

    ``readout_error`` is a probability p in [0, 1], not an angle. At p=0
    no Aer NoiseModel is attached: gates, reset, and readout are ideal,
    while programmed coherent rotations still execute. At p>0 the sole
    NoiseModel entry independently flips measurement reports on ancilla
    q7 with matrix [[1-p,p],[p,1-p]]. Final q0 readout remains ideal.

    This is a limited readout sensitivity model, not a calibrated device
    model: it has no relaxation, dephasing, gate error, or routing noise.
    Statevector is the simulation method; it does not mean that the supplied
    circuit is noiseless. Measurement, reset, and optional report errors
    are sampled during shots. Random seeds are set by _run_mixture.
    """
    kwargs = dict(method="statevector", max_parallel_threads=1)
    if readout_error:
        noise = NoiseModel()
        p = readout_error
        # Only the ancilla's reported syndrome is corrupted. Final readout,
        # encoding, decoding, data rotations, and reset are ideal here.
        noise.add_readout_error(ReadoutError([[1-p, p], [p, 1-p]]), [7])
        kwargs["noise_model"] = noise
    return AerSimulator(**kwargs)


@dataclass
class ShotRecords:
    """Measured circuit outputs plus run metadata, kept outside the decoder.

    Attributes:
        syndromes: uint8 array (shots, cycles), integers 0..7, chronological.
            Bare controls have shape (shots, 0), since they have no checks.
        outcome: uint8 array (shots,), the inverse-prepared final bit.
            0 means return to ``state`` before a final decoder action.
        state: Pauli-input label used only for scoring the measurement basis.
        seconds: Wall-clock seconds for transpilation, execution, and merging;
            not a simulated storage duration or hardware latency.
        circuit_count: Number of sign-assignment circuits that actually ran
            with positive shot counts, not a count of shots or QEC cycles.
    """
    syndromes: np.ndarray
    outcome: np.ndarray  # 0 = successful return in the chosen input basis
    state: str
    seconds: float
    circuit_count: int

    def decoder_input(self):
        """Return only T and per-shot counts of nonzero syndrome integers.

        This is sufficient for the symmetric ideal branch model. It counts
        cycles with s != 0, not set bits, physical faults, or logical failures.
        Input labels and final readout outcomes remain inaccessible to Decoder.
        Under added report noise this is still the deployed interface, but
        ideal-model sufficiency/optimality is no longer asserted.
        """
        return SyndromeRecord(self.syndromes.shape[1], np.count_nonzero(self.syndromes, axis=1))

    def score(self, actions):
        """Return float 0/1 survival scores after final logical Pauli decisions.

        ``actions`` is a boolean array (shots,): False=I and True=logical X.
        X reverses the Z/Y readout outcome and preserves X-basis survival.
        This destructive-readout frame is equivalent to the final Pauli for
        these observables; it does not leave a corrected live quantum state.
        Averaging scores estimates state fidelity for the prepared input.
        """
        # Destructive-readout Pauli frame; logical X flips Z/Y but not X.
        flip = np.asarray(actions) if self.state not in ("+", "-") else False
        return np.logical_not(np.logical_xor(self.outcome, flip)).astype(float)


def _run_mixture(circuits, weights, shots, seed, backend):
    """Draw iid mixture membership, group for execution, erase injection labels.

    Random multinomial counts (rather than fixed stratified allocations) keep
    calibration shots iid Bernoulli after the hidden schedules are forgotten.

    ``circuits`` and normalized ``weights`` have one entry per hidden sign
    assignment. ``shots`` is a nonnegative total count. ``seed`` initializes
    NumPy's allocations/shuffle and derived Aer seeds; compilation always
    uses seed_transpiler=7. ``backend`` selects whether report noise is added.
    Returns (Qiskit memory strings, elapsed seconds, executed circuit count).
    The returned strings contain observations only, never schedule labels.
    Zero shots returns ([], 0.0, 0). A failed backend result raises an error.
    """
    if shots == 0:
        return [], 0.0, 0
    rng = np.random.default_rng(seed)
    allocations = rng.multinomial(shots, weights)
    start = perf_counter()
    compiled = transpile(circuits, backend, optimization_level=0, seed_transpiler=7)
    memory, executed = [], 0
    for qc, count in zip(compiled, allocations):
        if count:
            result = backend.run(qc, shots=int(count), memory=True,
                                 seed_simulator=int(rng.integers(1, 2**31))).result()
            if not result.success:
                raise RuntimeError(result.status)
            memory.extend(result.get_memory())
            executed += 1
    # Forget grouping as well as signs, preserving independent shot membership.
    rng.shuffle(memory)
    return memory, perf_counter()-start, executed


def run_storage(theta, cycles, model, state, shots, seed, readout_error=0.0):
    """Execute short repeated-QEC circuits and return ShotRecords.

    ``theta`` is radians, ``cycles`` is interval count T, ``model`` is the
    injected F/R distribution, and ``state`` is one of code.STATES. ``shots``
    is the total over sign assignments and ``seed`` controls their sampling
    and Aer outcomes. ``readout_error`` is the optional ancilla-report flip
    probability (0 means no Aer NoiseModel, not no coherent injection).

    This runs actual dynamic Qiskit circuits. Refreshed-sign enumeration is
    limited by noise.sign_mixture; the released experiment runs T<=4.
    Parsing keeps each syndrome in chronological order and the final bit
    separately. Only decoder_input() is passed to final recovery rules.
    """
    schedules, weights = sign_mixture(model, cycles)
    circuits = [storage_circuit(theta, signs, state) for signs in schedules]
    memory, elapsed, count = _run_mixture(circuits, weights, shots, seed, simulator(readout_error))
    values = np.array([int(bits.replace(" ", ""), 2) for bits in memory], dtype=np.uint64)
    # The first-created classical register occupies the low bits despite
    # appearing at the right of Qiskit's displayed memory string.
    syn = np.column_stack([(values >> (3*t)) & 7 for t in range(cycles)]).astype(np.uint8) if cycles else np.empty((shots, 0), np.uint8)
    outcome = ((values >> (3*cycles)) & 1).astype(np.uint8)
    return ShotRecords(syn, outcome, state, elapsed, count)


def run_calibration(theta, model, shots, seed):
    """Execute independent two-interval probes without an Aer NoiseModel.

    ``theta`` is radians; ``model`` (F/R) is known only to the generator.
    ``shots`` is the probe budget N and ``seed`` controls shot generation.
    Returns (syndromes, seconds, circuit_count), where syndromes is uint8
    shape (N,) containing integers 0..7. Count nonzero entries to obtain K
    for calibration.fit_model(K, N, theta); do not pass the true model.
    Deliberate coherent rotations remain present when theta != 0. Callers
    must use calibration seeds independent of evaluation seeds.
    """
    schedules, weights = sign_mixture(model, 2)
    circuits = [calibration_circuit(theta, signs) for signs in schedules]
    memory, elapsed, count = _run_mixture(circuits, weights, shots, seed, simulator())
    syndromes = np.array([int(bits.replace(" ", ""), 2) for bits in memory], dtype=np.uint8)
    return syndromes, elapsed, count


def run_bare(theta, cycles, model, state, shots, seed):
    """Execute the unencoded control without an Aer NoiseModel.

    Inputs have the same meanings as run_storage, except that no readout
    noise option is applied. Returns ShotRecords with measured final bits
    and an empty syndrome axis. Explicit RX noise is present for theta != 0.
    These records are a separate control, not inputs to encoded decoders.
    """
    schedules, weights = sign_mixture(model, cycles)
    circuits = [bare_circuit(theta, signs, state) for signs in schedules]
    memory, elapsed, count = _run_mixture(circuits, weights, shots, seed, simulator())
    outcome = np.array([int(bits, 2) for bits in memory], dtype=np.uint8)
    return ShotRecords(np.empty((shots, 0), np.uint8), outcome, state, elapsed, count)


def resource_counts(circuit):
    """All-to-all compilation. Branch bodies counted separately, not executed.

    Depth is compiled instruction-layer depth and includes control-flow boxes;
    it is not a device duration. CX gates are unconditional for this circuit.

    ``circuit`` is a QuantumCircuit. The returned dict gives qubits,
    classical_bits, cx, resets, measurements, rx, conditional_x_branches,
    compiled_depth, and explanatory convention strings. Counts describe
    one execution of the supplied protocol; multiply by shots for a budget.
    Seven X branches are present per cycle but at most one executes.
    This function compiles only: long T resource rows are not long Qiskit
    runs, and no target device connectivity, timing, or error is supplied.
    """
    qc = transpile(circuit, basis_gates=["h", "s", "sdg", "x", "z", "rx", "cx", "reset", "measure", "if_else"],
                   optimization_level=0, seed_transpiler=7)
    ops = qc.count_ops()
    conditional_x = sum(block.count_ops().get("x", 0)
                        for item in qc.data if item.operation.name == "if_else"
                        for block in item.operation.blocks)
    return dict(qubits=qc.num_qubits, classical_bits=qc.num_clbits,
                cx=int(ops.get("cx", 0)), resets=int(ops.get("reset", 0)),
                measurements=int(ops.get("measure", 0)), rx=int(ops.get("rx", 0)),
                conditional_x_branches=int(conditional_x),
                compiled_depth=int(qc.depth()),
                depth_convention="instruction layers including control-flow boxes; no hardware timing",
                connectivity="all-to-all simulator; no hardware layout")
