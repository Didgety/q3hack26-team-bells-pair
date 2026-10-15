"""Steane [[7,1,3]] encoding and X-error correction conventions.

An integer support mask uses bit j for data qubit q_j; Qiskit displays the
corresponding basis ket as |q6 ... q0>. The three rows of H specify the Z
checks used to locate X errors, and also generate the X-stabilizer group.
This experiment measures those three Z checks, not all six checks needed
for arbitrary Pauli noise. See docs/math-note.md §1 and docs/code-walkthrough.md.
"""
from functools import lru_cache
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

# H has shape (3 checks, 7 data qubits), with arithmetic over GF(2).
H = np.array([[0, 0, 0, 1, 1, 1, 1],
              [0, 1, 1, 0, 0, 1, 1],
              [1, 0, 1, 0, 1, 0, 1]], dtype=np.uint8)
# ROWS = (120, 102, 85). CX denotes the binary code C_X, not a CX gate.
ROWS = tuple(sum(int(b) << j for j, b in enumerate(row)) for row in H)
CX = frozenset(a * ROWS[0] ^ b * ROWS[1] ^ c * ROWS[2]
               for a in (0, 1) for b in (0, 1) for c in (0, 1))
LOGICAL_X = 127  # X on all seven data qubits implements logical X.
STATES = ("0", "1", "+", "-", "+i", "-i")  # Z, X, Y eigenstates.


def syndrome(mask):
    """Return s = H e as an integer s0 + 2*s1 + 4*s2 in [0, 7].

    ``mask`` is the integer X-support e in [0, 127]; bit j means X acts on
    q_j. Each syndrome bit is the parity of overlap with one Z check.
    For example, X on q0 has mask 1 and syndrome 4 (displayed ``100``).
    """
    return sum(((int(mask) & row).bit_count() % 2) << i
               for i, row in enumerate(ROWS))


# R_s = X**RECOVERY[s] is the fixed per-cycle weight-zero/one correction.
# It is distinct from the later decoder's final logical I/X decision.
RECOVERY = {0: 0, **{syndrome(1 << q): 1 << q for q in range(7)}}
RECOVERY_QUBIT = {s: mask.bit_length() - 1 for s, mask in RECOVERY.items() if s}


def codewords():
    """Return the complex (128, 2) encoding isometry V = [|0_L>, |1_L>].

    |0_L> is the equal superposition over C_X; |1_L> uses its disjoint
    all-ones coset. Column labels are logical computational-basis inputs,
    while rows are little-endian physical basis indices. V†V = I_2.
    """
    v = np.zeros((128, 2), dtype=complex)
    for mask in CX:
        v[mask, 0] = 1 / np.sqrt(8)
        v[mask ^ LOGICAL_X, 1] = 1 / np.sqrt(8)
    return v


def _rank(columns):
    """Return the GF(2) rank of binary column vectors stored as int masks."""
    basis = {}
    for mask in columns:
        while mask:
            pivot = mask.bit_length() - 1
            if pivot in basis:
                mask ^= basis[pivot]
            else:
                basis[pivot] = mask
                break
    return len(basis)


@lru_cache(maxsize=1)
def encoder():
    """Return a cached seven-qubit Clifford encoder accepting any q0 input.

    The initial state must be |psi> on q0 and |0> on q1..q6. H on q1..q3
    supplies three uniform binary variables. A reversible binary map sends
    columns 0..3 to (1111111, h0, h1, h2), so each basis term becomes
    x*1111111 + u0*h0 + u1*h1 + u2*h2 (mod 2). This realizes codewords()
    for arbitrary a|0> + b|1>, using gates rather than initialize().

    Row elimination synthesizes the CNOT/SWAP network; reversing it builds
    the forward map. The cached QuantumCircuit has no measurements and is
    shared between calls: compose it into a new circuit instead of mutating
    it. Its inverse decodes onto q0. See math-note.md §1.
    """
    columns = [LOGICAL_X, *ROWS]
    for q in range(7):
        if _rank([*columns, 1 << q]) > len(columns):
            columns.append(1 << q)
    matrix = np.array([[(col >> row) & 1 for col in columns]
                       for row in range(7)], dtype=np.uint8)
    ops = []
    for col in range(7):
        pivot = next(row for row in range(col, 7) if matrix[row, col])
        if pivot != col:
            matrix[[col, pivot]] = matrix[[pivot, col]]
            ops.append(("swap", col, pivot))
        for row in range(7):
            if row != col and matrix[row, col]:
                matrix[row] ^= matrix[col]
                ops.append(("cx", col, row))
    assert np.array_equal(matrix, np.eye(7, dtype=np.uint8))
    qc = QuantumCircuit(7, name="Steane_encode")
    qc.h([1, 2, 3])
    for gate, a, b in reversed(ops):
        getattr(qc, gate)(a, b)
    return qc


def preparation(state):
    """Return a one-qubit circuit mapping |0> to the named Pauli eigenstate.

    ``state`` is one of STATES: "0"/"1", "+"/"-", or "+i"/"-i".
    The last pair are (|0> ± i|1>)/sqrt(2), the Y eigenstates. Applying
    this circuit's inverse before Z measurement makes outcome 0 mean
    survival of that chosen input. Unknown labels raise ValueError.
    """
    qc = QuantumCircuit(1)
    if state == "1":
        qc.x(0)
    elif state in ("+", "-", "+i", "-i"):
        qc.h(0)
        if state == "-":
            qc.z(0)
        elif state == "+i":
            qc.s(0)
        elif state == "-i":
            qc.sdg(0)
    elif state != "0":
        raise ValueError(f"Unknown logical input {state!r}")
    return qc


def encoded_state(state):
    """Return the ideal Qiskit Statevector of one STATES input after encoding.

    The output has 128 complex amplitudes ordered by physical integer mask.
    This is an exact statevector reference; it performs no noisy simulator
    shots, syndrome measurement, or calibration.
    """
    qc = QuantumCircuit(7)
    qc.compose(preparation(state), [0], inplace=True)
    qc.compose(encoder(), inplace=True)
    return Statevector.from_instruction(qc)
