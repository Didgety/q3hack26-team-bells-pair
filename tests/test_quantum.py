import numpy as np
import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Pauli, Statevector
from q3hack.code import (H, STATES, codewords, encoded_state,
                         preparation, syndrome)
from q3hack.circuits import simulator, storage_circuit, run_storage, run_calibration
from q3hack.logical_model import branches, physical_branches, exact_fidelity
from q3hack.calibration import nonzero_probability


@pytest.mark.parametrize("state", STATES)
def test_encoder_stabilizers_and_logical_operators(state):
    v = codewords()
    psi = encoded_state(state)
    target = v @ Statevector.from_instruction(preparation(state)).data
    np.testing.assert_allclose(psi.data, target, atol=2e-14)
    for row in H:
        for axis in ("X", "Z"):
            label = "".join(axis if bit else "I" for bit in reversed(row))
            assert abs(psi.expectation_value(Pauli(label))-1) < 2e-14
    np.testing.assert_allclose(v.conj().T @ v, np.eye(2), atol=1e-14)
    np.testing.assert_allclose(v.conj().T @ Pauli("XXXXXXX").to_matrix() @ v,
                               [[0, 1], [1, 0]], atol=1e-14)
    np.testing.assert_allclose(v.conj().T @ Pauli("ZZZZZZZ").to_matrix() @ v,
                               [[1, 0], [0, -1]], atol=1e-14)


@pytest.mark.parametrize("state", STATES)
def test_zero_noise_and_all_single_x_dynamic_circuits(state):
    backend = simulator()
    circuits = [storage_circuit(0, (1, 1), state, injected_x=q) for q in range(7)]
    circuits += [storage_circuit(0, (1,)*t, state) for t in (0, 1, 2, 4)]
    compiled = transpile(circuits, backend, optimization_level=0)
    result = backend.run(compiled, shots=8, seed_simulator=36).result()
    for i in range(len(circuits)):
        counts = result.get_counts(i)
        assert len(counts) == 1
        value = int(next(iter(counts)).replace(" ", ""), 2)
        assert value == (syndrome(1 << i) if i < 7 else 0)


@pytest.mark.parametrize("theta", [0, 0.05, 0.1, 0.2])
@pytest.mark.parametrize("xi", [-1, 1])
def test_coherent_branches_against_qiskit(theta, xi):
    br = branches(theta, xi)
    physical, leakage = physical_branches(theta, xi)
    assert leakage < 2e-14
    np.testing.assert_allclose(physical, br.kraus, atol=2e-14)
    for s in range(8):
        np.testing.assert_allclose(br.kraus[s].conj().T @ br.kraus[s],
                                   br.q[s] * np.eye(2), atol=2e-14)
    other = branches(theta, -xi)
    np.testing.assert_allclose(br.q, other.q, atol=2e-14)
    np.testing.assert_allclose(br.g.conj(), other.g, atol=2e-14)


@pytest.mark.parametrize("theta", [0.1, 0.2])
def test_repeated_physical_projections(theta):
    """Independent 7-qubit propagation for selected 2/4-cycle histories."""
    from q3hack.code import RECOVERY
    v = codewords()
    sectors = np.array([syndrome(i) for i in range(128)])
    for signs in [(1, -1), (1, 1, 1, 1), (1, -1, -1, 1)]:
        histories = [(0,)*len(signs), tuple([1, 7, 2, 0][:len(signs)]),
                     tuple([0, 1, 0, 7][:len(signs)])]
        for history in histories:
            physical, logical = v.copy(), np.eye(2, dtype=complex)
            for xi, s in zip(signs, history):
                noise = QuantumCircuit(7)
                noise.rx(2*theta*xi, range(7))
                physical = np.column_stack([Statevector(physical[:, i]).evolve(noise).data for i in (0, 1)])
                physical *= (sectors == s)[:, None]
                physical = physical[np.arange(128) ^ RECOVERY[s]]
                logical = branches(theta, xi).kraus[s] @ logical
            np.testing.assert_allclose(physical, v @ logical, atol=2e-14)


def test_repeated_measured_circuit_and_calibration_rates():
    records = run_storage(0.2, 4, "F", "0", 4096, 1729)
    expected = exact_fidelity(0.2, 4, "F")
    assert abs(1-records.outcome.mean()-expected) < 0.035
    q = branches(0.2).q
    hist = np.bincount(records.syndromes.ravel(), minlength=8)/records.syndromes.size
    assert np.max(np.abs(hist-q)) < 0.025
    for model in ("F", "R"):
        syn, _, _ = run_calibration(0.1, model, 4096, 1977)
        assert abs(np.mean(syn != 0)-nonzero_probability(0.1, model)) < 0.03


@pytest.mark.parametrize("state", STATES)
def test_final_frame_matches_physical_logical_x(state):
    backend = simulator()
    circuits = [storage_circuit(0.15, (1, 1), state, final_x=x) for x in (False, True)]
    result = backend.run(transpile(circuits, backend, optimization_level=0),
                         shots=4096, memory=True, seed_simulator=91).result()
    rates = []
    for i in (0, 1):
        memory = [int(x.replace(" ", ""), 2) for x in result.get_memory(i)]
        rates.append(np.mean([(value >> 6) & 1 for value in memory]))
    # Independent simulator streams: compare probabilities, not bitwise samples.
    predicted = rates[0] if state in ("+", "-") else 1-rates[0]
    assert abs(predicted-rates[1]) < 0.04
