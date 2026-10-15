"""Exact reduced logical channel for coherent X noise and ideal correction.

This module is classical algebra, apart from the independent Qiskit
statevector reference physical_branches(). It is not Aer shot execution.
It enumerates coherent physical amplitudes, validates their symmetry, and
reduces a T-cycle syndrome history to its count k of nonzero syndromes.
The resulting binomial sums support long storage without long circuit runs.

All formulas assume uniform programmed angles and ideal syndrome reports,
gates, and resets. Adding the optional Aer readout NoiseModel breaks these
assumptions; its short-circuit results must not be predicted by this model.
Derivations and citations: docs/math-note.md §§3–4, 7.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.stats import binom
from .code import CX, LOGICAL_X, RECOVERY, codewords, syndrome


@dataclass(frozen=True)
class Branches:
    """Eight corrected logical branches, indexed by syndrome integer s.

    ``kraus``: complex array (8, 2, 2), K_s=alpha_s I+beta_s X, before
    the optional irrelevant phase change. ``q``: float (8,), q_s from
    K_s†K_s=q_s I, so q is independent of logical input and hidden sign.
    ``a`` and ``b``: real (8,) after rephasing K_s to a_s I-i*xi*b_s X.
    ``g``: complex (8,), (a_s-i*xi*b_s)^2/q_s; the normalized branch's
    multiplier on |+><-| coherence. Its modulus is one for q_s>0;
    impossible branches use g_s=1 as a zero-weight placeholder.

    Arrays are shared by the cached branches() result; treat them as
    read-only even though the frozen dataclass does not freeze array data.
    """
    kraus: np.ndarray  # raw, without the optional phase adjustment
    q: np.ndarray
    g: np.ndarray
    a: np.ndarray
    b: np.ndarray


@lru_cache(maxsize=128)
def branches(theta, xi=1):
    """Sum all 128 physical X-support amplitudes into corrected branches.

    ``theta`` is the per-interval half-rotation angle in radians;
    ``xi`` is the hidden sign +1 or -1, not the F/R temporal-model label.
    Returns Branches for one interval, without random sampling or Aer noise.

    For support e of weight w, c_e=cos(theta)^(7-w)*(-i*xi*sin(theta))^w.
    Correction maps e to e XOR r_s. Membership in C_X versus its logical-X
    coset determines whether c_e contributes to alpha_s or beta_s.
    Amplitudes in each class MUST be added before taking magnitudes:
    summing |c_e|^2 instead would erase interference and change the model.
    See math-note.md §3 for explicit a_s/b_s polynomials and normalization.
    """
    if xi not in (-1, 1):
        raise ValueError("xi must be a sign")
    alpha, beta = np.zeros(8, complex), np.zeros(8, complex)
    for e in range(128):
        s, w = syndrome(e), e.bit_count()
        amplitude = np.cos(theta)**(7-w) * (-1j * xi * np.sin(theta))**w
        corrected = e ^ RECOVERY[s]
        # Different supports can have the same syndrome and logical action.
        # They interfere coherently within alpha or beta, not as Pauli draws.
        if corrected in CX:
            alpha[s] += amplitude
        else:
            assert corrected ^ LOGICAL_X in CX
            beta[s] += amplitude
    kraus = np.array([[[a, b], [b, a]] for a, b in zip(alpha, beta)])
    # A syndrome-dependent global phase leaves the physical channel intact
    # and exposes the real a,b representation used by the history formulas.
    phase = np.array([1, *([1j * xi] * 7)])
    a, b = (alpha * phase).real, -(beta * phase).imag / xi
    q = np.abs(alpha)**2 + np.abs(beta)**2
    g = np.ones(8, complex)
    np.divide((a - 1j * xi * b)**2, q, out=g, where=q > 0)
    np.testing.assert_allclose(q.sum(), 1, atol=2e-13)
    return Branches(kraus, q, g, a, b)


def symmetric_parameters(theta):
    """Return (q0, g0, gn) after checking all seven nonzero branches agree.

    ``theta`` is radians. q0 is the no-syndrome probability; g0 and gn are
    complex normalized coherence multipliers for xi=+1 and zero/nonzero
    syndromes. The total nonzero probability is 1-q0 = 7*q_n, not q_n.
    This equality permits count compression; an assertion stops the
    calculation if a changed branch model breaks it.
    """
    br = branches(float(theta))
    # This compression must fail rather than silently hide broken symmetry.
    np.testing.assert_allclose(br.q[1:], br.q[1], atol=1e-14)
    np.testing.assert_allclose(br.g[1:], br.g[1], atol=1e-13)
    return br.q[0], br.g[0], br.g[1]


def multiplier(theta, cycles, nonzero_count, model):
    """Return the conditional logical-channel multiplier lambda_h(k).

    Args:
        theta: Known programmed half-angle in radians.
        cycles: Number T of ideal correction intervals.
        nonzero_count: Integer scalar/array k, with 0 <= k <= T. Counts
            nonzero syndrome cycles, not physical faults or syndrome bits.
        model: Assumed/true temporal distribution "F" or "R" for this
            calculation; never an individual realized sign sequence.

    Returns:
        Real scalar/array with the shape of k and values in [-1,1]. The
        averaged conditional channel is ((1+lambda)/2)*rho plus
        ((1-lambda)/2)*X rho X. F averages one common sign after multiplying
        branch coherences; R averages signs independently before multiplying:
        lambda_F=Re(g0**(T-k)*gn**k),
        lambda_R=Re(g0)**(T-k)*Re(gn)**k. This ordering is the entire temporal
        distinction, even though both models have the same p(k). See §4.
    """
    _, g0, gn = symmetric_parameters(theta)
    k = np.asarray(nonzero_count, dtype=int)
    if np.any(k < 0) or np.any(k > cycles):
        raise ValueError("Nonzero count outside storage schedule")
    if model == "F":
        result = np.cos((cycles-k) * np.angle(g0) + k * np.angle(gn))
    elif model == "R":
        result = g0.real**(cycles-k) * gn.real**k
    else:
        raise ValueError("Expected F or R")
    return np.clip(result, -1, 1)


def count_distribution(theta, cycles):
    """Return all (k, p_T(k)) for ideal T-cycle storage, under either F or R.

    ``theta`` is radians and ``cycles`` is nonnegative integer T. The two
    arrays have shape (T+1,): k=0..T and Binomial(T,1-q0) probabilities.
    Every count is included and normalization is checked; no sampled shots
    or discarded tails enter the exact fidelity calculation.
    """
    q0, _, _ = symmetric_parameters(theta)
    k = np.arange(cycles + 1)
    p = binom.pmf(k, cycles, np.clip(1-q0, 0, 1))
    np.testing.assert_allclose(p.sum(), 1, atol=2e-12)
    return k, p


def exact_fidelity(theta, cycles, true_model, assumed_model=None):
    """Return exact entanglement fidelity Fe, a float in [0,1].

    ``theta`` is radians, ``cycles`` is T, and ``true_model`` is F/R for
    generating the ideal channel. ``assumed_model`` is None for final I,
    or F/R for that model's optimal final I/X rule. The distribution oracle
    uses assumed_model=true_model; it still has no realized hidden signs.

    Sum p_T(k)*(1+eta(k)*lambda_true(k))/2 over every count. eta=+1 means
    I; eta=-1 means logical X. This is an exact classical model value,
    not a measured Qiskit estimate and not a device prediction. Fe equals
    Z-input survival only for the ideal I/X channel derived here (§§4,6).
    """
    k, p = count_distribution(theta, cycles)
    lam = multiplier(theta, cycles, k, true_model)
    eta = 1 if assumed_model is None else np.where(
        multiplier(theta, cycles, k, assumed_model) < 0, -1, 1)
    return float(np.dot(p, (1 + eta * lam) / 2))


def bare_fidelity(theta, cycles, model):
    """Return exact unencoded entanglement fidelity for model "F" or "R".

    ``theta`` is radians and ``cycles`` is interval count T. The coherence
    is cos(2*T*theta) for a fixed sign and cos(2*theta)**T for refreshed
    signs; Fe=(1+coherence)/2. No shots, QEC, final recovery, or Aer
    NoiseModel enter this classical baseline (§7). Valid model labels are
    required from the caller. F can revive at longer times.
    """
    coherence = np.cos(2 * cycles * theta) if model == "F" else np.cos(2*theta)**cycles
    return float((1 + coherence) / 2)


def physical_branches(theta, xi=1):
    """Qiskit statevector reference: rotate, project, correct, decode.

    No physical-error amplitude enumeration is used in this calculation.

    ``theta`` is radians and ``xi`` is ±1. Returns (matrices, leakage),
    where matrices is complex shape (8,2,2) holding V† R_s Pi_s U_xi V
    and leakage is the largest residual norm ||R_s Pi_s U_xi V-V K_s||.
    This exact Qiskit Statevector calculation cross-checks branches();
    it uses ideal operators, mathematical syndrome projectors, and no
    Aer NoiseModel or sampled dynamic-circuit shots. Separate run_storage
    tests validate actual measurement/reset/conditional correction.
    """
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    v = codewords()
    qc = QuantumCircuit(7)
    qc.rx(2 * xi * theta, range(7))
    rotated = np.column_stack([Statevector(v[:, i]).evolve(qc).data for i in (0, 1)])
    sectors = np.array([syndrome(i) for i in range(128)])
    matrices, leakage = [], []
    for s in range(8):
        projected = rotated * (sectors == s)[:, None]
        corrected = projected[np.arange(128) ^ RECOVERY[s]]
        logical = v.conj().T @ corrected
        matrices.append(logical)
        leakage.append(np.linalg.norm(corrected - v @ logical))
    return np.array(matrices), max(leakage)
