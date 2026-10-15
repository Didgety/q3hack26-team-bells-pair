"""Classical evaluation scores, fidelity conversions, and uncertainty.

No function here executes Qiskit or adds an Aer NoiseModel. Long records
sample the ideal logical channel already checked against short quantum
circuits. Calibration remains separate Qiskit data. Identical held-out
records score each decoder, and paired resampling preserves that comparison.
See docs/math-note.md §6 for metric and uncertainty assumptions.
"""
import numpy as np
from scipy.stats import norm
from .decoders import SyndromeRecord
from .logical_model import multiplier, symmetric_parameters


def mean_fidelity(fe):
    """Convert qubit entanglement fidelity Fe to mean state fidelity.

    ``fe`` is a scalar or array in [0,1]; the same-shaped result is
    F_mean=(2*Fe+1)/3 for a trace-preserving qubit channel with a fixed,
    input-independent recovery procedure. A uniform mean over the six
    Pauli eigenstates equals this state average. Fe is not generally
    one-state survival; that shortcut is specific to our ideal I/X channel.
    """
    return (2 * np.asarray(fe) + 1) / 3


def wilson(successes, shots, confidence=0.95):
    """Return [lower, upper] Wilson bounds for a binomial success probability.

    ``successes`` is K and ``shots`` is N independent binary observations,
    0<=K<=N. ``confidence`` is a probability in (0,1), usually 0.95.
    Returns two floats; N=0 returns [nan,nan]. This interval describes
    shot uncertainty for one measured state, not model uncertainty or
    independent calibration-batch variation. See the hierarchical method
    below for long decoder comparisons.
    """
    if shots == 0:
        return [float("nan"), float("nan")]
    z = norm.ppf((1+confidence)/2)
    p = successes / shots
    center = (p + z*z/(2*shots)) / (1+z*z/shots)
    half = z * np.sqrt(p*(1-p)/shots + z*z/(4*shots*shots)) / (1+z*z/shots)
    return [float(center-half), float(center+half)]


def sample_logical_records(theta, cycles, model, shots, seed):
    """Classical sampling of the validated channel, not Qiskit execution.

    Count k is a sufficient summary for these uniform Steane branch channels.
    Z-input failures estimate Fe only because the ideal channel is I/X.

    ``theta`` is radians, ``cycles`` is T, ``model`` is the generator's
    true F/R distribution, ``shots`` is N, and ``seed`` controls NumPy
    sampling. Draw k~Binomial(T,1-q0), then failure~Bernoulli((1-lambda)/2).
    Returns (SyndromeRecord(T,k), failure), with integer k and boolean
    failure arrays of shape (N,). A True failure means logical Z-input
    non-survival before final recovery. Only the record goes to Decoder;
    failures remain in scoring. Use seeds independent of calibration.
    There is no additional gate/readout noise and no long circuit execution.
    """
    rng = np.random.default_rng(seed)
    q0, _, _ = symmetric_parameters(theta)
    k = rng.binomial(cycles, np.clip(1-q0, 0, 1), size=shots)
    lam = multiplier(theta, cycles, k, model)
    failure = rng.random(shots) < (1-lam)/2
    return SyndromeRecord(cycles, k), failure


def score(failure, actions):
    """Return float 0/1 Z-input survival after final logical I/X actions.

    ``failure`` and ``actions`` are same-shaped boolean arrays (shots,).
    An X action toggles a logical Z-readout failure, so score=not(failure
    XOR action). Mean score estimates Fe only for this ideal I/X channel.
    Never use evaluation failures to choose the supplied decoder actions.
    """
    return (~np.logical_xor(failure, actions)).astype(float)


def paired_hierarchical_bootstrap(scores, repetitions=1500, seed=0):
    """Two-level paired bootstrap: independent batches, then held-out shots.

    scores has shape (calibration batches, evaluation shots, decoders).
    A joint multinomial over the observed score vectors implements the exact
    empirical shot bootstrap efficiently. Each replicate uses identical shot
    draws for all decoders. Batch resampling carries calibration variation.

    ``scores`` contains only binary 0/1 scores with positive batch/shot
    counts. ``repetitions`` is the number of bootstrap replicates;
    ``seed`` controls NumPy resampling. Returns a dict: mean and
    difference_mean have shape (decoders,); ci and difference_ci have
    shape (decoders,2), with lower/upper 2.5th/97.5th percentiles. Differences
    are paired against column 0, the standard decoder in the experiment.

    Intervals are approximate 95% empirical intervals, can be conservative,
    and cannot include rare model-selection failures absent from the
    observed batches. Compare with exact finite-binomial expectations.
    They cover sampling under the specified model, not hardware drift,
    unknown angle, model misspecification, or systematic calibration bias.
    """
    scores = np.asarray(scores, dtype=np.uint8)
    batches, shots, decoders = scores.shape
    rng = np.random.default_rng(seed)
    patterns, probs = [], []
    for batch in scores:
        unique, counts = np.unique(batch, axis=0, return_counts=True)
        patterns.append(unique)
        probs.append(counts / shots)
    boot = np.zeros((repetitions, decoders))
    for r in range(repetitions):
        for b in rng.integers(batches, size=batches):
            boot[r] += rng.multinomial(shots, probs[b]) @ patterns[b]
    boot /= batches * shots
    means = scores.mean(axis=(0, 1))
    return dict(mean=means, ci=np.quantile(boot, [0.025, 0.975], axis=0).T,
                difference_mean=means-means[0],
                difference_ci=np.quantile(boot-boot[:, [0]], [0.025, 0.975], axis=0).T)
