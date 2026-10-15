"""Select F versus R using separate known-angle temporal probes.

circuits.run_calibration executes the probes in Qiskit, with programmed
coherent rotations and no additional Aer NoiseModel. This module only fits
their observed nonzero count. It receives no injection signs, true model,
logical input labels, or held-out evaluation outcomes. The programmed theta
is assumed known, not estimated; noisy-report/unknown-angle fitting is not
implemented. See docs/math-note.md §5 for the binomial likelihood.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import xlogy
from .logical_model import branches


@dataclass(frozen=True)
class CalibrationFit:
    """Observed calibration summary and selected decoder distribution.

    ``model`` is selected F/R, not ground truth. ``theta`` is the known
    half-angle in radians. ``shots`` is N independent two-interval probes;
    ``nonzero`` is K probes with syndrome != 0. ``log_odds`` is the natural
    logarithm of posterior F/R odds, including the supplied prior. Positive
    odds select F; zero/negative odds select R. None of these fields stores
    an evaluation outcome or the hidden sign of any probe.
    """
    model: str
    theta: float  # Known programmed parameter, supplied to all model decoders.
    shots: int
    nonzero: int
    log_odds: float


def nonzero_probability(theta, model):
    """Return the exact probability of a nonzero delayed-probe syndrome.

    ``theta`` is the half-angle per interval in radians and ``model`` is
    F/R. Two equal signs yield effective half-angle 2*theta, hence
    Q=1-q0(2*theta). F always has equal signs, giving Q; R cancels in half
    its assignments, giving Q/2. These ideal-operation predictions are
    dimensionless probabilities, not fitted data or Aer-noise rates.
    """
    q = float(np.clip(1 - branches(2 * theta).q[0], 0, 1))
    if model not in ("F", "R"):
        raise ValueError("Expected F or R")
    return q if model == "F" else q / 2


def fit_model(nonzero, shots, theta, prior_fixed=0.5):
    """Fit a temporal-model label from observed K, N and known theta.

    ``nonzero`` is integer K, ``shots`` is nonnegative integer N with
    0<=K<=N, ``theta`` is radians, and ``prior_fixed`` is prior P(F) in
    (0,1). Returns a CalibrationFit. For each candidate with probability p,
    log L = K*log(p)+(N-K)*log(1-p); the common binomial coefficient cancels.
    xlogy handles zero-count times log(0) without introducing spurious NaNs.

    Positive posterior log odds select F; ties select R. With N=0 or
    indistinguishable candidate rates, the prior alone decides (default R
    for an equal prior). This interface deliberately excludes the true
    injected model, sign sequence, and evaluation data.
    """
    if not (0 <= nonzero <= shots) or not 0 < prior_fixed < 1:
        raise ValueError("Invalid counts or prior")
    p_f, p_r = (nonzero_probability(theta, m) for m in ("F", "R"))
    if shots == 0 or p_f == p_r:
        odds = np.log(prior_fixed / (1-prior_fixed))
    else:
        ll_f = xlogy(nonzero, p_f) + xlogy(shots-nonzero, 1-p_f)
        ll_r = xlogy(nonzero, p_r) + xlogy(shots-nonzero, 1-p_r)
        odds = float(ll_f - ll_r + np.log(prior_fixed / (1-prior_fixed)))
    return CalibrationFit("F" if odds > 0 else "R", float(theta), shots, nonzero, float(odds))
