"""Hidden sign mixtures for deliberately programmed coherent data rotations.

Each interval applies U_xi(theta) = tensor_j exp(-i*xi*theta*X_j), with one
common sign across all seven data qubits. F (fixed) retains one fair sign
through a shot; R (refreshed) draws a new fair sign each interval. Both
models therefore have spatial correlation; only F retains temporal memory.
R is a model label, not a register. There is no bath qubit in this model.

These schedules specify explicit RX gates, not an Aer NoiseModel. Optional
simulator readout errors are configured separately in circuits.simulator.
Hidden schedules are erased before a decoder sees records. See math-note.md §2.
"""
from itertools import product
import numpy as np


def sign_mixture(model, cycles):
    """Return all hidden schedules and their probabilities for short circuits.

    Args:
        model: "F" or "R"; a generator setting, never a measured label.
        cycles: Nonnegative integer number of noise intervals T. R explicitly
            enumerates 2**T schedules and is limited to T <= 10 here.

    Returns:
        (schedules, weights): a list of length-T tuples containing ±1 and
        a same-length float probability array summing to one. T=0 returns
        one empty schedule with weight 1. No random sampling happens here;
        circuits._run_mixture samples independent shot memberships using
        its seed, then runs the corresponding Qiskit circuits.
    """
    if model not in ("F", "R") or cycles < 0:
        raise ValueError("Expected F/R and nonnegative cycles")
    if cycles == 0:
        return [()], np.array([1.0])
    if model == "F":
        return [(1,) * cycles, (-1,) * cycles], np.array([0.5, 0.5])
    if cycles > 10:
        raise ValueError("Use the validated reduced model for long storage")
    schedules = list(product((1, -1), repeat=cycles))
    return schedules, np.full(len(schedules), 1 / len(schedules))
