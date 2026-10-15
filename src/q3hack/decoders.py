"""Restricted final logical I/X decisions after fixed per-cycle correction.

The storage circuit always uses the same physical single-X correction
table. These decoders differ only in the distribution they assume when
choosing the final Pauli. They receive a syndrome summary and stored model
parameters; outcomes, input-state labels, and realized signs are excluded.
The ideal I/X-channel optimality proof is in docs/math-note.md §4. It does
not establish noisy-circuit optimality when Aer report errors are added.
"""
from dataclasses import dataclass
import numpy as np
from .calibration import CalibrationFit
from .logical_model import multiplier


@dataclass(frozen=True)
class SyndromeRecord:
    """Allowed observations: interval count T and per-shot nonzero counts k.

    ``cycles`` is nonnegative integer T. ``nonzero_count`` is an integer
    array (shots,) with values in [0,T]. It counts nonzero syndrome cycles,
    not syndrome Hamming weights or physical-error counts. Seven nonzero
    branches are symmetric in this ideal model, making k sufficient.
    ``SyndromeRecord`` contains neither outcome nor input or hidden model.
    """
    cycles: int
    nonzero_count: np.ndarray


@dataclass(frozen=True)
class Decoder:
    """Known-angle final recovery rule with an optional assumed F/R model.

    ``theta`` is the known programmed half-angle in radians. ``assumed_model``
    None means standard final I; "R" means the independence assumption;
    a separately selected F/R gives the calibrated rule. Supplying the true
    distribution gives the comparison oracle, but no per-shot latent sign.
    The last two roles use the same class and action interface.
    """
    theta: float
    assumed_model: str | None = None

    def actions(self, record: SyndromeRecord):
        """Return bool array (shots,), False=I and True=final logical X.

        For conditional channel multiplier lambda, choose X iff lambda<0;
        this maximizes (1+eta*lambda)/2 with eta=+1 for I and -1 for X.
        Ties choose I. Only record.cycles and record.nonzero_count enter.
        Applying these actions to measured outputs belongs to the scorer,
        so the rule cannot choose a recovery after seeing success/failure.
        """
        if self.assumed_model is None:
            return np.zeros_like(record.nonzero_count, dtype=bool)
        return multiplier(self.theta, record.cycles, record.nonzero_count,
                          self.assumed_model) < 0

    @classmethod
    def calibrated(cls, fit: CalibrationFit):
        """Create the final rule from a separate CalibrationFit's theta/model.

        Only these two inferred/known parameters are retained. Probe counts
        justify model selection but are not per-shot evaluation information.
        """
        return cls(fit.theta, fit.model)
