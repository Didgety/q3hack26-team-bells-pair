from dataclasses import fields
import numpy as np
from q3hack.logical_model import branches, multiplier, exact_fidelity, count_distribution
from q3hack.calibration import fit_model
from q3hack.decoders import SyndromeRecord, Decoder
from q3hack.metrics import paired_hierarchical_bootstrap, sample_logical_records


def test_decision_sensitive_reference_values():
    for theta, t, f, r, wrong in [
        (0.1, 128, 0.8230515273, 0.6111174360, 0.3845121731),
        (0.05, 929, 0.8296646494, 0.7438220766, 0.2548793741),
    ]:
        np.testing.assert_allclose([exact_fidelity(theta, t, "F", "F"),
                                    exact_fidelity(theta, t, "R", "R"),
                                    exact_fidelity(theta, t, "F", "R")],
                                   [f, r, wrong], atol=6e-11, rtol=0)


def test_note_polynomials_against_enumeration():
    for theta in (0, 0.05, 0.1, 0.2, 0.4):
        c, u = np.cos(theta), np.sin(theta)
        a0, b0 = c**7+7*c**3*u**4, -(7*c**4*u**3+u**7)
        an = c**6*u-4*c**4*u**3+3*c**2*u**5
        bn = 3*c**5*u**2-4*c**3*u**4+c*u**6
        br = branches(theta)
        np.testing.assert_allclose(br.a, [a0]+[an]*7, atol=1e-14)
        np.testing.assert_allclose(br.b, [b0]+[bn]*7, atol=1e-14)


def test_syndrome_law_and_single_interval_channel():
    for theta in [0, 0.05, 0.1, 0.2]:
        k, p = count_distribution(theta, 1)
        np.testing.assert_allclose(multiplier(theta, 1, k, "F"),
                                   multiplier(theta, 1, k, "R"), atol=1e-14)
        assert abs(p.sum()-1) < 1e-13
    assert abs(exact_fidelity(0, 929, "F", "F")-1) < 1e-14


def test_calibration_ties_and_information_boundary():
    assert fit_model(0, 0, 0.1).model == "R"
    assert fit_model(220, 1000, 0.1).model == "F"
    assert fit_model(110, 1000, 0.1).model == "R"
    assert [f.name for f in fields(SyndromeRecord)] == ["cycles", "nonzero_count"]
    assert [f.name for f in fields(Decoder)] == ["theta", "assumed_model"]
    rec, failure = sample_logical_records(0.1, 128, "F", 400, 5)
    other, _ = sample_logical_records(0.1, 128, "F", 400, 6)
    assert not np.array_equal(rec.nonzero_count, other.nonzero_count)


def test_bootstrap_preserves_pairing():
    rng = np.random.default_rng(90)
    base = rng.binomial(1, 0.6, (8, 64))
    scored = np.stack([base, base, 1-base], axis=2)
    result = paired_hierarchical_bootstrap(scored, repetitions=200, seed=3)
    np.testing.assert_array_equal(result["difference_ci"][1], [0, 0])
    assert np.all(result["difference_ci"][2] < 0)
