"""Class effective-length ratio and where the class-distinguishing windows sit.

In simulation Salmon split ambiguous mass between the classes by effective length under
positional coverage skew, so the class effective-length ratio tracked the size of the
error and, with the direction of the skew, its sign.  The skew itself cannot be seen from
sequence; the report gives the two quantities that can, and states the direction only
inside the band where it was measured.
"""
import json
import math
import random

import pytest

from isoform_dominance import cli, identifiability as I


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


# --------------------------------------------------------------------------- #
# effective-length ratio
# --------------------------------------------------------------------------- #
def test_the_ratio_reproduces_the_simulation_harness():
    """ALDOA in class_lengths.json: class means 588 and 1776.8 nt at frag_mean 200 give
    log2 ratio -2.0200722823562707, the value the -0.547 correlation was computed on.
    These lengths have different harmonic means, so a harmonic mean fails here."""
    a = [512, 600, 652]
    b = [1487] + [1809] * 9
    mean_a, mean_b, ratio = I.class_efflen_ratio(a, b, frag_mean=200.0)
    assert (mean_a, mean_b) == pytest.approx((389.0, 1577.8))
    assert ratio == pytest.approx(-2.0200722823562707, abs=1e-12)


def test_transcript_count_does_not_enter_the_ratio():
    """One transcript against five of the same length is a ratio of 0, not log2(1/5).
    A class total picks up class size: it read -0.344 at uniform coverage, where the
    answer has to be null."""
    _, _, ratio = I.class_efflen_ratio([1199], [1199] * 5, frag_mean=200.0)
    assert ratio == 0.0


def test_the_ratio_uses_the_design_fragment_length():
    cfg = {"groups": {"A": ["A1"], "B": ["B1"]}, "primary_comparison": ["A", "B"]}
    seqs = {"A1": _seq(600, 1), "B1": _seq(1800, 2)}
    at = {fm: I.analyze(cfg, sequences=seqs, frag_mean=fm)["contrast"]["log2_efflen_ratio"]
          for fm in (100.0, 200.0)}
    assert at[100.0] == pytest.approx(math.log2(501 / 1701))
    assert at[200.0] == pytest.approx(math.log2(401 / 1601))


# --------------------------------------------------------------------------- #
# positions of the distinguishing windows
# --------------------------------------------------------------------------- #
_SHARED = _seq(1500, 21)


def test_positions_say_which_end_distinguishes_each_class():
    cfg = {"groups": {"A": ["A1"], "B": ["B1"]}, "primary_comparison": ["A", "B"]}
    seqs = {"A1": _SHARED + _seq(300, 22),          # A's own sequence at its 3' end
            "B1": _seq(300, 23) + _SHARED}          # B's at its 5' end
    pos = I.analyze(cfg, sequences=seqs)["contrast"]["distinguishing_window_position"]
    assert pos["A"]["median"] > 0.8 and pos["B"]["median"] < 0.2
    for p in pos.values():
        assert 0.0 <= p["q1"] <= p["median"] <= p["q3"] <= 1.0
        assert p["n"] > 0


def test_positions_are_in_each_transcripts_own_coordinates():
    """A short and a long transcript, each distinguished at its own 3' end.  In its own
    coordinates every distinguishing window of class A sits past 0.7; measured against
    the longest transcript, A1's would sit near 0.2 and pull the lower quartile down."""
    s1, s2 = _seq(500, 31), _seq(3000, 32)
    cfg = {"groups": {"A": ["A1", "A2"], "B": ["B1"]}, "primary_comparison": ["A", "B"]}
    seqs = {"A1": s1 + _seq(200, 33), "A2": s2 + _seq(200, 34), "B1": s1 + s2}
    pos = I.analyze(cfg, sequences=seqs)["contrast"]["distinguishing_window_position"]
    assert pos["A"]["q1"] > 0.7


def test_a_class_with_no_distinguishing_window_reports_none():
    cfg = {"groups": {"sub": ["S1"], "sup": ["P1"]}, "primary_comparison": ["sub", "sup"]}
    seqs = {"S1": _SHARED, "P1": _SHARED + _seq(400, 41)}
    pos = I.analyze(cfg, sequences=seqs)["contrast"]["distinguishing_window_position"]
    assert pos["sub"] == {"n": 0, "q1": None, "median": None, "q3": None}
    assert pos["sup"]["n"] > 0


# --------------------------------------------------------------------------- #
# the CLI: ratio and positions always, direction only inside the measured band
# --------------------------------------------------------------------------- #
EVIDENCE = "in a 49-gene simulation (Salmon, one quantifier, monotone positional skew)"


def _run(tmp_path, capsys, len_a, len_b, primary=("A", "B"), json_out=False):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"groups": {"A": ["A1"], "B": ["B1"]},
                               "primary_comparison": list(primary),
                               "ensembl_release": 116}))
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps({"A1": _seq(len_a, 51), "B1": _seq(len_b, 52)}))
    argv = ["identify", "--config", str(cfg), "--sequences", str(sq)]
    cli.main(argv + (["--json"] if json_out else []))
    return capsys.readouterr()


def test_below_the_band_the_ratio_and_positions_are_printed_without_a_direction(
        tmp_path, capsys):
    out = _run(tmp_path, capsys, 1199, 1422)          # efflen 1000 vs 1223: log2 0.290
    assert "effective length" in out.out
    assert "log2 ratio -0.29" in out.out
    assert "distinguishing windows" in out.out
    assert "5'-skewed" not in out.err and EVIDENCE not in out.err


def test_at_the_band_edge_the_direction_is_printed_with_its_evidence(tmp_path, capsys):
    out = _run(tmp_path, capsys, 1199, 1439)          # efflen 1000 vs 1240: log2 0.310
    assert "log2 ratio -0.31" in out.out
    assert EVIDENCE in out.err
    assert "shorter class (A)" in out.err and "longer (B)" in out.err
    assert "35-36 of the 39 genes with a ratio above 1.23x" in out.err
    assert "30-32" in out.err
    # the rule assumes ambiguity where the skew piles reads up; LEPR breaks it, so say so
    assert ("When the shorter class is itself distinguished at that end -- see its "
            "distinguishing-window position above -- the direction can reverse; the "
            "interaction was not measured.") in out.err


def test_the_direction_names_the_shorter_class_not_the_first_one(tmp_path, capsys):
    out = _run(tmp_path, capsys, 3000, 1000, primary=("A", "B"))
    assert "shorter class (B)" in out.err and "longer (A)" in out.err


def test_the_json_report_carries_ratio_positions_and_band(tmp_path, capsys):
    rep = json.loads(_run(tmp_path, capsys, 1199, 1439, json_out=True).out)["contrast"]
    assert rep["log2_efflen_ratio"] == pytest.approx(math.log2(1000 / 1240))
    assert rep["class_mean_efflen"] == {"A": 1000.0, "B": 1240.0}
    assert set(rep["distinguishing_window_position"]) == {"A", "B"}
    assert rep["efflen_direction_in_band"] is True
    rep = json.loads(_run(tmp_path, capsys, 1199, 1422, json_out=True).out)["contrast"]
    assert rep["efflen_direction_in_band"] is False
