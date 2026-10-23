"""An estimand past the linearisation limit counts as not resolved.

`paper.md` and the 2.3.0 changelog say a figure past ``LINEARISATION_LIMIT`` reads as "not
resolvable at this design".  Through 2.4.0 ``effect_resolvable`` compared the figure with
``--min-log2fc`` and ignored ``beyond_linear``: a contrast at relative SE 0.32 gave a
figure of 0.91, under 1.5, and the run said "resolved" and exited 0.

The EFFECT SIZE line says why each estimand is not resolved.  Through 2.8.0 it put every one
flagged ``beyond_linear`` "beyond the linearisation limit, relative SE > 0.3" -- among them
those with no figure at all, whose relative SE is infinite: NTRK2's two class totals and
contrast at GENCODE 50, none of them estimable, as the walkthrough showed.
"""
import json
import pathlib
import random

from isoform_dominance import cli, identifiability
from test_strict_json import _not_estimable_case

SAVED = pathlib.Path(__file__).resolve().parents[1] / "example" / "saved"


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_BASE = _seq(1500, 1)
SEQS = {"T1": _BASE[:1200] + _seq(300, 2), "T2": _BASE[:1200] + _seq(300, 3),
        "T3": _BASE[300:] + _seq(200, 4), "T4": _seq(200, 5) + _BASE[300:]}
CONFIG = {"groups": {"A": ["T1", "T2"], "B": ["T3", "T4"]}, "primary_comparison": ["A", "B"]}


def _run(tmp_path, capsys, *extra):
    cfg, seqs = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CONFIG))
    seqs.write_text(json.dumps(SEQS))
    argv = ["identifiability", "--config", str(cfg), "--sequences", str(seqs),
            "--no-gene-background"] + list(extra)
    code = cli.main(argv + ["--json"])
    res = json.loads(capsys.readouterr().out)
    cli.main(argv)
    return code, res, capsys.readouterr().out


def test_a_contrast_past_the_limit_is_not_resolved_whatever_its_figure(tmp_path, capsys):
    code, res, out = _run(tmp_path, capsys, "--tpm", "1.0", "--min-log2fc", "1.5")
    c = res["contrast"]
    # the case: past the limit, with a figure under the requested effect
    assert c["beyond_linear"] and c["gls_relative_se"] > identifiability.LINEARISATION_LIMIT
    assert c["min_resolvable_log2fc"] < 1.5
    assert res["effect_resolvable"] is False
    assert code == cli.EXIT_EFFECT_NOT_RESOLVED
    line = [ln for ln in out.splitlines() if "EFFECT SIZE" in ln]
    assert line == ["  EFFECT SIZE: |log2FC| 1.500 NOT resolved at this design (the contrast "
                    "beyond the linearisation limit, relative SE > 0.3)"]
    # the figure and the flag are still reported as they were
    assert c["min_resolvable_log2fc"] == identifiability.min_resolvable_log2fc(
        c["gls_relative_se"], 1)


def test_inside_the_limit_nothing_changes(tmp_path, capsys):
    code, res, out = _run(tmp_path, capsys, "--tpm", "10", "--min-log2fc", "1.5")
    assert not any(res["groups"][g]["beyond_linear"] for g in "AB")
    assert not res["contrast"]["beyond_linear"]
    assert res["effect_resolvable"] is True and code == cli.EXIT_OK
    assert "linearisation" not in [ln for ln in out.splitlines() if "EFFECT SIZE" in ln][0]


def test_without_an_effect_size_the_exit_status_is_still_zero(tmp_path, capsys):
    code, res, _ = _run(tmp_path, capsys, "--tpm", "1.0")
    assert res["contrast"]["beyond_linear"]
    assert res["effect_resolvable"] is None and code == cli.EXIT_OK


def _effect_line(capsys):
    out, err = capsys.readouterr()
    return [ln for ln in out.splitlines() if "EFFECT SIZE" in ln], out, err


def test_the_starred_figure_has_its_footnote(tmp_path, capsys):
    cfg, seqs = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CONFIG))
    seqs.write_text(json.dumps(SEQS))
    cli.main(["identifiability", "--config", str(cfg), "--sequences", str(seqs),
              "--no-gene-background", "--tpm", "1.0", "--min-log2fc", "1.5"])
    _, out, err = _effect_line(capsys)
    assert [ln for ln in out.splitlines() if ln.startswith("  contrast ")
            and "*;" in ln], out
    assert "* first-order figure past the linearisation limit" in err


def test_no_figure_is_not_past_the_limit(capsys):
    # the walkthrough's second run: no estimand of NTRK2 at GENCODE 50 is estimable, and
    # 2.8.0 said "(iso_477aa, iso_822aa, the contrast beyond the linearisation limit,
    # relative SE > 0.3)", with the footnote on starred figures where none was starred
    code = cli.main(["identifiability", "--config", str(SAVED / "NTRK2.gencode50.config.json"),
                     "--inputs", str(SAVED / "NTRK2.gencode50.inputs.json"),
                     "--min-log2fc", "0.5"])
    assert code == cli.EXIT_EFFECT_NOT_RESOLVED
    line, out, err = _effect_line(capsys)
    assert line == ["  EFFECT SIZE: |log2FC| 0.500 NOT resolved at this design (iso_477aa, "
                    "iso_822aa, the contrast not estimable)"]
    assert "*;" not in out and "first-order figure" not in err


def test_an_estimable_contrast_with_no_figure_says_so(tmp_path, capsys):
    # the class totals are not estimable and the contrast is, with no finite figure
    argv = [a for a in _not_estimable_case(tmp_path) if a != "--json"]
    assert cli.main(argv) == cli.EXIT_EFFECT_NOT_RESOLVED
    line, out, err = _effect_line(capsys)
    assert line == ["  EFFECT SIZE: |log2FC| 0.500 NOT resolved at this design (A, B not "
                    "estimable; the contrast estimable, but with no finite figure)"]
    assert "first-order figure" not in err
