"""The package's own code calls no deprecated alias of its own.

The old names stay until 3.0 for callers outside the package (``counting_noise_floor``
for :func:`identifiability.unique_read_counting_se`); a call to one from inside it would
warn every user on every run.  A ``DeprecationWarning`` is made an error only when it is
attributed to a module of this package: the floor job's matplotlib 3.7 warns about
pyparsing 3.3's old names while drawing a figure, which is no call of ours.
"""
import json
import pathlib
import warnings

import pytest

from isoform_dominance import cli, identifiability, io

ROOT = pathlib.Path(__file__).resolve().parents[1]
SAVED = ROOT / "example" / "saved"

#: warnings attributed to this package's modules, and to nothing else
OURS = r"isoform_dominance(\.|$)"


def _ours_are_errors():
    warnings.filterwarnings("error", category=DeprecationWarning, module=OURS)


def test_the_filter_catches_this_package_and_nothing_else():
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        _ours_are_errors()
        for module in ("isoform_dominance", "isoform_dominance.identifiability"):
            with pytest.raises(DeprecationWarning):
                warnings.warn_explicit("ours", DeprecationWarning, "x.py", 1, module=module)
        for module in ("matplotlib._mathtext", "isoform_dominance_extra", "tests.test_x"):
            warnings.warn_explicit(module, DeprecationWarning, "x.py", 1, module=module)
    # warned, not raised
    assert [str(w.message) for w in seen] == [
        "matplotlib._mathtext", "isoform_dominance_extra", "tests.test_x"]
    with warnings.catch_warnings(), pytest.raises(DeprecationWarning):
        _ours_are_errors()
        _call_the_alias_as_package_code()


def _call_the_alias_as_package_code():
    """``counting_noise_floor`` called from code whose module is the package's: the alias
    warns with ``stacklevel=2``, so its warning is attributed to that caller."""
    exec(compile("identifiability.counting_noise_floor(10, 10)", "<package>", "exec"),
         {"__name__": "isoform_dominance.identifiability", "identifiability": identifiability})


def test_analyze_on_a_saved_input_calls_no_deprecated_alias():
    cfg = io.load_config(str(SAVED / "LEPR.gencode50.config.json"))
    saved = io.load_inputs(str(SAVED / "LEPR.gencode50.inputs.json"))
    with warnings.catch_warnings():
        _ours_are_errors()
        res = identifiability.analyze(
            cfg, sequences=saved["sequences"],
            background_sequences=(saved["background_sequences"]
                                  if saved["analysis"]["gene_background"] else None),
            background_gene_transcripts=False)
    assert res["unique_read_counting_se"]["defined"] is True


def test_identifiability_json_on_a_saved_input_calls_no_deprecated_alias(capsys):
    with warnings.catch_warnings():
        _ours_are_errors()
        code = cli.main(["identifiability", "--config", str(SAVED / "FLT1.gencode44.config.json"),
                         "--inputs", str(SAVED / "FLT1.gencode44.inputs.json"), "--json"])
    assert code == 0
    assert json.loads(capsys.readouterr().out)["verdict"] == "identifiable"


def test_selftest_calls_no_deprecated_alias(capsys):
    with warnings.catch_warnings():
        _ours_are_errors()
        code = cli.main(["selftest", "--json"])
    assert code == 0
    assert json.loads(capsys.readouterr().out)["ok"] is True
