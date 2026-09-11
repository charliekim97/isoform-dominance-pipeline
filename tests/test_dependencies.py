"""The dependency floor in pyproject.toml is a claim about correctness; check it."""
import pathlib
import re

PYPROJECT = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"


def _floor(name):
    m = re.search(r'"%s>=([0-9.]+)"' % re.escape(name), PYPROJECT.read_text())
    assert m, "no %s>= floor in pyproject.toml" % name
    return tuple(int(x) for x in m.group(1).split(".")[:2])


def test_scipy_floor_excludes_the_anti_conservative_auto_rule():
    """Before SciPy 1.15, ``wilcoxon(method="auto")`` answered any zero difference with
    the normal approximation.  For the pairs in
    ``test_zeros_shrink_the_effective_n_and_raise_the_floor`` that is P = 0.0455 where
    the exact value is 0.125: an anti-conservative per-cohort p.  A floor that admits
    those versions admits that answer."""
    assert _floor("scipy") >= (1, 15)
