"""What the README, the paper and the docstrings say a method does, held to the code.

Stouffer combines each cohort's signed-rank p-value as SciPy's ``method="auto"`` gives it:
exact at n <= 13, and up to n = 50 with neither a zero difference nor a tie, and the normal
approximation otherwise.  Through 2.11.0 the README, the paper and the docstrings of
``stats`` called its inputs the per-cohort *exact* tests (2.4.1's changelog corrected the
same sentence in 2.3.0's entry).  These tests hold the sentence that replaced it to the
dispatch at its edges and to the self-test, and keep the old one out.
"""
import pathlib
import re

import numpy as np
import pytest

from isoform_dominance import _selftest, stats

ROOT = pathlib.Path(__file__).resolve().parents[1]


def flat(text):
    """``text`` with ``*`` taken out and every run of white space made one space."""
    return re.sub(r"\s+", " ", text.replace("*", ""))


def readme():
    return flat((ROOT / "README.md").read_text(encoding="utf-8"))


#: the README's sentence on what Stouffer combines
STOUFFER_SENTENCE = ("a cohort's p-value is exact at n ≤ 13, and up to n = 50 with neither a "
                     "zero difference nor a tie (SciPy's exact null, or the exhaustive "
                     "permutation test), while the stratified statistic always rests on a "
                     "normal approximation. Otherwise a cohort's p is SciPy's normal "
                     "approximation, and Stouffer takes it as it is")


def _clean(n):
    """``n`` differences with no zero and no tie among their absolute values."""
    return np.arange(1, n + 1, dtype=float) * np.where(np.arange(n) % 3 == 0, -1.0, 1.0)


def _tied(n):
    d = _clean(n)
    d[1] = -d[2]
    return d


def _zero(n):
    d = _clean(n)
    d[0] = 0.0
    return d


@pytest.mark.parametrize("d, method", [
    (_clean(5), "exact"), (_tied(13), "permutation"), (_zero(13), "permutation"),
    (_tied(14), "asymptotic"), (_zero(14), "asymptotic"), (_clean(50), "exact"),
    (_clean(51), "asymptotic"),
], ids=["5 clean", "13 tied", "13 zero", "14 tied", "14 zero", "50 clean", "51 clean"])
def test_the_readme_gives_the_dispatch_stouffer_takes(d, method):
    assert STOUFFER_SENTENCE in readme()
    got = stats.paired_stat_detail(np.abs(d) + 10.0 + d, np.abs(d) + 10.0, n_boot=0)
    assert got["wilcoxon_method"] == method
    assert stats._wilcoxon_auto_method(d) == method


def test_both_cohorts_of_the_self_test_are_exact():
    assert "Both cohorts of the self-test are exact." in readme()
    cohorts = {c["name"]: c for c in _selftest.result()["checks"]
               if c["name"] in ("GSE228458", "GSE137619")}
    assert {k: v["wilcoxon_method"] for k, v in cohorts.items()} == {
        "GSE228458": "exact", "GSE137619": "exact"}


#: what was said of Stouffer's inputs through 2.11.0
OLD = ("per-cohort exact tests", "exact p-value is trustworthy",
       "keeps each cohort's test exact", "keeps each stratum's test exact")


@pytest.mark.parametrize("path", ["README.md", "docs/index.md", "paper.md",
                                  "src/isoform_dominance/stats.py"])
def test_no_text_calls_stouffers_inputs_exact(path):
    text = flat((ROOT / path).read_text(encoding="utf-8"))
    for old in OLD:
        assert old not in text, old
