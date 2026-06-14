"""End-to-end: synthetic data must reproduce the published LEPR result."""
from isoform_dominance import _selftest


def test_reproduces_published_result(tmp_path):
    ok, msgs = _selftest.run(str(tmp_path))
    assert ok, "self-test failed:\n" + "\n".join(msgs)


def test_combined_pvalue(tmp_path):
    info = _selftest.generate(str(tmp_path))
    from isoform_dominance import extract, stats
    perdonor = {}
    for cohort, p in info.items():
        out = str(tmp_path / ("pd_%s.csv" % cohort))
        extract.run(_selftest.CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    res = stats.run(_selftest.CONFIG, "control", perdonor, str(tmp_path / "r"))
    cn, cgt, cp, _ = res["combined"]
    assert cn == 11 and cgt == 11
    assert abs(cp - 9.7656e-4) < 1e-4
