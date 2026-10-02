"""`extract` writes its rows in quant.sf path order, as 2.3.0 did.

2.4.0 sorted by donor name instead.  The two orders part where one donor name is a prefix
of another and the next character sorts before ``/``: ``D1-2/quant.sf`` < ``D1/quant.sf``
but ``D1`` < ``D1-2``.  `stats` bootstraps the fold interval by row index, so the order
changes the interval it reports for the same donors.

``tests/data/extract_order/perdonor_v230.csv`` was written by the 2.3.0 code
(commit 9c1f37d) from the quant.sf files beside it::

    git worktree add /tmp/v230 9c1f37d
    cd tests/data/extract_order
    PYTHONPATH=/tmp/v230/src python -m isoform_dominance.cli extract \\
        --config config.json --quantdir quant --samplemap samplemap.csv \\
        --cohort C --out perdonor_v230.csv
"""
import os
import shutil

from isoform_dominance import cli

DATA = os.path.join(os.path.dirname(__file__), "data", "extract_order")


def test_rows_are_in_quant_path_order_byte_for_byte_as_in_2_3_0(tmp_path):
    # copied, so that the sidecar and the CSV are written outside the source tree
    work = tmp_path / "case"
    shutil.copytree(DATA, work)
    out = work / "perdonor.csv"
    assert cli.main(["extract", "--config", str(work / "config.json"),
                     "--quantdir", str(work / "quant"),
                     "--samplemap", str(work / "samplemap.csv"),
                     "--cohort", "C", "--out", str(out)]) == 0
    with open(os.path.join(DATA, "perdonor_v230.csv"), "rb") as f:
        golden = f.read()
    assert out.read_bytes() == golden
