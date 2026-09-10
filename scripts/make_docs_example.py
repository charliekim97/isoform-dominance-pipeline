#!/usr/bin/env python3
"""Regenerate docs/example_output.{png,pdf,svg} and docs/example_output_stats.csv.

The committed example is the output of the bundled, download-free self-test, so it
is reproducible on any machine with the package installed and needs no data access:

    python scripts/make_docs_example.py

It must be re-run whenever the figure layout or the reported statistics change, so
that what the README shows is what the current code produces.

The output is byte-reproducible: ``SOURCE_DATE_EPOCH`` suppresses the creation
timestamp Matplotlib would otherwise embed in the SVG and PDF, and a fixed
``svg.hashsalt`` makes the generated element ids stable.  ``git status`` is therefore
a valid check that the committed artefacts match the current code.
"""
import os

os.environ.setdefault("SOURCE_DATE_EPOCH", "0")  # must precede the Matplotlib import

import shutil  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
mpl.rcParams["svg.hashsalt"] = "isoform-dominance-docs-example"

from isoform_dominance import _selftest, extract, stats  # noqa: E402

# Pin the family. The library's default stack prefers Arial, so the same code renders a
# different figure on a machine that has it (most macOS installs) than on one that does
# not (most Linux, and CI). DejaVu Sans ships with Matplotlib, so it is present
# everywhere and the committed example is the same file wherever it is rebuilt.
stats.FONT_STACK = ["DejaVu Sans"]

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(os.path.dirname(HERE), "docs")
STEM = "example_output"


def main():
    work = tempfile.mkdtemp(prefix="idp_docs_example_")
    try:
        info = _selftest.generate(work)
        perdonor = {}
        for cohort, paths in info.items():
            out = os.path.join(work, "perdonor_%s.csv" % cohort)
            extract.run(_selftest.CONFIG, paths["quantdir"], paths["samplemap"], cohort, out)
            perdonor[cohort] = out

        res = stats.run(_selftest.CONFIG, "control", perdonor, os.path.join(work, STEM))

        cn, cgt, cp, _ = res["combined"]
        en, eng, ep = _selftest.EXPECT["COMBINED"]
        if not (cn == en and cgt == eng and abs(cp - ep) < 1e-4):
            print("refusing to write: the reference result did not reproduce "
                  "(got n=%d, %d/%d, P=%.4g)" % (cn, cgt, cn, cp), file=sys.stderr)
            return 1

        os.makedirs(DOCS, exist_ok=True)
        for suffix in (".png", ".pdf", ".svg", "_stats.csv"):
            shutil.copyfile(os.path.join(work, STEM + suffix),
                            os.path.join(DOCS, STEM + suffix))
            print("wrote docs/%s%s" % (STEM, suffix))
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
