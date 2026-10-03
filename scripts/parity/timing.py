"""Time the local-annotation mode on whole files, as a user runs it.

    python scripts/parity/timing.py --gtf gencode.v50.annotation.gtf.gz \\
        --fasta gencode.v50.transcripts.fa.gz [--gene LEPR] [--repeat 3]

Runs ``isoform-dominance annotate --gtf`` for the gene, then for the gene's name in lower
case (the case-ignoring second read of the file), and ``identifiability --gtf
--transcripts-fasta`` on the config the first wrote, each ``--repeat`` times, through the
command line in a fresh interpreter, and prints the wall-clock seconds of each run.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time


def run(argv):
    t0 = time.perf_counter()
    r = subprocess.run([sys.executable, "-m", "isoform_dominance.cli", *argv],
                       capture_output=True, text=True)
    took = time.perf_counter() - t0
    if r.returncode not in (0, 2, 3):
        raise SystemExit("failed (%d): %s\n%s" % (r.returncode, " ".join(argv), r.stderr))
    return took


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--gtf", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gene", default="LEPR")
    p.add_argument("--repeat", type=int, default=3)
    a = p.parse_args(argv)
    with tempfile.TemporaryDirectory() as where:
        cfg = os.path.join(where, "cfg.json")
        steps = [
            ("annotate --gtf", ["annotate", "--gene", a.gene, "--gtf", a.gtf, "--out", cfg]),
            ("annotate --gtf, case ignored", ["annotate", "--gene", a.gene.lower(), "--gtf",
                                              a.gtf, "--out", os.path.join(where, "lc.json")]),
            ("identifiability --gtf --transcripts-fasta",
             ["identifiability", "--config", cfg, "--gtf", a.gtf, "--transcripts-fasta",
              a.fasta, "--json"]),
        ]
        for name, argv_ in steps:
            took = [run(argv_) for _ in range(a.repeat)]
            print("%-45s %s s (min %.2f)" % (name, ", ".join("%.2f" % t for t in took),
                                             min(took)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
