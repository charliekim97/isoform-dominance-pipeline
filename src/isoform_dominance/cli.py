"""Unified command-line interface: isoform-dominance <subcommand>."""
import argparse, json, sys
from urllib.error import URLError, HTTPError
from . import __version__, io, extract, stats, contamination, annotate, identifiability


def _kv(items):
    out = {}
    for s in (items or []):
        if "=" not in s:
            raise SystemExit("argument error: expected NAME=path, got %r" % s)
        k, v = s.split("=", 1)
        if not k or not v:
            raise SystemExit("argument error: expected NAME=path, got %r" % s)
        out[k] = v
    return out


def _net_fail(e):
    print("Ensembl request failed (%s). Check your network connection, the gene symbol, "
          "and species, or supply sequences offline." % e, file=sys.stderr)
    return 1


def cmd_annotate(a):
    try:
        cfg = annotate.run(a.gene, a.out, species=a.species)
    except (URLError, HTTPError) as e:
        return _net_fail(e)
    print("Proposed groups for %s -> %s" % (a.gene, a.out))
    for g, ids in cfg["groups"].items():
        print("  %s: %d transcripts" % (g, len(ids)))
    print("  primary_comparison:", cfg["primary_comparison"])
    print("  REVIEW _proposed/_clusters and rename groups before use.")


def cmd_identifiability(a):
    cfg = io.load_config(a.config)
    seqs = json.load(open(a.sequences)) if a.sequences else None
    try:
        res = identifiability.analyze(cfg, k=a.k, sequences=seqs)
    except (URLError, HTTPError) as e:
        return _net_fail(e)
    except ValueError as e:
        print("config error: %s" % e, file=sys.stderr)
        return 1
    print("Identifiability (k=%d):" % res["k"])
    for g, r in res["groups"].items():
        flag = "OK" if r["distinguishable"] else "NOT DISTINGUISHABLE"
        print("  [%s] %s: %d unique k-mers (%d transcripts)"
              % (flag, g, r["n_unique_kmers"], r["n_transcripts"]))
    print("  primary_comparison distinguishable by short reads:", res["primary_distinguishable"])
    if not res["primary_distinguishable"]:
        print("  WARNING: a primary group has no unique k-mers - short-read quantification "
              "cannot resolve it. Reconsider the grouping or use long reads.", file=sys.stderr)
        return 2
    return 0


def cmd_extract(a):
    n = extract.run(io.load_config(a.config), a.quantdir, a.samplemap, a.cohort, a.out)
    print("wrote %s (n=%d donors)" % (a.out, n))


def cmd_stats(a):
    res = stats.run(io.load_config(a.config), a.condition, _kv(a.perdonor), a.out)
    for name, n, ngt, p, fold in res["per_cohort"]:
        print("  %-12s n=%d  %d/%d  fold=%.1fx  P=%.4g" % (name, n, ngt, n, fold, p))
    cn, cgt, cp, cfold = res["combined"]
    print("  COMBINED     n=%d  %d/%d  fold=%.1fx  P=%.4g" % (cn, cgt, cn, cfold, cp))
    print("wrote %s.{png,pdf,svg} + %s_stats.csv" % (a.out, a.out))


def cmd_qc(a):
    rows = contamination.run(io.load_config(a.config), _kv(a.markers), _kv(a.target), a.out)
    for name, n, rho, p, ratio in rows:
        print("  %-12s n=%d  rho=%+.3f  P=%.3f  contam/tissue=%.3f" % (name, n, rho, p, ratio))
    print("wrote %s.{png,pdf,svg} + %s_scores.csv" % (a.out, a.out))


def cmd_selftest(a):
    from . import _selftest
    return _selftest.main()


def build_parser():
    p = argparse.ArgumentParser(prog="isoform-dominance",
                                description="Isoform-usage quantification and discrimination from bulk RNA-seq.")
    p.add_argument("--version", action="version", version="isoform-dominance %s" % __version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("annotate", help="gene symbol -> proposed isoform groups (Ensembl)")
    s.add_argument("--gene", required=True); s.add_argument("--species", default="homo_sapiens")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_annotate)

    s = sub.add_parser("identifiability", aliases=["identify"],
                       help="are the groups distinguishable by short reads?")
    s.add_argument("--config", required=True); s.add_argument("--k", type=int, default=31)
    s.add_argument("--sequences", help="optional JSON {transcript_id: cdna} (offline)")
    s.set_defaults(func=cmd_identifiability)

    s = sub.add_parser("extract", help="quant.sf -> per-donor isoform-group TPM")
    for x in ("config", "quantdir", "samplemap", "cohort", "out"):
        s.add_argument("--" + x, required=True)
    s.set_defaults(func=cmd_extract)

    s = sub.add_parser("stats", help="paired Wilcoxon + figure")
    s.add_argument("--config", required=True); s.add_argument("--condition", default="control")
    s.add_argument("--perdonor", action="append", required=True, help="NAME=perdonor.csv (repeatable)")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_stats)

    s = sub.add_parser("qc", help="contamination control")
    s.add_argument("--config", required=True)
    s.add_argument("--markers", action="append", required=True, help="NAME=marker_tpm.csv")
    s.add_argument("--target", action="append", required=True, help="NAME=perdonor.csv")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_qc)

    s = sub.add_parser("selftest", help="run the download-free reproducibility test")
    s.set_defaults(func=cmd_selftest)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
