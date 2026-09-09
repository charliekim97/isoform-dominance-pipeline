"""Unified command-line interface: isoform-dominance <subcommand>.

Every subcommand accepts ``--json``, which writes the full result object to stdout
instead of the human summary, so the tool composes inside a workflow manager without
anyone having to parse its printed text.

``identifiability`` distinguishes three outcomes in its exit status, because the
middle one is the interesting case and a boolean cannot carry it:

===== ==============================================================
  0   identifiable
  3   weakly identifiable -- estimable, but ill-conditioned or
      starved of informative fragments at the stated design
  2   not identifiable -- no depth fixes this; regroup, or use
      long reads
===== ==============================================================
"""
import argparse
import json
import sys
from urllib.error import HTTPError, URLError

from . import (__version__, annotate, contamination, extract, identifiability, io,
               stats)

EXIT_OK = 0
EXIT_NOT_IDENTIFIABLE = 2
EXIT_WEAK = 3

_VERDICT_EXIT = {
    "identifiable": EXIT_OK,
    "weakly_identifiable": EXIT_WEAK,
    "not_identifiable": EXIT_NOT_IDENTIFIABLE,
}


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


def _emit(payload):
    json.dump(payload, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")


def _net_fail(e):
    print("Ensembl request failed (%s). Check your network connection, the gene symbol, "
          "and species, or supply sequences offline." % e, file=sys.stderr)
    return 1


def cmd_annotate(a):
    try:
        cfg = annotate.run(a.gene, a.out, species=a.species)
    except (URLError, HTTPError) as e:
        return _net_fail(e)
    if a.json:
        _emit(cfg)
        return EXIT_OK
    print("Proposed groups for %s -> %s" % (a.gene, a.out))
    for g, ids in cfg["groups"].items():
        print("  %s: %d transcripts" % (g, len(ids)))
    print("  primary_comparison:", cfg["primary_comparison"])
    print("  REVIEW _proposed/_clusters and rename groups before use.")
    return EXIT_OK


def cmd_identifiability(a):
    cfg = io.load_config(a.config)
    seqs = json.load(open(a.sequences)) if a.sequences else None
    background = json.load(open(a.background_sequences)) if a.background_sequences else None
    try:
        res = identifiability.analyze(
            cfg, k=a.k, sequences=seqs,
            canonical=not a.strand_aware,
            window=a.window,
            background_sequences=background,
            background_fasta=a.background_fasta,
            background_gene_transcripts=False if a.no_gene_background else "auto",
            read_length=a.read_length, frag_mean=a.frag_mean, frag_sd=a.frag_sd,
            paired=not a.single_end, depth=a.depth, tpm=a.tpm, n_donors=a.donors,
            conditioning_tau=a.tau, min_informative_reads=a.min_informative_reads)
    except (URLError, HTTPError) as e:
        return _net_fail(e)
    except ValueError as e:
        print("config error: %s" % e, file=sys.stderr)
        return 1

    if a.json:
        _emit(res)
        return _VERDICT_EXIT[res["verdict"]]

    bg = res["background"]
    print("Identifiability (window=%d, k=%d, %s k-mers)"
          % (res["window"], res["k"], "canonical" if res["canonical"] else "strand-aware"))
    print("  background: %d same-gene transcript(s)%s"
          % (bg["n_background_transcripts"],
             ", FASTA %s" % bg["fasta"] if bg["fasta"] else ""))
    if not bg["fasta"]:
        scope = ("this gene's other transcripts"
                 if bg["n_background_transcripts"] else "the configured groups only")
        print("  NOTE: uniqueness judged against %s. A quantifier resolves fragments "
              "against the whole index, so pseudogenes, paralogues and homologous "
              "loci outside this gene are not accounted for here. Pass "
              "--background-fasta <the FASTA the Salmon index was built from> for the "
              "answer that matches what the quantifier actually sees; that is the "
              "recommended way to run this command." % scope, file=sys.stderr)
    d = res["design"]
    print("  design: %s %dbp reads, fragments %.0f+-%.0f, depth %.0fM, TPM %.3g, n=%d"
          % ("paired" if d["paired"] else "single", d["read_length"],
             d["frag_mean"], d["frag_sd"], d["depth"] / 1e6, d["tpm"], d["n_donors"]))
    for g, r in res["groups"].items():
        print("  [%s] %s: %d unique k-mers, %d bp in %d block(s), "
              "~%.0f informative reads, conditioning %.2f"
              % (r["verdict"], g, r["n_unique_kmers"], r["unique_length"],
                 r["n_blocks"], r["expected_informative_reads"], r["conditioning_factor"]))
    c = res["contrast"]
    print("  contrast %s vs %s: %s (conditioning %.2f)"
          % (res["primary_comparison"][0], res["primary_comparison"][1],
             c["verdict"], c["conditioning_factor"]))
    noise = res["counting_noise"]
    print("  counting-noise floor on log2 ratio: SE %.3f per donor, "
          "min resolvable |log2FC| %.3f at n=%d"
          % (noise["log2_ratio_se"], noise["min_resolvable_log2fc"], d["n_donors"]))
    print("  VERDICT:", res["verdict"])
    for reason in res["reasons"]:
        print("    - %s" % reason, file=sys.stderr)
    if res["verdict"] == "not_identifiable":
        print("  No sequencing depth fixes this: the class contrast is outside the "
              "row space of the compatibility system. Regroup, or use long reads.",
              file=sys.stderr)
    return _VERDICT_EXIT[res["verdict"]]


def cmd_extract(a):
    n = extract.run(io.load_config(a.config), a.quantdir, a.samplemap, a.cohort, a.out)
    if a.json:
        _emit({"out": a.out, "n_donors": n, "cohort": a.cohort})
    else:
        print("wrote %s (n=%d donors)" % (a.out, n))
    return EXIT_OK


def cmd_stats(a):
    res = stats.run(io.load_config(a.config), a.condition, _kv(a.perdonor), a.out,
                    n_boot=a.n_boot, seed=a.seed)
    if a.json:
        _emit(res)
        return EXIT_OK
    for det in res["detail"]:
        line = ("  %-12s n=%d  %d/%d  fold=%.1fx [%.1f-%.1f]  P=%.4g"
                % (det["cohort"], det["n"], det["n_greater"], det["n"],
                   det["median_fold"], det["fold_ci"][0], det["fold_ci"][1], det["p"]))
        if det["underpowered"]:
            line += "  (floor %.4g: cannot reach 0.05)" % det["p_floor"]
        print(line)
    cn, cgt, cp, cfold = res["combined"]
    print("  %-12s n=%d  %d/%d  fold=%.1fx  P=%.4g" % ("POOLED", cn, cgt, cn, cfold, cp))
    for key in ("stouffer", "stratified_signed_rank"):
        c = res["combination"][key]
        print("  %-12s k=%d cohorts  P=%.4g" % (key.upper(), c["k"], c["p"]))
    print("  headline combination: %s (pooling donors across independent cohorts "
          "ignores the cohort factor)" % res["headline_combination"])
    print("wrote %s.{png,pdf,svg} + %s_stats.csv" % (a.out, a.out))
    return EXIT_OK


def cmd_qc(a):
    rows = contamination.run(io.load_config(a.config), _kv(a.markers), _kv(a.target), a.out)
    if a.json:
        _emit([{"cohort": r[0], "n": r[1], "rho": r[2], "p": r[3], "ratio": r[4]}
               for r in rows])
        return EXIT_OK
    for name, n, rho, p, ratio in rows:
        print("  %-12s n=%d  rho=%+.3f  P=%.3f  contam/tissue=%.3f" % (name, n, rho, p, ratio))
    print("wrote %s.{png,pdf,svg} + %s_scores.csv" % (a.out, a.out))
    return EXIT_OK


def cmd_selftest(a):
    from . import _selftest
    return _selftest.main()


def build_parser():
    p = argparse.ArgumentParser(prog="isoform-dominance",
                                description="Isoform-usage quantification and discrimination from bulk RNA-seq.")
    p.add_argument("--version", action="version", version="isoform-dominance %s" % __version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def _json(sp):
        sp.add_argument("--json", action="store_true",
                        help="emit the full result as JSON on stdout")
        return sp

    s = _json(sub.add_parser("annotate", help="gene symbol -> proposed isoform groups (Ensembl)"))
    s.add_argument("--gene", required=True); s.add_argument("--species", default="homo_sapiens")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_annotate)

    s = _json(sub.add_parser(
        "identifiability", aliases=["identify"],
        help="are the classes measurable by short reads, and how precisely?"))
    s.add_argument("--config", required=True)
    s.add_argument("--k", type=int, default=identifiability.DEFAULT_K)
    s.add_argument("--window", type=int, default=None,
                   help="window length for the compatibility system (default: k; "
                        "set to the read length for a sharper, still conservative, system)")
    s.add_argument("--sequences", help="optional JSON {transcript_id: cdna} (offline)")
    s.add_argument("--background-sequences", help="optional JSON {transcript_id: cdna} of background transcripts")
    s.add_argument("--background-fasta",
                   help="FASTA (optionally gzipped) to judge uniqueness against -- "
                        "ideally the one the Salmon index was built from")
    s.add_argument("--no-gene-background", action="store_true",
                   help="do not fetch the gene's other transcripts as background (v2.1 behaviour)")
    s.add_argument("--strand-aware", action="store_true",
                   help="do not fold k-mers to their canonical form (v2.1 behaviour)")
    s.add_argument("--read-length", type=int, default=identifiability.DEFAULT_READ_LENGTH)
    s.add_argument("--frag-mean", type=float, default=identifiability.DEFAULT_FRAG_MEAN)
    s.add_argument("--frag-sd", type=float, default=identifiability.DEFAULT_FRAG_SD)
    s.add_argument("--single-end", action="store_true")
    s.add_argument("--depth", type=float, default=identifiability.DEFAULT_DEPTH,
                   help="mapped fragments per library")
    s.add_argument("--tpm", type=float, default=identifiability.DEFAULT_TPM,
                   help="class abundance to condition the read model on")
    s.add_argument("--donors", type=int, default=1)
    s.add_argument("--tau", type=float, default=identifiability.DEFAULT_CONDITIONING_TAU,
                   help="structural conditioning factor above which a class is only weakly identifiable")
    s.add_argument("--min-informative-reads", type=float,
                   default=identifiability.DEFAULT_MIN_INFORMATIVE_READS)
    s.set_defaults(func=cmd_identifiability)

    s = _json(sub.add_parser("extract", help="quant.sf -> per-donor isoform-group TPM"))
    for x in ("config", "quantdir", "samplemap", "cohort", "out"):
        s.add_argument("--" + x, required=True)
    s.set_defaults(func=cmd_extract)

    s = _json(sub.add_parser("stats", help="paired Wilcoxon, cohort combination + figure"))
    s.add_argument("--config", required=True); s.add_argument("--condition", default="control")
    s.add_argument("--perdonor", action="append", required=True, help="NAME=perdonor.csv (repeatable)")
    s.add_argument("--n-boot", type=int, default=stats.DEFAULT_N_BOOT,
                   help="bootstrap replicates for the fold-change interval (0 disables)")
    s.add_argument("--seed", type=int, default=stats.DEFAULT_SEED)
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_stats)

    s = _json(sub.add_parser("qc", help="contamination control"))
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
