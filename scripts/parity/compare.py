"""Hold the GTF path to REST: every config, and every transcript's cDNA.

    python scripts/parity/compare.py --gtf gencode.v50.annotation.gtf.gz \\
        --fasta gencode.v50.transcripts.fa.gz --rest rest116.json.gz [--genes genes.txt] \\
        [--json result.json]

``--rest`` is what record_rest.py records (or ``tests/data/gencode_mini/rest116_mini.json.gz``):
REST's answers at one release.  For each gene, ``annotate.build_config`` answered from it
is compared with ``annotate.build_config_from_gtf`` on ``--gtf``, ``annotation_source``
aside; the transcripts of every biotype the two give the gene the GTF chose; and the cDNA
of every one of them, from ``--fasta``, with REST's by MD5.  The GTF and the FASTA must be the release ``--rest`` was recorded at
(GENCODE 44, 48 and 50 are Ensembl 110, 114 and 116).

The lines of the genes asked about, and their FASTA records, are first cut from the files
in one pass each, as ``tests/data/gencode_mini`` was; ``--whole`` reads the files
themselves for every gene instead, three passes a gene.  Without ``--genes`` the genes are
the symbols ``--rest`` was recorded for.  Exit status 0 when every config and every
sequence is equal.
"""
import argparse
import hashlib
import json
import os
import sys
import tempfile
import time

from isoform_dominance import annotate
from isoform_dominance import annotation_files as af
from isoform_dominance import io


def load(path):
    with io.open_text(path) as fh:
        return json.load(fh)


def answer_from(rest):
    """``annotate``'s ``_get`` and ``_post``, answered from a recorded release."""
    names = {}
    for g in rest["lookups"].values():
        names.setdefault(g["display_name"], []).append(g["id"])

    def get(path, **kw):
        part = path.split("?")[0].split("/")
        if path.startswith("/info/data"):
            return {"releases": [rest["release"]]}
        if path.startswith("/lookup/symbol/"):
            return rest["lookups"][rest["lookup_symbol"].get(part[4]) or names[part[4]][0]]
        if path.startswith("/xrefs/symbol/"):
            ids = (rest["xrefs"].get(part[4], []) if "xrefs" in rest
                   else names.get(part[4], []))
            return [{"type": "gene", "id": i} for i in ids]
        if path.startswith("/lookup/id/"):
            return rest["lookups"][part[3]]
        raise KeyError("no recorded answer for %s" % path)

    def post(path, body, **kw):
        return {i: rest["lookups"][i] for i in body["ids"] if i in rest["lookups"]}
    return get, post


def cut(gtf, fasta, symbols, gene_ids, where):
    """The header and every line of ``gtf`` that names one of ``symbols`` or ``gene_ids``,
    and the records of ``fasta`` for those lines' transcripts: ``(gtf, fasta)`` paths."""
    needles = [('gene_name "%s";' % s).encode() for s in symbols]
    needles += [('gene_id "%s.' % g).encode() for g in gene_ids]
    needles += [('gene_id "%s"' % g).encode() for g in gene_ids]

    def found(text):
        return any(n in text for n in needles)
    out_gtf = os.path.join(where, "cut.gtf")
    tids = set()
    with open(out_gtf, "w") as out, io.open_text(gtf) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            out.write(line)
        for line in af._matching_lines(gtf, found, af.BLOCK, {"lines": 0, "transcripts": 0,
                                                                "basic": 0}):
            out.write(line.rstrip("\r\n") + "\n")
            if 'transcript_id "' in line:
                tids.add(line.split('transcript_id "', 1)[1].split('"', 1)[0].split(".")[0])
    out_fa = os.path.join(where, "cut.fa")
    keep = False
    with open(out_fa, "w") as out, io.open_text(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                first = af.fasta_id(line[1:])
                keep = bool(first) and first.split(".")[0] in tids
            if keep:
                out.write(line)
    return out_gtf, out_fa


def compare(gene, gtf, fasta, rest):
    """One gene: whether the configs are equal, whether the GTF and REST give the gene the
    same transcripts -- of every biotype, the gene background -- and the cDNA compared and
    differing."""
    row = {"gene": gene}
    try:
        want = annotate.build_config(gene)
    except Exception as e:                       # recorded, not raised
        want = "%s: %s" % (type(e).__name__, e)
    t0 = time.perf_counter()
    try:
        got = annotate.build_config_from_gtf(gene, gtf)
        got.pop("annotation_source")
    except Exception as e:
        got = "%s: %s" % (type(e).__name__, e)
    row["seconds"] = round(time.perf_counter() - t0, 2)
    row["config_equal"] = got == want
    if not row["config_equal"]:
        row["differ"] = (sorted(k for k in set(got) | set(want) if got.get(k) != want.get(k))
                         if isinstance(got, dict) and isinstance(want, dict)
                         else {"rest": want if isinstance(want, str) else "config",
                               "gtf": got if isinstance(got, str) else "config"})
    row["cdna_compared"], row["cdna_differ"] = 0, []
    if isinstance(got, dict):
        rec = next(g for g in af.scan(gtf, gene_id=got["gene_id"]) if not g["_par_y"])
        mine = {t["id"] for t in rec["Transcript"]}
        theirs = {t["id"].split(".")[0]
                  for t in rest["lookups"].get(got["gene_id"], {}).get("Transcript", [])}
        if mine != theirs:
            row["transcripts_differ"] = {"rest_only": sorted(theirs - mine),
                                         "gtf_only": sorted(mine - theirs)}
        try:
            seqs = af.sequences_for(rec, fasta)
        except af.AnnotationFileError as e:
            row["cdna_error"] = str(e)
            seqs = {}
        for t, s in sorted(seqs.items()):
            if t in rest["cdna_md5"]:
                row["cdna_compared"] += 1
                if hashlib.md5(s.encode()).hexdigest() != rest["cdna_md5"][t]:
                    row["cdna_differ"].append(t)
    return row


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--gtf", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--rest", required=True, help="REST answers recorded by record_rest.py")
    p.add_argument("--genes", help="one symbol per line (default: those --rest holds)")
    p.add_argument("--whole", action="store_true",
                   help="read the whole files for every gene rather than cutting them first")
    p.add_argument("--json", help="write every gene's row here")
    a = p.parse_args(argv)
    rest = load(a.rest)
    if a.genes:
        with open(a.genes) as fh:
            genes = [w for w in (x.split("#")[0].strip() for x in fh) if w]
    else:
        genes = sorted(rest.get("xrefs") or {g["display_name"]
                                              for g in rest["lookups"].values()})
    get, post = answer_from(rest)
    annotate._get, annotate._post = get, post
    annotate.ensembl.resolve_server = lambda release=None, **kw: None
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory() as where:
        gtf, fasta = (a.gtf, a.fasta) if a.whole else cut(
            a.gtf, a.fasta, genes, sorted(rest["lookups"]), where)
        cut_s = time.perf_counter() - t0
        rows = [compare(g, gtf, fasta, rest) for g in genes]
    for r in rows:
        moved = r.get("transcripts_differ")
        print("%-10s config %s  %scDNA %d compared%s%s"
              % (r["gene"], "equal" if r["config_equal"] else "DIFFERS %s" % r["differ"],
                 "transcripts DIFFER: %s  " % "; ".join(
                     "%s only %s" % (side, ", ".join(moved[key]))
                     for side, key in (("REST", "rest_only"), ("GTF", "gtf_only"))
                     if moved[key]) if moved else "",
                 r["cdna_compared"],
                 ", %d DIFFER: %s" % (len(r["cdna_differ"]), ", ".join(r["cdna_differ"]))
                 if r["cdna_differ"] else "", "; %s" % r["cdna_error"]
                 if r.get("cdna_error") else ""))
    n_eq = sum(r["config_equal"] for r in rows)
    n_tx = sum("transcripts_differ" not in r for r in rows)
    n_seq = sum(r["cdna_compared"] for r in rows)
    n_bad = sum(len(r["cdna_differ"]) for r in rows)
    print("release %s: configs equal %d/%d (annotation_source aside); transcript sets equal "
          "%d/%d; cDNA byte-identical %d/%d; %s %.1f s"
          % (rest["release"], n_eq, len(rows), n_tx, len(rows), n_seq - n_bad, n_seq,
             "reading the whole files" if a.whole else "cutting the files", cut_s))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({"release": rest["release"], "rows": rows}, fh, indent=1)
    return 0 if n_eq == n_tx == len(rows) and not n_bad and not any(
        r.get("cdna_error") for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
