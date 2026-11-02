"""Per-transcript detail of the 3.0 gate's run, for the genes its write-up discusses.

    python validation/tpm_run_detail.py RUN_DIR INDEX_DIR CDNA_FASTA \
        --out validation/data/sim/c1_d1_tpm_detail.json

``RUN_DIR`` is where the run was made (``run_simulation.py run --work``): it holds the
simulator's truth of each replicate, ``rep{r}.truth.json``, and Salmon's ``quant.sf`` of
each, kept beside the run as ``quant_kept/quant{r}.quant.sf.gz`` (the harness itself
deletes them); its ``c1_d1_tpm.json`` must be the stored run, by SHA-256.  ``INDEX_DIR`` is
the stored Salmon index, whose ``duplicate_clusters.tsv`` names the records it dropped for
an identical copy, and ``CDNA_FASTA`` the FASTA it was built from, whose headers give the
gene symbol of a record.

For each gene of ``--genes`` (default LDHA, NDUFS1 and SMN2) and each replicate: the
fragments the simulator wrote for each transcript of the two compared classes, and
``quant.sf``'s ``EffectiveLength``, ``TPM`` and ``NumReads`` of each of them, and of each
copy the index kept in place of a dropped one, with the fragments written for the whole
gene.  For every replicate, the fragments written for the library, and the ``NumReads``
``quant.sf`` gives the records outside the 49 genes, summed by the gene symbol of their
FASTA header.  Some of those records are not another gene's at all: the index kept them in
place of a transcript of the 49 genes with the same sequence, which it dropped, so their
reads are that transcript's under another name and another symbol.  For each replicate,
``outside_stand_ins`` gives each such record's ``NumReads`` and symbol, and the transcripts
it stands for, with each one's gene and the fragments written for it.  Nothing here grades
or summarises; ``tests/test_gate_run.py`` recomputes these genes' rows from this file with
the harness's own truth and holds them to the run's.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RUN = "sim/c1_d1_tpm.json"
OUT = os.path.join(DATA, "sim", "c1_d1_tpm_detail.json")
GENES = ("LDHA", "NDUFS1", "SMN2")
QUANT_HEADER = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]


def _strip(t):
    return t.split(".")[0]


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_quant(path):
    """``{id without version: (EffectiveLength, TPM, NumReads)}`` of a gzipped quant.sf."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        header = next(f).rstrip("\n").split("\t")
        if header[:5] != QUANT_HEADER:
            raise SystemExit("%s: header %r is not Salmon's" % (path, header))
        for line in f:
            p = line.rstrip("\n").split("\t")
            k = _strip(p[0])
            if k in out:
                raise SystemExit("%s: two records are %s without their version" % (path, k))
            out[k] = (float(p[2]), float(p[3]), float(p[4]))
    return out


def symbols(fasta):
    """``{record id without version: gene symbol}`` from the FASTA's headers; a record
    without ``gene_symbol:`` goes by its ``gene:`` id."""
    out = {}
    with gzip.open(fasta, "rt", encoding="utf-8") as f:
        for line in f:
            if line[0] == ">":
                m = re.search(r" gene_symbol:(\S+)", line) or re.search(r" gene:(\S+)", line)
                out[_strip(line[1:].split()[0])] = m.group(1) if m else None
    return out


def kept_copies(index):
    """``{dropped: kept}`` from ``duplicate_clusters.tsv``, as ``run_simulation`` reads it."""
    out = {}
    with open(os.path.join(index, "duplicate_clusters.tsv"), encoding="utf-8") as f:
        next(f)
        for line in f:
            kept, dropped = line.rstrip("\n").split("\t")[:2]
            out[_strip(dropped)] = _strip(kept)
    return out


def stand_ins(kept, panel):
    """``{record: [transcripts]}``: each record outside ``panel`` that the index kept in
    place of transcripts of ``panel`` it dropped for the same sequence, with those
    transcripts.  ``kept`` is :func:`kept_copies`'s ``{dropped: kept}``."""
    out = {}
    for dropped, k in kept.items():
        if dropped in panel and k not in panel:
            out.setdefault(k, []).append(dropped)
    return {k: sorted(v) for k, v in sorted(out.items())}


def detail(run_dir, index, fasta, genes=GENES, data=DATA):
    with open(os.path.join(data, "README.md"), encoding="utf-8") as f:
        want = re.search(r"^\| `%s` \|.*\| `([0-9a-f]{64})` \|$" % re.escape(RUN), f.read(),
                         re.M)
    got = file_sha256(os.path.join(run_dir, os.path.basename(RUN)))
    if not want or got != want.group(1):
        raise SystemExit("%s/%s is not the run data/README.md names"
                         % (run_dir, os.path.basename(RUN)))
    with open(os.path.join(run_dir, os.path.basename(RUN)), encoding="utf-8") as f:
        run = json.load(f)
    with gzip.open(os.path.join(data, "sweep116.json.gz"), "rt", encoding="utf-8") as f:
        source = json.load(f)["genes"]
    gene_of = {_strip(t): g for g, e in source.items() for t in e["sequences"]}
    panel = set(gene_of)
    kept = kept_copies(index)
    carriers = stand_ins(kept, panel)
    names = symbols(fasta)
    out = {"run": {"file": RUN, "sha256": got}, "genes": {}, "fragments_written": [],
           "outside_num_reads": [], "outside_stand_ins": []}
    for gene in genes:
        cfg = source[gene]["config"]
        pc = list(cfg["primary_comparison"])
        groups = {c: [_strip(t) for t in cfg["groups"][c]] for c in pc}
        out["genes"][gene] = {"primary_comparison": pc, "classes": groups,
                              "kept_copies": {}, "replicates": []}
    for rep in range(run["design"]["replicates"]):
        with open(os.path.join(run_dir, "rep%d.truth.json" % rep), encoding="utf-8") as f:
            truth = json.load(f)
        quant = read_quant(os.path.join(run_dir, "quant_kept", "quant%d.quant.sf.gz" % rep))
        for gene, g in out["genes"].items():
            members = [t for c in g["primary_comparison"] for t in g["classes"][c]]
            for t in members:
                if t not in quant and t in kept:
                    g["kept_copies"][t] = kept[t]
            listed = sorted({t for t in members if t in quant} | set(g["kept_copies"].values()))
            written = truth["genes"][gene]["fragments_per_transcript"]
            g["replicates"].append({
                "fragments": {t: written.get(t, 0) for t in members},
                "gene_fragments_written": truth["genes"][gene]["fragments_written"],
                "quant": {t: list(quant[t]) for t in listed}})
        out["fragments_written"].append(truth["fragments_written"])
        outside = {}
        for k, (_e, _t, n) in quant.items():
            if n > 0 and k not in panel:
                s = names.get(k) or k
                outside[s] = outside.get(s, 0.0) + n
        out["outside_num_reads"].append({s: round(n, 3) for s, n in sorted(outside.items())})
        out["outside_stand_ins"].append({
            k: {"num_reads": round(quant[k][2], 3), "symbol": names.get(k) or k,
                "stands_for": {t: {"gene": gene_of[t], "fragments_written":
                                   truth["genes"][gene_of[t]]["fragments_per_transcript"]
                                   .get(t, 0)}
                               for t in stands}}
            for k, stands in carriers.items()})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("index")
    ap.add_argument("fasta")
    ap.add_argument("--genes", default=",".join(GENES))
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    doc = detail(a.run_dir, a.index, a.fasta, a.genes.split(","))
    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n")
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
