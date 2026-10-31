"""The 49 genes' contrasts against the whole index the simulation was quantified with.

    python validation/whole_index.py Homo_sapiens.GRCh38.cdna.all.fa.gz \
        --out validation/data/whole_index.json

For each gene of ``data/sweep116.json.gz``, three backgrounds:

- ``gene_background``: the call the harness made before any read existed
  (``run_simulation.predicted``: the configured transcripts, the gene's other transcripts as
  the background, the design the simulation drew with);
- ``whole_index``: the same call with ``background_fasta`` set to the transcript FASTA the
  stored Salmon index was built from, which ``identifiability.analyze`` calls the honest
  whole-index answer -- Ensembl's ``cdna.all``, its records on scaffolds (alternate loci,
  patches and a few others) included;
- ``reference_chromosomes``: the same with only that FASTA's records on chromosomes 1-22, X, Y
  and MT, the records the README's "Index scope" keeps to build an index.

Both estimands, as ``analyze`` reports them, and for the two FASTA backgrounds the records
that became columns, each with its gene symbol and location from its FASTA header.  The
FASTA must be the one ``data/README.md`` names for the index, by SHA-256: against another
file the comparison describes another index.

A gene takes about 50 s, each FASTA being read twice for it.  ``--genes`` runs some genes
only, and ``--from`` takes the genes an earlier run of this script wrote with the same
FASTA, inputs, design and package version, so a run can be made in parts: the output is
the same whichever way it was made.

CI does not run this on the index's FASTA, which it does not have; ``tests/test_gate_run.py``
holds the stored output's gene-background half to what the package gives now, and runs the
script on a FASTA of three records.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
import tempfile

from isoform_dominance import __version__
from isoform_dominance import identifiability as I
from isoform_dominance import io

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(DATA, "whole_index.json")
INPUTS = "sweep116.json.gz"
#: the design the simulation drew with (data/README.md), as run_simulation.predicted passes it
DESIGN = {"tpm": 10.0, "depth": 30e6, "mean_efflen": 1500.0}
ESTIMANDS = ("fragment", "abundance")
#: the regions an index of reference-chromosome transcripts keeps: the README's "Index scope"
#: filter of Ensembl's cdna.all, on the third field of the header's region
REFERENCE = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])


def _summarize():
    spec = importlib.util.spec_from_file_location("_summarize",
                                                  os.path.join(HERE, "summarize.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def index_fasta(data=DATA):
    """``(name, bytes, SHA-256)`` of the FASTA ``data/README.md`` says the index was built
    from."""
    with open(os.path.join(data, "README.md"), encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"of Ensembl 116's `([^`]+)` \(([\d,]+)\s+bytes, SHA-256 `([0-9a-f]{64})`\)",
                  text)
    if not m:
        raise SystemExit("%s/README.md names no FASTA for the index" % data)
    return m.group(1), int(m.group(2).replace(",", "")), m.group(3)


def reference_fasta(fasta, out):
    """Write to ``out``, gzipped, the records of ``fasta`` on a reference chromosome: those
    whose header region (``chromosome:GRCh38:<name>:...``) names 1-22, X, Y or MT.  Returns
    ``{"records", "sha256"}``, the SHA-256 of the records' text, as written before
    compression."""
    h, n, keep = hashlib.sha256(), 0, False
    with gzip.open(fasta, "rt", encoding="utf-8") as f, \
            gzip.open(out, "wt", encoding="utf-8", compresslevel=1) as g:
        for line in f:
            if line[0] == ">":
                fields = line[1:].split()
                region = fields[2].split(":") if len(fields) > 2 else []
                keep = len(region) > 2 and region[2] in REFERENCE
                n += keep
            if keep:
                g.write(line)
                h.update(line.encode("utf-8"))
    return {"records": n, "sha256": h.hexdigest()}


def header(fasta, data=DATA):
    """What a run is made with: the FASTA (refused unless it is the index's), the inputs,
    the design and the package version.  Runs whose headers agree give the same genes."""
    s = _summarize()
    name, size, sha = index_fasta(data)
    got = s.file_sha256(fasta)
    if got != sha:
        raise SystemExit("%s has SHA-256 %s, not the %s of %s, the FASTA data/README.md says "
                         "the index was built from" % (fasta, got, sha, name))
    inputs = s.file_sha256(os.path.join(data, INPUTS))
    if inputs != s.sha256_table(os.path.join(data, "README.md")).get(INPUTS):
        raise SystemExit("%s is not the file data/README.md names" % INPUTS)
    return {"fasta": {"file": name, "bytes": size, "sha256": sha},
            "inputs": {"file": INPUTS, "sha256": inputs},
            "design": dict(DESIGN), "isoform_dominance": __version__}


def figures(res):
    """Per estimand, what ``analyze`` says of the contrast, and the predicted SE the harness
    records beside it (``min_resolvable_log2fc`` / 1.96)."""
    out = {}
    for name in ESTIMANDS:
        e = res["estimands"][name]
        c = e["contrast"]
        fc = c["min_resolvable_log2fc"]
        out[name] = {"verdict": e["verdict"], "estimable": c["estimable"],
                     "conditioning_factor": c["conditioning_factor"],
                     "min_resolvable_log2fc": fc,
                     "predicted_se": fc / 1.96 if math.isfinite(fc) else math.inf,
                     "beyond_linear": c["beyond_linear"]}
    return out


def _call(entry):
    """``(config, keyword arguments)`` of ``run_simulation.predicted``'s call for a gene."""
    cfg = entry["config"]
    seqs = {t.split(".")[0]: s for t, s in entry["sequences"].items()}
    need = {t.split(".")[0] for v in cfg["groups"].values() for t in v}
    return cfg, {"sequences": {t: seqs[t] for t in need},
                 "background_sequences": {t: s for t, s in seqs.items() if t not in need},
                 "background_gene_transcripts": False, "tpm": DESIGN["tpm"],
                 "depth": DESIGN["depth"], "mean_efflen": DESIGN["mean_efflen"]}


def gene_background(entry):
    """The gene's figures with its other transcripts as the background."""
    cfg, kw = _call(entry)
    return figures(I.analyze(cfg, **kw))


def against(entry, fasta):
    """The gene's figures with ``fasta`` as the background, and the records that became
    columns, ``{id: windows shared with the configured transcripts}``."""
    cfg, kw = _call(entry)
    res = I.analyze(cfg, background_fasta=fasta, **kw)
    b = res["background"]
    out = figures(res)
    out.update(fasta_competitors=b["fasta_competitors"], fasta_left_out=b["fasta_left_out"],
               windows_dropped={k: b["windows_dropped"][k]
                                for k in ("total", "configured", "configured_of")})
    return out


def describe(fasta, ids):
    """``{id: {"gene_symbol", "location"}}`` of FASTA records, from their headers: the
    symbol after ``gene_symbol:`` and the first three fields of the location
    (``chromosome:GRCh38:5``, ``scaffold:GRCh38:HSCHR5_1_CTG1_1``)."""
    want, out = set(ids), {}
    if not want:
        return out
    with gzip.open(fasta, "rt", encoding="utf-8") as f:
        for line in f:
            if line[0] != ">":
                continue
            fields = line[1:].split()
            rid = fields[0].split(".")[0]
            if rid in want:
                sym = re.search(r" gene_symbol:(\S+)", line)
                out[rid] = {"gene_symbol": sym.group(1) if sym else None,
                            "location": ":".join(fields[2].split(":")[:3])
                            if len(fields) > 2 else None}
    return out


def run(fasta, genes=None, parts=(), data=DATA):
    """The output document: :func:`header`, the reference-chromosome records, and per gene
    its three backgrounds' figures."""
    head = header(fasta, data)
    with gzip.open(os.path.join(data, INPUTS), "rt", encoding="utf-8") as f:
        source = json.load(f)["genes"]
    with tempfile.TemporaryDirectory() as tmp:
        primary = os.path.join(tmp, "reference_chromosomes.fa.gz")
        head["reference_chromosomes"] = reference_fasta(fasta, primary)
        done = {}
        for path in parts:
            with open(path, encoding="utf-8") as f:
                doc = json.load(f)
            if {k: v for k, v in doc.items() if k != "genes"} != head:
                raise SystemExit("%s was made with another FASTA, inputs, design or package "
                                 "version" % path)
            done.update(doc["genes"])
        todo = [g for g in (genes or sorted(source)) if g not in done]
        unknown = [g for g in todo if g not in source]
        if unknown:
            raise SystemExit("not in %s: %s" % (INPUTS, ", ".join(unknown)))
        new = {}
        for gene in todo:
            new[gene] = {"gene_background": gene_background(source[gene]),
                         "whole_index": against(source[gene], fasta),
                         "reference_chromosomes": against(source[gene], primary)}
            print("%s: abundance %s -> %s / %s, fragment %s -> %s / %s" % (
                gene, *(new[gene][b][e]["estimable"] for e in ("abundance", "fragment")
                        for b in ("gene_background", "whole_index", "reference_chromosomes"))),
                file=sys.stderr, flush=True)
    names = describe(fasta, {r for g in new.values()
                             for b in ("whole_index", "reference_chromosomes")
                             for r in g[b]["fasta_competitors"]})
    for g in new.values():
        for b in ("whole_index", "reference_chromosomes"):
            w = g[b]
            w["fasta_competitors"] = {r: dict(windows=n, **names.get(r, {}))
                                      for r, n in sorted(w["fasta_competitors"].items())}
    done.update(new)
    return dict(head, genes=done)


def as_json(doc):
    """The output's text: keys sorted, a number that is not finite as null."""
    return json.dumps(io.finite_json(doc), indent=1, sort_keys=True, allow_nan=False,
                      ensure_ascii=False) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="the transcript FASTA the stored index was built from")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--genes", help="comma-separated genes to run (default: all 49)")
    ap.add_argument("--from", dest="parts", action="append", default=[], metavar="JSON",
                    help="an earlier output of this script whose genes are taken as they are")
    a = ap.parse_args(argv)
    doc = run(a.fasta, a.genes.split(",") if a.genes else None, a.parts)
    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        f.write(as_json(doc))
    print("wrote %s: %d genes" % (a.out, len(doc["genes"])), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
