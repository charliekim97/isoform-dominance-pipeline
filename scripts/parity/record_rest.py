"""Record what Ensembl REST answers for a list of genes at one release, for compare.py.

    python scripts/parity/record_rest.py --release 116 --genes genes.txt --out rest116.json.gz

``--genes`` names one symbol per line (``#`` starts a comment).  For each symbol the file
holds what ``annotate`` reads: ``lookup/symbol`` (which gene it gives, ``lookup_symbol``),
``xrefs/symbol`` (the gene ids it lists, ``xrefs``), and ``lookup/id?expand=1`` of every
one of those genes (``lookups``), trimmed to the fields the package reads; and the MD5 of
the cDNA of every transcript of every such gene (``cdna_md5``).  The format is that of
``tests/data/gencode_mini/rest116_mini.json.gz``.

This needs the network, and a release that REST or its archive still serves.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import sys
import urllib.parse

from isoform_dominance import ensembl

GENE_KEYS = ("id", "version", "display_name", "seq_region_name", "start", "end", "strand",
             "biotype", "canonical_transcript")
TRANSCRIPT_KEYS = ("id", "version", "biotype", "is_canonical")


def trim(g):
    """A ``lookup/id?expand=1`` gene as the package reads it."""
    out = {k: g.get(k) for k in GENE_KEYS}
    out["Transcript"] = [
        dict({k: t.get(k) for k in TRANSCRIPT_KEYS},
             Translation={"length": t["Translation"]["length"]} if t.get("Translation")
             else None,
             Exon=[{"start": e["start"], "end": e["end"]} for e in t.get("Exon", [])])
        for t in g.get("Transcript", [])]
    return out


def symbols(path):
    with open(path) as fh:
        return [w for w in (x.split("#")[0].strip() for x in fh) if w]


def record(genes, release, species="homo_sapiens", **net):
    server = ensembl.resolve_server(release, **net)
    net = dict(net, server=server)
    got = ensembl.release_number(ensembl.get_json("/info/data", **net))
    doc = {"release": got, "lookup_symbol": {}, "xrefs": {}, "lookups": {}, "cdna_md5": {},
           "source": "%s lookup/symbol|id?expand=1, xrefs/symbol and sequence/id?type=cdna, "
                     "fetched %s by scripts/parity/record_rest.py; trimmed to the fields the "
                     "package reads" % (server, datetime.date.today().isoformat())}
    q = urllib.parse.quote
    for sym in genes:
        g = ensembl.get_json("/lookup/symbol/%s/%s?expand=1" % (q(species), q(sym, safe="")),
                             **net)
        doc["lookup_symbol"][sym] = g["id"]
        doc["lookups"][g["id"]] = trim(g)
        xr = ensembl.get_json("/xrefs/symbol/%s/%s?object_type=gene"
                              % (q(species), q(sym, safe="")), **net)
        ids = sorted({x["id"] for x in xr if x.get("type") == "gene"
                      and str(x.get("id", "")).startswith("ENS")})
        doc["xrefs"][sym] = ids
        others = [i for i in ids if i not in doc["lookups"]]
        if others:
            more = ensembl.request_json("/lookup/id", {"ids": others, "expand": 1}, **net)
            doc["lookups"].update((i, trim(r)) for i, r in more.items() if r)
        print("%s: %s" % (sym, ", ".join(ids) or g["id"]), file=sys.stderr)
    tids = sorted({t["id"].split(".")[0] for g in doc["lookups"].values()
                   for t in g["Transcript"]})
    cdna = ensembl.fetch_cdna_batch(tids, **net)
    doc["cdna_md5"] = {t: hashlib.md5(s.encode()).hexdigest() for t, s in sorted(cdna.items())}
    return doc


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--release", type=int, required=True)
    p.add_argument("--genes", required=True, help="one symbol per line")
    p.add_argument("--out", required=True, help="JSON, gzipped when it ends in .gz")
    p.add_argument("--species", default="homo_sapiens")
    a = p.parse_args(argv)
    doc = record(symbols(a.genes), a.release, a.species)
    data = json.dumps(doc, sort_keys=True).encode()
    with open(a.out, "wb") as fh:
        fh.write(gzip.compress(data, mtime=0) if a.out.endswith(".gz") else data)
    print("release %s: %d genes, %d lookups, %d cDNA -> %s"
          % (doc["release"], len(doc["lookup_symbol"]), len(doc["lookups"]),
             len(doc["cdna_md5"]), a.out), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
