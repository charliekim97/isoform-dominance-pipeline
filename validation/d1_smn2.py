"""D1: did the index's records beyond SMN2's own transcripts produce SMN2's failure of the
3.0 gate?  ``validation/README.md``, "What produced SMN2's failure (D1)", fixes the design
and the reading; this file fixes the same in code, before the run.

    python validation/d1_smn2.py fasta CDNA_FASTA {ref,gene} OUT.fa.gz
    python validation/d1_smn2.py index {ref,gene} INDEX_DIR
    python validation/d1_smn2.py quant {ref,gene} QUANT_SF[.gz]
    python validation/d1_smn2.py grade RUN_REF RUN_GENE

``fasta`` writes Ensembl 116's ``cdna.all``, the FASTA the stored index was built from,
without the records an arm leaves out (:data:`LEFT_OUT`):

- ``ref``: the SMN1 and SMN2 records off chromosome 5, their copies on the alternate loci
  HSCHR5_1_CTG1_1 and HSCHR5_2_CTG1_1, which the README's "Index scope" leaves out of an
  index;
- ``gene``: every SMN1 record, and the SMN2 records off chromosome 5, so that of the two genes
  the index holds SMN2's 23 transcripts on chromosome 5 and nothing else, as the prediction
  did.

It refuses another input FASTA, and a run that would leave out other records or write other
text, by SHA-256.  ``index`` and ``quant`` check that an index built from that FASTA, and a
``quant.sf`` made against it, hold what they should: Salmon 2.7.0, k 31, no decoys,
duplicates not kept, :data:`NUM_REFS` records, none of them left out.

``grade`` reads the two runs, ``run_simulation.py run`` made as the gate's run was against
each arm's index, and grades SMN2's summed-TPM row in each as step 3 of the gate grades it:
an arm meets the predicted SE the gate fixed when its |bias| and its replicate SD are both
below it.  The exit status is 0 when ``gene`` meets it (the failure was the index's), 1 when
it does not (D1 does not attribute the failure to the index), and 2 when the runs are not
the ones the design names; then nothing is graded.  The ``NumReads`` rows, ``ref``'s grade
beside ``gene``'s, and the other 48 genes' rows against the stored run are reported and
decide nothing.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STORED_TPM_RUN = os.path.join(HERE, "data", "sim", "c1_d1_tpm.json")
GENE = "SMN2"
ARMS = ("ref", "gene")

#: Ensembl 116's Homo_sapiens.GRCh38.cdna.all.fa.gz, as validation/data/README.md names it
CDNA_SHA256 = "683eb19310c40bf1396e4718f45afa2ce86755717c0990f47a171f535d248ea1"
#: the records each arm leaves out, without their versions
LEFT_OUT = {
    "ref": (
        "ENST00000610814", "ENST00000611442", "ENST00000611906", "ENST00000612446",
        "ENST00000612974", "ENST00000614610", "ENST00000614773", "ENST00000617648",
        "ENST00000618251", "ENST00000618661", "ENST00000619658", "ENST00000620650",
        "ENST00000622158", "ENST00000622739", "ENST00000624213", "ENST00000624254",
        "ENST00000624634", "ENST00000625462", "ENST00000625687", "ENST00000625953",
        "ENST00000626580", "ENST00000626988", "ENST00000627306", "ENST00000627341",
        "ENST00000628000", "ENST00000628043", "ENST00000628353", "ENST00000628642",
        "ENST00000629122", "ENST00000629741", "ENST00000630432", "ENST00000631104",
        "ENST00000631131"),
}
LEFT_OUT["gene"] = tuple(sorted(LEFT_OUT["ref"] + (
    "ENST00000351205", "ENST00000380707", "ENST00000503079", "ENST00000506163",
    "ENST00000506239", "ENST00000507905", "ENST00000510679", "ENST00000513228",
    "ENST00000514951", "ENST00000518504", "ENST00000625245", "ENST00000870900",
    "ENST00000932548", "ENST00000932549", "ENST00000932550", "ENST00000932551",
    "ENST00000932552", "ENST00000932553", "ENST00000932554", "ENST00000932555",
    "ENST00000979241", "ENST00001022050", "ENST00001133339")))
#: SHA-256 of each arm's left-out ids, sorted, one per line
LEFT_OUT_SHA256 = {
    "ref": "e1961e0f15a0f7fd59b1c2a21e8cb1395f7880bb29c52f1cfbf50a0d91e5b434",
    "gene": "d62aa9718e9e8644382fe27932a8c1a696c3264bb88eb9642c3c4be522c858fe",
}
#: SHA-256 of each arm's FASTA text (decompressed): cdna.all's, less the left-out records
ARM_FASTA_SHA256 = {
    "ref": "697899f6d73e1125808be11bd4ca873f13b185283287c6ac328b82c143d46075",
    "gene": "687ae4df6d456a5574951e08b3b68573f329e7761ac429e8125d62f7e3d79572",
}
#: the records of each arm's index: the stored index's 453,553, less the left-out records it
#: kept (a left-out record it dropped for an identical copy counts nothing, and one it kept
#: in place of a record that stays is replaced by it)
NUM_REFS = {"ref": 453527, "gene": 453507}
INDEX_INFO = {"salmon_version": "2.7.0", "k": 31, "keep_duplicates": False, "num_decoys": 0}


class NotGradeable(Exception):
    """The runs, the index, the FASTA or the predictions are not the ones D1 names."""


def _gate():
    spec = importlib.util.spec_from_file_location("_gate_3_0",
                                                  os.path.join(HERE, "gate_3_0.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _ids_sha256(ids):
    return hashlib.sha256(("\n".join(sorted(ids)) + "\n").encode("utf-8")).hexdigest()


def leaves_out(header, arm):
    """Does ``arm`` leave out the record of this FASTA header (without its ``>``)?"""
    parts = header.split()
    on_5 = len(parts) > 2 and parts[2].startswith("chromosome:GRCh38:5:")
    m = re.search(r"\sgene_symbol:(\S+)", " " + header)
    symbol = m.group(1) if m else None
    if arm == "ref":
        return symbol in ("SMN1", "SMN2") and not on_5
    if arm == "gene":
        return symbol == "SMN1" or (symbol == "SMN2" and not on_5)
    raise ValueError("no arm %r" % (arm,))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_fasta(cdna, arm, out):
    """Write ``arm``'s FASTA, gzipped, to ``out`` and return the ids it left out; refuse an
    input, a left-out set or an output text that is not the one fixed here."""
    got = _sha256(cdna)
    if got != CDNA_SHA256:
        raise NotGradeable("%s has SHA-256 %s, not cdna.all's %s" % (cdna, got, CDNA_SHA256))
    left, keep, text = [], True, hashlib.sha256()
    tmp = out + ".part"
    with gzip.open(cdna, "rt", encoding="utf-8", newline="") as f, \
            gzip.GzipFile(tmp, "wb", compresslevel=1, mtime=0) as o:
        for line in f:
            if line.startswith(">"):
                keep = not leaves_out(line[1:].rstrip("\n"), arm)
                if not keep:
                    left.append(line[1:].split()[0].split(".")[0])
            if keep:
                b = line.encode("utf-8")
                o.write(b)
                text.update(b)
    try:
        if _ids_sha256(left) != LEFT_OUT_SHA256[arm] or sorted(left) != sorted(LEFT_OUT[arm]):
            raise NotGradeable("%s left out %d records, not the %d fixed for it"
                               % (arm, len(left), len(LEFT_OUT[arm])))
        if text.hexdigest() != ARM_FASTA_SHA256[arm]:
            raise NotGradeable("%s's FASTA text has SHA-256 %s, not %s"
                               % (arm, text.hexdigest(), ARM_FASTA_SHA256[arm]))
    except NotGradeable:
        os.remove(tmp)
        raise
    os.replace(tmp, out)
    return left


def check_index(arm, index):
    """Refuse an index that is not one Salmon 2.7.0 built from ``arm``'s FASTA as the stored
    one was built."""
    with open(os.path.join(index, "info.json"), encoding="utf-8") as f:
        info = json.load(f)
    for key, want in INDEX_INFO.items():
        if info.get(key) != want:
            raise NotGradeable("%s's index has %s %r, not %r" % (arm, key, info.get(key), want))
    if info.get("num_refs") != NUM_REFS[arm]:
        raise NotGradeable("%s's index has %r records, not %d"
                           % (arm, info.get("num_refs"), NUM_REFS[arm]))
    return info


def check_quant(arm, path):
    """Refuse a ``quant.sf`` that does not list :data:`NUM_REFS` records, or lists a record
    ``arm`` leaves out."""
    opener = gzip.open if path.endswith(".gz") else open
    names = []
    with opener(path, "rt", encoding="utf-8") as f:
        next(f)
        for line in f:
            names.append(line.split("\t", 1)[0].split(".")[0])
    there = sorted(set(names) & set(LEFT_OUT[arm]))
    if there:
        raise NotGradeable("%s lists %d of the records %s leaves out (%s)"
                           % (path, len(there), arm, ", ".join(there[:3])))
    if len(names) != NUM_REFS[arm]:
        raise NotGradeable("%s lists %d records, not %d" % (path, len(names), NUM_REFS[arm]))
    return len(names)


def _row(run, key="genes_tpm", gene=GENE):
    return next((r for r in run.get(key) or () if r["gene"] == gene), None)


def meets(row, se):
    """Step 3's test of one row: |bias| and replicate SD both below ``se``; no row fails."""
    return row is not None and abs(row["bias"]) < se and row["realised_sd"] < se


def grade(runs, predictions, stored):
    """``runs`` is ``{arm: run}``.  Returns the predicted SE, each arm's rows and grade, the
    stored run's, the other genes' largest change, and the exit status."""
    gate = _gate()
    for arm in ARMS:
        try:
            gate.check_run(runs[arm], predictions)
        except gate.NotGradeable as e:
            raise NotGradeable("%s: %s" % (arm, e)) from None
        if runs[arm].get("design", {}).get("genes_simulated") != \
                stored["design"]["genes_simulated"]:
            raise NotGradeable("%s's libraries did not hold the stored run's genes" % arm)
    se = predictions["primary"][GENE]
    out = {"predicted_se": se, "arms": {}, "stored": {}}
    for name, run in [("stored", stored)] + [(a, runs[a]) for a in ARMS]:
        tpm, num = _row(run), _row(run, "genes")
        entry = {"tpm": tpm, "num_reads": num, "meets": meets(tpm, se)}
        if name == "stored":
            out["stored"] = entry
        else:
            out["arms"][name] = entry
    others = {}
    for arm in ARMS:
        want = {r["gene"]: r for r in stored["genes_tpm"] if r["gene"] != GENE}
        got = {r["gene"]: r for r in runs[arm]["genes_tpm"] if r["gene"] != GENE}
        worst = {"bias": (0.0, None), "realised_sd": (0.0, None)}
        for g in sorted(set(want) | set(got)):
            for key in worst:
                if g not in want or g not in got:
                    worst[key] = (float("inf"), g)
                    continue
                d = abs(got[g][key] - want[g][key])
                if d > worst[key][0]:
                    worst[key] = (d, g)
        others[arm] = worst
    out["others"] = others
    out["exit"] = 0 if out["arms"]["gene"]["meets"] else 1
    return out


def _fmt(row):
    if row is None:
        return "no row (fewer than 3 replicates count)"
    return "|bias| %.3f, SD %.3f, n %d" % (abs(row["bias"]), row["realised_sd"], row["n"])


def report(result):
    se = result["predicted_se"]
    lines = ["D1 (validation/README.md, \"What produced SMN2's failure (D1)\")",
             "SMN2's predicted SE, fixed with the gate: %.3f" % se,
             "summed TPM, the row graded; NumReads beside it, deciding nothing:"]
    for name, e in [("stored run (all records)", result["stored"])] + [
            ("%s (%d records left out)" % (a, len(LEFT_OUT[a])), result["arms"][a])
            for a in ARMS]:
        lines.append("  %-30s TPM %s: %s; NumReads %s" % (
            name, _fmt(e["tpm"]), "meets it" if e["meets"] else "does not meet it",
            _fmt(e["num_reads"])))
    for arm in ARMS:
        w = result["others"][arm]
        lines.append("  the other 48 genes' summed-TPM rows in %s against the stored run: "
                     "largest bias change %.4f (%s), largest SD change %.4f (%s)"
                     % (arm, w["bias"][0], w["bias"][1] or "none changes",
                        w["realised_sd"][0], w["realised_sd"][1] or "none changes"))
    ref = "met" if result["arms"]["ref"]["meets"] else "did not meet"
    if result["exit"] == 0:
        lines.append("D1 exit 0: THE INDEX. With SMN1 and SMN2 held to SMN2's own transcripts, "
                     "SMN2's summed TPM met its predicted SE; with the README's index scope "
                     "(ref) it %s it." % ref)
    else:
        lines.append("D1 exit 1: NOT THE INDEX. With SMN1 and SMN2 held to SMN2's own "
                      "transcripts, SMN2's summed TPM still did not meet its predicted SE "
                      "(with ref it %s it); what produced the failure stays unmeasured." % ref)
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fasta")
    p.add_argument("cdna")
    p.add_argument("arm", choices=ARMS)
    p.add_argument("out")
    p = sub.add_parser("index")
    p.add_argument("arm", choices=ARMS)
    p.add_argument("index")
    p = sub.add_parser("quant")
    p.add_argument("arm", choices=ARMS)
    p.add_argument("quant")
    p = sub.add_parser("grade")
    p.add_argument("run_ref")
    p.add_argument("run_gene")
    p.add_argument("--stored", default=STORED_TPM_RUN)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "fasta":
            left = write_fasta(a.cdna, a.arm, a.out)
            print("%s: wrote %s, %d records left out (SHA-256 of their ids %s)"
                  % (a.arm, a.out, len(left), LEFT_OUT_SHA256[a.arm]))
        elif a.cmd == "index":
            info = check_index(a.arm, a.index)
            print("%s: index %s, Salmon %s, %d records, seq_hash %s"
                  % (a.arm, a.index, info["salmon_version"], info["num_refs"],
                     info.get("seq_hash")))
        elif a.cmd == "quant":
            n = check_quant(a.arm, a.quant)
            print("%s: %s lists %d records, none left out" % (a.arm, a.quant, n))
        else:
            gate = _gate()
            try:
                predictions = gate.frozen_predictions()
            except gate.NotGradeable as e:
                raise NotGradeable(str(e)) from None
            runs = {}
            for arm, path in (("ref", a.run_ref), ("gene", a.run_gene)):
                with open(path, encoding="utf-8") as f:
                    runs[arm] = json.load(f)
            with open(a.stored, encoding="utf-8") as f:
                stored = json.load(f)
            result = grade(runs, predictions, stored)
            print(report(result))
            return result["exit"]
    except NotGradeable as e:
        if a.cmd == "grade":
            print("D1 exit 2: NOT GRADED. %s" % e)
        else:
            print("d1_smn2.py %s refused: %s" % (a.cmd, e), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
