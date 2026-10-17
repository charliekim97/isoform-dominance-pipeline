"""The harness reads the stored inputs bundle as it read the per-gene files it was made from.

The simulator draws each gene's transcripts in the order its input gives them, so the bundle
must give them in the files' order: a replicate simulated from the bundle and one simulated
from the bundle unpacked into per-gene files, at the same seed, are the same fragments.  The
bundle's genes and each gene's transcripts happen to be in sorted order, so on it keeping the
order and sorting give the same draws; the order is tested on a reordered copy.
"""
import gzip
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
BUNDLE = VALIDATION / "data" / "sweep116.json.gz"

#: the two smallest genes of the bundle by sequence, at a small depth: a few fragments each
GENES = ("MCL1", "SMN2")
DESIGN = ["--depth-scale", "0.05", "--replicate", "3", "--bias-5p=1.0"]


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sim = _module("simulate_reads")
harness = _module("run_simulation")


def _bundle():
    with gzip.open(BUNDLE, "rt") as f:
        return json.load(f)


def _unpack(genes, rundir):
    """The bundle's genes as per-gene files, each in the bundle's order."""
    rundir.mkdir()
    for gene, g in genes.items():
        (rundir / ("%s.cfg.json" % gene)).write_text(json.dumps(g["config"], indent=1))
        (rundir / ("%s.seqs.json" % gene)).write_text(json.dumps(g["sequences"]))


def _replicate(out, genes, rundir=None, inputs=None):
    """``(truth, read 1, read 2)`` of replicate 3, the FASTQ decompressed."""
    where = [str(rundir), str(out)] if rundir else [str(out), "--inputs", str(inputs)]
    sim.main(where + ["--genes", ",".join(genes), *DESIGN])
    truth = json.loads((out / "rep3.truth.json").read_text())
    reads = [gzip.decompress((out / ("rep3_%d.fq.gz" % i)).read_bytes()) for i in (1, 2)]
    return truth, reads[0], reads[1]


def test_the_bundle_and_its_unpacked_files_simulate_the_same_replicate(tmp_path):
    _unpack(_bundle()["genes"], tmp_path / "sweep")
    a = _replicate(tmp_path / "from_bundle", GENES, inputs=BUNDLE)
    b = _replicate(tmp_path / "from_files", GENES, rundir=tmp_path / "sweep")
    assert set(a[0]["genes"]) == set(GENES)
    assert all(a[0]["genes"][g]["fragments_written"] > 0 for g in GENES)
    assert a == b
    # and the prediction the harness compares against is the same from either
    src = {"bundle": sim.read_inputs(inputs=BUNDLE),
           "files": sim.read_inputs(rundir=str(tmp_path / "sweep"))}
    assert src["bundle"] == src["files"]
    for gene in GENES:
        assert (harness.predicted(src["bundle"], gene, 10.0, 1.0, 1500.0)
                == harness.predicted(src["files"], gene, 10.0, 1.0, 1500.0))


def test_each_gene_is_drawn_in_the_order_its_input_gives(tmp_path):
    genes = _bundle()["genes"]
    gene = GENES[0]
    turned = dict(genes[gene], sequences=dict(reversed(list(genes[gene]["sequences"].items()))))
    assert list(turned["sequences"]) != sorted(turned["sequences"])
    reordered = tmp_path / "reordered.json.gz"
    with gzip.open(reordered, "wt") as f:
        json.dump({"genes": {gene: turned}}, f)
    assert list(sim.read_inputs(inputs=reordered)[gene]["sequences"]) == list(turned["sequences"])
    _unpack({gene: turned}, tmp_path / "sweep")
    a = _replicate(tmp_path / "a", [gene], inputs=reordered)
    assert a == _replicate(tmp_path / "b", [gene], rundir=tmp_path / "sweep")
    # the order matters: the bundle's own order draws other fragments
    assert a[1] != _replicate(tmp_path / "c", [gene], inputs=BUNDLE)[1]
