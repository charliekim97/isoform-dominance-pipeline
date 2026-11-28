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
import math
import pathlib
import random

import pytest

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


def _replicate(out, genes, rundir=None, inputs=None, extra=()):
    """``(truth, read 1, read 2)`` of replicate 3, the FASTQ decompressed."""
    where = [str(rundir), str(out)] if rundir else [str(out), "--inputs", str(inputs)]
    sim.main(where + ["--genes", ",".join(genes), *DESIGN, *extra])
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


@pytest.mark.parametrize("target", [0.25, 1.0, 4.0])
def test_abundance_ratio_preserves_expected_fragments(target):
    gene = {"seqs": {"a1": "A" * 1000, "a2": "C" * 1600, "b": "G" * 3400,
                     "bg": "T" * 2100},
            "groups": {"A": ["a1", "a2"], "B": ["b"]}, "pc": ["A", "B"]}
    base = sim.transcript_fragment_counts(gene, 10.0, 1.0, None, None, expectations=True)
    counts = sim.transcript_fragment_counts(
        gene, 10.0, 1.0, None, None, expectations=True, abundance_ratio=target)
    abundance = {t: counts[t] / (len(seq) - 199) for t, seq in gene["seqs"].items()}
    assert (abundance["a1"] + abundance["a2"]) / abundance["b"] == pytest.approx(target)
    assert abundance["a1"] == pytest.approx(abundance["a2"])
    assert sum(counts.values()) == pytest.approx(sum(base.values()))
    assert counts["bg"] == base["bg"]
    assert list(counts) == list(base)

    legacy = sim.transcript_fragment_counts(gene, 10.0, 1.0, target, None, expectations=True)
    assert (legacy["a1"] + legacy["a2"]) / legacy["b"] == pytest.approx(target)


def test_abundance_target_is_recorded_without_changing_fragment_truth(tmp_path):
    gene = GENES[0]
    truth, _, _ = _replicate(tmp_path / "abundance", [gene], inputs=BUNDLE,
                            extra=["--abundance-ratio", "4"])
    assert truth["abundance_ratio"] == 4.0
    assert truth["ratio"] is None
    row = truth["genes"][gene]
    a, b = row["primary_comparison"]
    assert row["true_log2_ratio"] == pytest.approx(
        math.log2(row["class_fragments"][a] / row["class_fragments"][b]))


@pytest.mark.parametrize("flags", [
    ["--abundance-ratio", "0"],
    ["--abundance-ratio", "nan"],
    ["--abundance-ratio", "4", "--ratio", "4"],
])
def test_invalid_abundance_target_is_rejected_before_writing(tmp_path, flags):
    output = tmp_path / "invalid"
    with pytest.raises(SystemExit) as error:
        sim.main([str(output), "--inputs", str(BUNDLE), *flags])
    assert error.value.code == 2
    assert not output.exists()


def test_invalid_abundance_rerun_keeps_existing_fastqs_and_truth(tmp_path, capsys):
    output = tmp_path / "replicate"
    flags = ["--abundance-ratio", "4"]
    _replicate(output, [GENES[0]], inputs=BUNDLE, extra=flags)
    saved = {p: p.read_bytes() for p in output.iterdir()}
    assert len(saved) == 3

    with pytest.raises(SystemExit) as error:
        sim.main([str(output), "--inputs", str(BUNDLE), "--genes", GENES[0],
                  *DESIGN, *flags, "--tpm", "0"])
    assert error.value.code == 2
    assert "positive expected fragments" in capsys.readouterr().err
    assert {p: p.read_bytes() for p in output.iterdir()} == saved


@pytest.mark.parametrize("comparison", [None, ["A", "B"]])
def test_matched_abundance_se_uses_window_weights_at_the_simulated_mix(comparison):
    # Independent sequences give a diagonal observation model. Class A has two
    # different lengths, so changing from window to fragment lengths changes its SE.
    randomizer = random.Random(878)
    lengths = {"a1": 250, "a2": 550, "b": 1000}
    sequences = {t: "".join(randomizer.choices("ACGT", k=length))
                 for t, length in lengths.items()}
    config = {"groups": {"B": ["b"], "A": ["a1", "a2"]}}
    if comparison is not None:
        config["primary_comparison"] = comparison
    source = {"gene": {"config": config, "sequences": sequences}}
    got = harness._predictions_at_abundance_ratio(source, 10.0, 0.1, 1500.0, 4.0)

    effective = {t: length - 199 for t, length in lengths.items()}
    total = (3e6 * 10.0 / (1e6 * 1500.0)) * sum(effective.values())
    # Two transcripts at relative abundance 2 each versus one at 1.
    mass = {"a1": 2 * effective["a1"], "a2": 2 * effective["a2"], "b": effective["b"]}
    expected = {t: total * value / sum(mass.values()) for t, value in mass.items()}
    for name in ("fragment", "abundance"):
        denom = {t: length - 30 if name == "abundance" else 1 for t, length in lengths.items()}
        variance = sum(
            sum(expected[t] / denom[t] ** 2 for t in group)
            / sum(expected[t] / denom[t] for t in group) ** 2
            for group in config["groups"].values())
        assert got[name]["gene"]["predicted_se"] == pytest.approx(
            math.sqrt(variance) / math.log(2.0))


def test_matched_prediction_rejects_a_transcript_that_cannot_emit_reads():
    sequences = {"a": "A" * 200, "b": "C" * 99}
    config = {"groups": {"A": ["a"], "B": ["b"]}}
    source = {"gene": {"config": config, "sequences": sequences}}
    with pytest.raises(ValueError, match="cannot emit reads"):
        harness._predictions_at_abundance_ratio(source, 10.0, 0.1, 1500.0, 4.0)


def test_reference_library_shares_competitors_and_uses_their_own_expected_mix(tmp_path, monkeypatch):
    randomizer = random.Random(491)
    sequences = {t: "".join(randomizer.choices("ACGT", k=210))
                 for t in ("a", "b", "c", "d", "extra")}
    sequences["c"] = sequences["a"][:50] + sequences["c"][50:]
    sequences["extra"] = sequences["a"][:45] + sequences["d"][:45] + sequences["extra"][90:]
    panel = {}
    for gene, pair in (("G1", ("a", "b")), ("G2", ("c", "d"))):
        panel[gene] = {"config": {"groups": {"A": [pair[0]], "B": [pair[1]]}},
                       "sequences": {t: sequences[t] for t in pair}}
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{t}\n{seq}\n" for t, seq in sequences.items()))
    backgrounds, competitors = harness._reference_background(panel, str(fasta))
    assert set(competitors) == {"extra"}
    assert set(backgrounds["G1"]["added"]) == {"c", "extra"}
    assert set(backgrounds["G2"]["added"]) == {"a", "extra"}

    seen_rates = []
    covariance = harness.I.gls_covariance

    def capture(matrix, theta):
        seen_rates.append(theta.copy())
        return covariance(matrix, theta)

    monkeypatch.setattr(harness.I, "gls_covariance", capture)
    harness._predictions_at_abundance_ratio(
        panel, 10.0, 0.1, 1500.0, 4.0, backgrounds=backgrounds, competitors=competitors)
    rates = sim.transcript_fragment_counts(
        sim.simulation_gene(panel["G2"]), 10.0, 0.1, None, None, 1500.0,
        expectations=True, abundance_ratio=4.0)
    ids = [*panel["G1"]["sequences"], *backgrounds["G1"]["added"]]
    assert seen_rates[0][ids.index("c")] == pytest.approx(rates["c"])


def test_recorded_reference_refuses_changed_fasta_and_index(tmp_path):
    index = tmp_path / "index"
    index.mkdir()
    metadata = {"k": 31, "canonical": True, "num_decoys": 0, "keep_duplicates": False,
                "seq_hash": "seq", "name_hash": "names", "salmon_version": "2.7.0"}
    info = index / "info.json"
    info.write_text(json.dumps(metadata))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">a\nACGT\n")
    manifest = {"fasta_sha256": harness.io.file_sha256(fasta),
                "index": harness._index_identity(index)}
    (index / harness.REFERENCE_MANIFEST).write_text(json.dumps(manifest))
    assert harness._reference_identity(index, fasta) == manifest

    fasta.write_text(">a\nTCGT\n")
    with pytest.raises(ValueError, match="FASTA hash differs"):
        harness._reference_identity(index, fasta)
    fasta.write_text(">a\nACGT\n")
    info.write_text(json.dumps(dict(metadata, name_hash="other names")))
    with pytest.raises(ValueError, match="index hashes differ"):
        harness._reference_identity(index, fasta)
