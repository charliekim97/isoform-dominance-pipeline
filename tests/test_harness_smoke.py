"""The harness end to end without Salmon, and the simulator across processes.

CI runs no quantifier, so ``run_simulation.py run`` is driven here with only the ``salmon``
command intercepted: the stub writes the files Salmon writes -- ``quant.sf`` with each
transcript's ``NumReads`` the number of fragments the simulator drew from it, read back from
the FASTQ, an ``EffectiveLength`` and the ``TPM`` those give, ``aux_info/meta_info.json`` and
``lib_format_counts.json`` -- and everything else is the harness's own: the simulator as a
subprocess, ``_quantify_salmon``'s reading of those files, ``summarise`` and
``summarise_tpm``.  A quantifier that returns the truth has no error, so every row's bias and
RMSE is 0, the TPM rows' too (to the six places ``quant.sf`` gives TPM).

The harness runs the simulator as a subprocess, so an order that depends on Python's hash
seed would show only between processes; ``test_validation_harness`` compares two runs in
one process.  Two subprocesses at ``PYTHONHASHSEED`` 0 and 1 write the same truth and the
same reads (the FASTQ decompressed: gzip's header holds the time).  Both comparisons are
within one environment; nothing is claimed across numpy versions.
"""
import collections
import gzip
import hashlib
import importlib.util
import json
import os
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
BUNDLE = VALIDATION / "data" / "sweep116.json.gz"

#: two genes whose smaller class expects >= 50 fragments at DEPTH (129 and 59), so that no
#: replicate loses a class to a Poisson zero under any numpy's stream
GENES = ("LEPR", "TGFB1")
DEPTH = 0.04
REPLICATES = 3              # summarise drops a gene with fewer usable replicates


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sim = _module("simulate_reads")
harness = _module("run_simulation")


@pytest.fixture(scope="module")
def two_genes(tmp_path_factory):
    """The bundle of GENES alone, in the bundle's order."""
    genes = sim.read_inputs(inputs=str(BUNDLE))
    path = tmp_path_factory.mktemp("bundle") / "two.json.gz"
    with gzip.open(path, "wt") as f:
        json.dump({"genes": {g: genes[g] for g in GENES}}, f)
    return path


def test_the_two_genes_expect_enough_fragments_in_the_smaller_class(two_genes):
    for gene, g in sim.load_genes(inputs=str(two_genes)).items():
        e = sim.transcript_fragment_counts(g, 10.0, DEPTH, None, None, expectations=True)
        assert min(sum(e[t] for t in g["groups"][c]) for c in g["pc"]) >= 50, gene


def _fragments(r1):
    """Fragments per transcript in a read-1 FASTQ, from the read names
    (``<gene>_r<rep>:<transcript>:<i>/1``)."""
    n = collections.Counter()
    with gzip.open(r1, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                n[line[1:].split(":")[1]] += 1
    return n


class SalmonStub:
    """``salmon quant`` answered with the truth; every other command run for real.

    ``bump`` adds that many fragments to one transcript of each gene's first compared
    class, to show that a quantifier off by one fragment is seen."""

    def __init__(self, sources, bump=0):
        self.sources, self.bump, self.calls = sources, bump, 0
        self.real = harness.sh

    def __call__(self, cmd, **kw):
        if cmd[0] != "salmon":
            return self.real(cmd, **kw)
        assert cmd[1] == "quant" and cmd[cmd.index("-l") + 1] == "A", cmd
        self.calls += 1
        r1, out = cmd[cmd.index("-1") + 1], cmd[cmd.index("-o") + 1]
        n = _fragments(r1)
        for entry in self.sources.values():
            first = entry["config"]["groups"][entry["config"]["primary_comparison"][0]][0]
            n[first.split(".")[0]] += self.bump
        os.makedirs(os.path.join(out, "aux_info"))
        # an effective length that is not the simulator's max(1, L - 199), as Salmon's is not:
        # the TPM truth must take the quantifier's own
        listed = [(t, s) for entry in self.sources.values() for t, s in entry["sequences"].items()]
        eff = {t: max(1.0, len(s) - 180.0) for t, s in listed}
        rate = {t: n[t.split(".")[0]] / eff[t] for t, _ in listed}
        scale = 1e6 / sum(rate.values())
        with open(os.path.join(out, "quant.sf"), "w") as f:
            f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
            for t, s in listed:
                f.write("%s\t%d\t%.3f\t%.6f\t%.3f\n" % (t, len(s), eff[t], rate[t] * scale,
                                                       n[t.split(".")[0]]))
        total = sum(_fragments(r1).values())
        with open(os.path.join(out, "aux_info", "meta_info.json"), "w") as f:
            json.dump({"salmon_version": "stub", "percent_mapped": 100.0,
                       "num_processed": total, "num_mapped": total,
                       "frag_length_mean": 200.0, "frag_length_sd": 60.0}, f)
        with open(os.path.join(out, "lib_format_counts.json"), "w") as f:
            json.dump({"expected_format": "ISF"}, f)
        return subprocess.CompletedProcess(cmd, 0)


def _run(tmp_path, monkeypatch, two_genes, bump=0):
    stub = SalmonStub(sim.read_inputs(inputs=str(two_genes)), bump)
    monkeypatch.setattr(harness, "sh", stub)
    out = tmp_path / "run.json"
    with pytest.raises(SystemExit) as stop:
        harness.main(["run", "--inputs", str(two_genes), "--index", "no_index",
                      "--work", str(tmp_path / "work"), "--depth-scale", str(DEPTH),
                      "--replicates", str(REPLICATES), "--threads", "1", "--out", str(out)])
    assert stop.value.code == 0
    assert stub.calls == REPLICATES
    return json.loads(out.read_text())


def test_a_quantifier_that_returns_the_truth_has_no_error(tmp_path, monkeypatch, two_genes):
    res = _run(tmp_path, monkeypatch, two_genes)
    rows = res["genes"]
    assert sorted(r["gene"] for r in rows) == sorted(GENES)
    for r in rows:
        assert r["n"] == REPLICATES, r
        assert (r["bias"], r["rmse"]) == (0.0, 0.0), r
        assert r["realised_sd"] > 0                    # the truth itself varies
    # the TPM rows: no error but the six places of quant.sf, and the abundance prediction
    tpm = res["genes_tpm"]
    assert sorted(r["gene"] for r in tpm) == sorted(GENES)
    for r in tpm:
        assert r["n"] == REPLICATES, r
        assert abs(r["bias"]) < 1e-6 and r["rmse"] < 1e-6, r
        assert r["realised_sd"] > 0 and r["predicted_se"] > 0
    assert res["effective_length_taken_from"] == {}
    assert res["design"]["genes_simulated"] == sorted(GENES)
    assert res["design"]["inputs_sha256"] == hashlib.sha256(two_genes.read_bytes()).hexdigest()
    assert res["salmon"] == {"version": "stub", "kallisto_version": None,
                             "percent_mapped_min": 100.0, "library_type": "ISF"}
    assert res["design"]["depth_scale"] == DEPTH
    # the replicate libraries are removed, the record of what salmon did is kept
    work = tmp_path / "work"
    assert not list(work.glob("*.fq.gz"))
    assert sorted(p.name for p in work.glob("*.salmon.json")) == [
        "rep%d.salmon.json" % i for i in range(REPLICATES)]


def test_a_quantifier_one_fragment_off_has_an_error(tmp_path, monkeypatch, two_genes):
    res = _run(tmp_path, monkeypatch, two_genes, bump=1)
    rows, tpm = res["genes"], res["genes_tpm"]
    assert len(rows) == 2 and len(tpm) == 2
    assert all(r["bias"] > 0 and r["rmse"] > 0 for r in rows)
    assert all(r["bias"] > 1e-4 and r["rmse"] > 1e-4 for r in tpm)


# --------------------------------------------------------------------------- #
# the simulator across processes and hash seeds
# --------------------------------------------------------------------------- #
def _simulate(out, seed):
    """Replicate 0 of GENES in a subprocess at ``PYTHONHASHSEED=seed``: ``(truth, read 1,
    read 2)``, the FASTQ decompressed."""
    env = dict(os.environ, PYTHONHASHSEED=str(seed))
    subprocess.run([sys.executable, str(VALIDATION / "simulate_reads.py"), str(out),
                    "--inputs", str(BUNDLE), "--genes", ",".join(GENES),
                    "--depth-scale", "0.01", "--replicate", "0", "--bias-5p=1.0"],
                   check=True, env=env, capture_output=True)
    truth = json.loads((out / "rep0.truth.json").read_text())
    reads = [gzip.decompress((out / ("rep0_%d.fq.gz" % i)).read_bytes()) for i in (1, 2)]
    return truth, reads[0], reads[1]


def test_the_simulator_does_not_depend_on_the_hash_seed(tmp_path):
    a = _simulate(tmp_path / "seed0", 0)
    b = _simulate(tmp_path / "seed1", 1)
    assert set(a[0]["genes"]) == set(GENES)
    assert all(a[0]["genes"][g]["fragments_written"] > 0 for g in GENES)
    assert a[1].count(b"\n") == 4 * a[0]["fragments_written"]
    assert a == b
