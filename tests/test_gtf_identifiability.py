"""``identifiability --gtf --transcripts-fasta``: the run REST makes, from local files (issue #13).

The gene's transcripts come from the GTF and their sequence from the transcript FASTA of
the same release.  That background then goes the way REST's does -- the identical-copy
rule, the id merge with ``--background-fasta``, ``--max-window-records`` and ``--decoys``
-- so that a background the two sources give alike gives one report.
"""
import gzip
import json
import pathlib
import urllib.request

import pytest

from isoform_dominance import annotate, cli, io
from isoform_dominance import annotation_files as af
from isoform_dominance import identifiability as I

DATA = pathlib.Path(__file__).parent / "data" / "gencode_mini"
GTF = str(DATA / "gencode.v50.mini.gtf.gz")
FASTA = str(DATA / "gencode.v50.mini.transcripts.fa.gz")
REST = json.loads(gzip.decompress((DATA / "rest116_mini.json.gz").read_bytes()))
SERVER = "https://rest.ensembl.org"
LEPR, CD99_X, CD99_Y = "ENSG00000116678", "ENSG00000002586", "ENSG00000292348"


def _records(path=FASTA):
    """``[(header, sequence)]`` of a FASTA, in file order."""
    out = []
    with io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append([line[1:].rstrip("\n"), []])
            else:
                out[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in out]


RECORDS = _records()
BY_ID = {h.split("|")[0].split(".")[0]: (h.split("|")[0], s) for h, s in RECORDS}


def _write_fasta(path, records):
    path.write_text("".join(">%s\n%s\n" % (h, s) for h, s in records))
    return str(path)


class Rest116:
    """REST at release 116 for the genes of the extract: the recorded lookups, and each
    transcript's cDNA as the FASTA holds it (the parity tests hold the two equal)."""

    def __init__(self):
        self.calls = []

    def __call__(self, req, timeout=None, **_):
        from test_ensembl_http import _Resp
        path = req.full_url[len(SERVER):]
        body = json.loads(req.data) if req.data else None
        self.calls.append((req.get_method(), path))
        parts = path.split("?")[0].split("/")
        if path.startswith("/info/data"):
            return _Resp(json.dumps({"releases": [116]}))
        if path.startswith("/lookup/id/"):
            return _Resp(json.dumps(REST["lookups"][parts[3].split(".")[0]]))
        if path.startswith("/lookup/symbol/"):
            sym = parts[4]
            gid = REST["lookup_symbol"].get(sym) or next(
                g for g, r in REST["lookups"].items() if r["display_name"] == sym)
            return _Resp(json.dumps(REST["lookups"][gid]))
        if path.startswith("/sequence/id") and req.get_method() == "POST":
            return _Resp(json.dumps([
                {"query": i, "id": BY_ID[i][0], "version": int(BY_ID[i][0].split(".")[1]),
                 "seq": BY_ID[i][1]} for i in body["ids"] if i in BY_ID]))
        raise AssertionError("unexpected request %s" % path)


@pytest.fixture
def rest(monkeypatch):
    fake = Rest116()
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    return fake


@pytest.fixture
def offline(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("a request was made")
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _cfg(gene="LEPR", **kw):
    return annotate.build_config_from_gtf(gene, GTF, **kw)


def _without_provenance(res):
    """A report less what says where its sequence came from, and less what only a GTF can
    say -- where a same-name gene of the index lies -- as JSON: a figure that is not a
    number (NaN) is equal to itself there."""
    bg = dict(res["background"], same_name_copies=[
        {k: c[k] for k in ("gene_id", "gene_name", "transcripts")}
        for c in res["background"]["same_name_copies"]])
    return json.dumps(dict(res, annotation=None, background=bg), sort_keys=True, default=str)


# --------------------------------------------------------------------------- #
# one background, one report
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("gene", ["LEPR", "CD99"])
@pytest.mark.parametrize("extra", [
    {},
    {"background_fasta": FASTA},
    {"background_fasta": FASTA, "max_window_records": 1},
    {"background_fasta": FASTA, "keep_duplicates": True},
    {"keep_duplicates": True},
], ids=["gene", "fasta", "fasta-M1", "fasta-keep", "keep"])
def test_the_files_give_the_report_rest_gives(rest, gene, extra):
    cfg = _cfg(gene)
    by_rest = I.analyze(cfg, min_log2fc=0.5, **extra)
    n = len(rest.calls)
    by_files = I.analyze(cfg, gtf=GTF, transcripts_fasta=FASTA, min_log2fc=0.5, **extra)
    assert len(rest.calls) == n                         # the files: no request
    assert _without_provenance(by_files) == _without_provenance(by_rest)
    # and the background is there to be equal: the gene's other transcripts, by id
    assert by_files["background"]["gene_id"] == cfg["gene_id"]
    assert by_files["background"]["n_background_transcripts"] > 5


def test_a_copy_in_the_gene_background_is_labelled_as_rests_is(rest):
    """CD99: the chrX gene's background, with the FASTA as the index -- the identical-copy
    rule (issue #15) and its source label, the id merge, the M rule."""
    cfg = _cfg("CD99")
    by_files = I.analyze(cfg, gtf=GTF, transcripts_fasta=FASTA, background_fasta=FASTA)
    by_rest = I.analyze(cfg, background_fasta=FASTA)
    bg = by_files["background"]
    assert bg["identical_source"] == by_rest["background"]["identical_source"]
    assert set(bg["identical_source"].values()) <= {"gene", "fasta"}


@pytest.fixture
def copy_fasta(tmp_path, monkeypatch):
    """The mini FASTA with a LEPR transcript of the gene background given a configured
    transcript's sequence -- a copy in the gene background -- and REST serving the same."""
    cfg = _cfg()
    configured = {t for ids in cfg["groups"].values() for t in ids}
    a, b = sorted(configured)[0], sorted(_lepr_ids() - configured)[0]
    seq = BY_ID[a][1]
    monkeypatch.setitem(BY_ID, b, (BY_ID[b][0], seq))
    path = _write_fasta(tmp_path / "copy.fa", [(h, seq if h.startswith(b + ".") else s)
                                               for h, s in RECORDS])
    return path, a, b


def test_a_copy_in_the_gene_background_is_the_genes_from_either_source(rest, copy_fasta):
    """Issue #15's rule, with the source REST's run gives the copy: "gene"."""
    path, a, b = copy_fasta
    by_rest = I.analyze(_cfg())
    by_files = I.analyze(_cfg(), gtf=GTF, transcripts_fasta=path)
    for res in (by_rest, by_files):
        assert res["background"]["identical_to_configured"] == {b: a}
        assert res["background"]["identical_source"] == {b: "gene"}
    assert _without_provenance(by_files) == _without_provenance(by_rest)


def test_the_report_says_where_the_sequence_came_from(offline):
    res = I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA)
    ann = res["annotation"]
    assert (ann["ensembl_release"], ann["fetched_release"], ann["file_release"]) == (116, None,
                                                                                     116)
    assert ann["source"] == {
        "kind": "gtf",
        "gtf": dict(af.provenance(GTF), **af.header(GTF)),
        "transcripts_fasta": af.provenance(FASTA)}
    assert res["background"]["gene_id"] == LEPR


def test_a_rest_report_has_no_file(rest):
    ann = I.analyze(_cfg())["annotation"]
    assert (ann["fetched_release"], ann["file_release"], ann["source"]) == (116, None, None)


# --------------------------------------------------------------------------- #
# which gene
# --------------------------------------------------------------------------- #
def _bare(cfg, **kw):
    """A config as a hand would write it: groups, and none of annotate's records."""
    return dict({"gene": cfg["gene"], "groups": cfg["groups"],
                 "primary_comparison": cfg["primary_comparison"]}, **kw)


def test_without_gene_id_the_one_gene_that_holds_every_configured_transcript(offline):
    """CD99 names a chrX and a chrY gene; the configured transcripts are the chrX
    gene's, so that is the gene, whatever the symbol rule would say."""
    cfg = _cfg("CD99")
    res = I.analyze(_bare(cfg), gtf=GTF, transcripts_fasta=FASTA)
    assert res["background"]["gene_id"] == CD99_X
    assert _without_provenance(res) == _without_provenance(I.analyze(
        dict(_bare(cfg), gene_id=CD99_X), gtf=GTF, transcripts_fasta=FASTA))


def test_without_gene_id_and_one_holder_the_symbol_rule_is_rests(offline, tmp_path):
    """No one gene holds them all -- one transcript is unknown to the GTF -- so the gene
    is the symbol's, and the unknown transcript is named."""
    cfg = _cfg("CD99")
    bare = _bare(cfg)
    bare["groups"] = dict(bare["groups"])
    first = next(iter(bare["groups"]))
    bare["groups"][first] = bare["groups"][first] + ["ENST00000000999"]
    with pytest.raises(ValueError, match="gives gene %s no transcript ENST00000000999, "
                                         "which the config names" % CD99_X):
        I.analyze(bare, gtf=GTF, transcripts_fasta=FASTA)


def test_without_gene_id_the_holder_is_taken_where_the_symbol_rule_would_stop(offline):
    """HERC3 names two genes: by symbol, AmbiguousGene; the older gene's transcripts are
    held by that gene alone, and so it is the gene."""
    old = _cfg("HERC3", gene_id="ENSG00000138641")
    res = I.analyze(_bare(old), gtf=GTF, transcripts_fasta=FASTA)
    assert res["background"]["gene_id"] == "ENSG00000138641"
    # the chrY gene's transcripts of CD99, which the symbol rule would not take
    y = _cfg("CD99", gene_id=CD99_Y)
    assert I.analyze(_bare(y), gtf=GTF, transcripts_fasta=FASTA)["background"]["gene_id"] \
        == CD99_Y


def test_an_ambiguous_symbol_with_no_holder_stops_as_annotate_does(offline):
    herc3 = _cfg("HERC3", gene_id="ENSG00000138641")
    bare = _bare(herc3)
    # one configured transcript of each HERC3 gene: neither holds them all
    other = annotate.build_config_from_gtf("HERC3", GTF, gene_id="ENSG00000287542")
    bare["groups"] = {"A": herc3["groups"][herc3["primary_comparison"][0]][:1],
                      "B": [next(iter(other["groups"].values()))[0]]}
    bare["primary_comparison"] = ["A", "B"]
    with pytest.raises(annotate.AmbiguousGene, match="ENSG00000138641.*ENSG00000287542"):
        I.analyze(bare, gtf=GTF, transcripts_fasta=FASTA)


def test_a_gene_id_the_gtf_lacks_is_refused(offline):
    with pytest.raises(ValueError, match="has no gene ENSG00000000009, which the config's "
                                         "gene_id names"):
        I.analyze(dict(_cfg(), gene_id="ENSG00000000009"), gtf=GTF, transcripts_fasta=FASTA)


def test_a_configured_transcript_of_another_gene_is_refused(offline):
    cfg = _cfg()
    foxo1 = _cfg("FOXO1")
    cfg["groups"] = dict(cfg["groups"])
    first = cfg["primary_comparison"][0]
    alien = foxo1["groups"][foxo1["primary_comparison"][0]][0]
    cfg["groups"][first] = cfg["groups"][first] + [alien]
    with pytest.raises(ValueError, match="gives gene %s no transcript %s, which the config "
                                         "names" % (LEPR, alien)):
        I.analyze(cfg, gtf=GTF, transcripts_fasta=FASTA)


# --------------------------------------------------------------------------- #
# the FASTA
# --------------------------------------------------------------------------- #
def _lepr_ids():
    return {t["id"] for t in af.scan(GTF, gene_id=LEPR)[0]["Transcript"]}


def test_a_fasta_short_of_a_gene_transcript_is_refused_with_what_to_pass(offline, tmp_path):
    """GENCODE's basic FASTA lacks 1,917 of the 4,396 survey transcripts at release 116:
    a background short of them is another answer."""
    cfg = _cfg()
    configured = {t for ids in cfg["groups"].values() for t in ids}
    gone = sorted(_lepr_ids() - configured)[0]
    short = _write_fasta(tmp_path / "short.fa", [r for r in RECORDS
                                                 if not r[0].startswith(gone)])
    with pytest.raises(af.AnnotationFileError, match="no record for 1 of the 18 transcripts "
                                                     "the GTF gives %s \\(%s\\)" % (LEPR, gone)):
        I.analyze(cfg, gtf=GTF, transcripts_fasta=short)
    with pytest.raises(af.AnnotationFileError, match="basic or pc_transcripts"):
        I.analyze(cfg, gtf=GTF, transcripts_fasta=short)
    # without the gene background only the configured transcripts are needed
    res = I.analyze(cfg, gtf=GTF, transcripts_fasta=short, background_gene_transcripts=False)
    assert res["background"]["n_background_transcripts"] == 0


def test_a_fasta_of_another_version_is_refused(offline, tmp_path):
    t = sorted(_lepr_ids())[0]
    vid = BY_ID[t][0]
    bumped = "%s.%d" % (t, int(vid.split(".")[1]) + 1)
    other = _write_fasta(tmp_path / "other.fa", [
        (h.replace(vid + "|", bumped + "|"), s) for h, s in RECORDS])
    with pytest.raises(af.AnnotationFileError,
                       match="disagree on the version of 1 transcript\\(s\\) of %s \\(%s is "
                             "%s in the FASTA and %s in the GTF\\)" % (LEPR, t, bumped, vid)):
        I.analyze(_cfg(), gtf=GTF, transcripts_fasta=other)


def test_a_fasta_that_holds_a_transcript_twice_is_refused(offline, tmp_path):
    t = sorted(_lepr_ids())[0]
    twice = _write_fasta(tmp_path / "twice.fa", RECORDS + [r for r in RECORDS
                                                          if r[0].startswith(t)])
    with pytest.raises(af.AnnotationFileError, match="holds %s twice" % t):
        I.analyze(_cfg(), gtf=GTF, transcripts_fasta=twice)


def test_ensembl_headers_and_par_y_records_are_read(offline, tmp_path):
    """Ensembl's cDNA FASTA: the id is the first word.  GENCODE 25-43 hold the chrY copy
    of a pseudoautosomal transcript as ``<id>_PAR_Y``, which is skipped, not a duplicate."""
    ens = []
    for h, s in RECORDS:
        f = h.split("|")
        ens.append(("%s cdna chromosome:GRCh38:1:1:2:1 gene:%s gene_biotype:%s "
                    "transcript_biotype:%s gene_symbol:%s" % (f[0], f[1], f[7], f[7], f[5]), s))
    t = sorted(_lepr_ids())[0]
    ens += [(BY_ID[t][0] + "_PAR_Y cdna chromosome:GRCh38:Y:1:2:1", "ACGT" * 50)]
    fa = _write_fasta(tmp_path / "ensembl.fa", ens)
    assert _without_provenance(I.analyze(_cfg(), gtf=GTF, transcripts_fasta=fa)) \
        == _without_provenance(I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA))


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def _cli(tmp_path, capsys, *argv, cfg=None):
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(cfg if cfg is not None else _cfg()))
    rc = cli.main(["identifiability", "--config", str(path), "--json", *argv])
    out, err = capsys.readouterr()
    return rc, (json.loads(out) if out.strip() else None), err


@pytest.mark.parametrize("argv, says", [
    (["--gtf", GTF], "--gtf gives no sequence: pass --transcripts-fasta with it"),
    (["--transcripts-fasta", FASTA], "--transcripts-fasta without --gtf"),
    (["--gtf", GTF, "--transcripts-fasta", FASTA, "--sequences", "x.json"],
     "--gtf and --transcripts-fasta replace --inputs, --sequences and "
     "--background-sequences"),
    (["--gtf", GTF, "--transcripts-fasta", FASTA, "--background-sequences", "x.json"],
     "--gtf and --transcripts-fasta replace --inputs, --sequences and "
     "--background-sequences"),
    (["--gtf", GTF, "--transcripts-fasta", FASTA, "--inputs", "x.json"],
     "--gtf and --transcripts-fasta replace --inputs, --sequences and "
     "--background-sequences"),
    (["--gtf", GTF, "--transcripts-fasta", FASTA, "--ensembl-release", "115"],
     "--ensembl-release 115 contradicts it"),
])
def test_flags_that_do_not_go_together_are_one_line_and_exit_1(tmp_path, capsys, offline,
                                                               argv, says):
    rc, res, err = _cli(tmp_path, capsys, *argv)
    assert rc == 1 and res is None
    assert len(err.strip().splitlines()) == 1 and says in err, err


@pytest.mark.parametrize("kw", [{"sequences": {}}, {"background_sequences": {}}])
def test_the_library_refuses_sequence_beside_the_files(offline, kw):
    with pytest.raises(ValueError, match=I.FILES_REPLACE):
        I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA, **kw)


def test_a_missing_transcript_is_one_line_and_exit_1(tmp_path, capsys, offline):
    t = sorted(_lepr_ids())[0]
    short = _write_fasta(tmp_path / "short.fa", [r for r in RECORDS if not r[0].startswith(t)])
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", short)
    assert rc == 1 and len(err.strip().splitlines()) == 1 and "no record for" in err


def test_the_fasta_alone_serves_a_run_without_the_gene_background(tmp_path, capsys, offline):
    rc, res, err = _cli(tmp_path, capsys, "--transcripts-fasta", FASTA, "--no-gene-background")
    assert rc == 0
    want = I.analyze(_cfg(), transcripts_fasta=FASTA, background_gene_transcripts=False)
    assert _without_provenance(want) == _without_provenance(res)
    # the answer the GTF gives without its background: the same, with no gene named
    with_gtf = I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA,
                         background_gene_transcripts=False)
    assert with_gtf["background"]["gene_id"] is None is res["background"]["gene_id"]
    assert _without_provenance(with_gtf) == _without_provenance(res)
    assert res["annotation"]["source"]["kind"] == "fasta"
    assert res["annotation"]["source"]["gtf"] is None


def test_the_same_file_can_be_the_transcripts_and_the_index(tmp_path, capsys, offline):
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--background-fasta", FASTA, "--keep-duplicates",
                        "--max-window-records", "5")
    assert rc == 0
    bg = res["background"]
    assert bg["gene_id"] == LEPR and bg["fasta"] == FASTA and bg["max_window_records"] == 5
    assert bg["fasta_other_versions"] == {}
    assert "another release" not in err


def test_the_cli_and_the_library_agree(tmp_path, capsys, offline):
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--min-log2fc", "0.5")
    want = I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA, min_log2fc=0.5)
    assert json.dumps(res, sort_keys=True) == json.dumps(want, sort_keys=True, default=str)
    assert rc == cli._identifiability_exit(want)


def test_a_matching_ensembl_release_is_accepted(tmp_path, capsys, offline):
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--ensembl-release", "116")
    assert rc == 0 and res["annotation"]["file_release"] == 116


# --------------------------------------------------------------------------- #
# NOTEs
# --------------------------------------------------------------------------- #
def test_a_config_of_another_release_or_another_gtf_is_said(tmp_path, capsys, offline):
    cfg = _cfg()
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA, cfg=cfg)
    assert "NOTE: the config" not in err
    # the same lines, gzipped again: other bytes
    regz = tmp_path / "gencode.v50.again.gtf.gz"
    regz.write_bytes(gzip.compress(gzip.decompress(pathlib.Path(GTF).read_bytes()), mtime=1))
    rc, res, err = _cli(tmp_path, capsys, "--gtf", str(regz), "--transcripts-fasta", FASTA,
                        cfg=cfg)
    assert rc == 0
    assert ("NOTE: the config was proposed from %s (sha256 %s...), not from this --gtf "
            "(sha256 %s...)" % (cfg["annotation_source"]["file"],
                                cfg["annotation_source"]["sha256"][:12],
                                io.file_sha256(regz)[:12])) in err
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        cfg=dict(cfg, ensembl_release=115))
    assert ("NOTE: the config was annotated against Ensembl release 115, and --gtf %s is "
            "release 116" % GTF) in err


def _bumped_index(tmp_path):
    """The mini FASTA with one configured LEPR transcript at another version: an index of
    another release."""
    cfg = _cfg()
    t = sorted(cfg["groups"][cfg["primary_comparison"][0]])[0]
    vid = BY_ID[t][0]
    new = "%s.%d" % (t, int(vid.split(".")[1]) + 1)
    path = _write_fasta(tmp_path / "index.fa", [(h.replace(vid + "|", new + "|"), s)
                                                for h, s in RECORDS])
    return path, t, vid, new


def test_an_index_of_another_release_is_said_from_the_files(tmp_path, capsys, offline):
    index, t, vid, new = _bumped_index(tmp_path)
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--background-fasta", index)
    assert res["background"]["fasta_other_versions"] == {t: {"used": vid, "fasta": [new]}}
    assert "the index was built from another release" in err and new in err


def test_an_index_of_another_release_is_said_from_rest(tmp_path, capsys, rest):
    index, t, vid, new = _bumped_index(tmp_path)
    rc, res, err = _cli(tmp_path, capsys, "--background-fasta", index)
    assert res["background"]["fasta_other_versions"] == {t: {"used": vid, "fasta": [new]}}
    assert "the index was built from another release" in err


# --------------------------------------------------------------------------- #
# saved inputs
# --------------------------------------------------------------------------- #
def test_saved_inputs_record_the_files_release_and_a_rerun_is_the_file_run(
        tmp_path, capsys, offline, copy_fasta):
    path, a, b = copy_fasta
    saved = tmp_path / "in.json"
    rc, by_files, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", path,
                             "--min-log2fc", "0.5", "--save-inputs", str(saved))
    assert by_files["background"]["identical_source"] == {b: "gene"}
    doc = io.load_inputs(str(saved))
    assert doc["ensembl_release"] == 116
    assert doc["annotation_source"] == by_files["annotation"]["source"]
    assert doc["gene_id"] == LEPR and doc["analysis"]["gene_background"] is True
    assert "saved this run's sequence (Ensembl release 116, from copy.fa)" in err
    rc2, again, err = _cli(tmp_path, capsys, "--inputs", str(saved), "--min-log2fc", "0.5")
    assert rc2 == rc
    assert again["annotation"]["inputs_release"] == 116
    assert _without_provenance(again) == _without_provenance(by_files)


# --------------------------------------------------------------------------- #
# index scope: what a same-name gene in the index is, by the GTF
# --------------------------------------------------------------------------- #
HERC3_OLD, HERC3_NEW = "ENSG00000138641", "ENSG00000287542"


def _herc3(tmp_path, capsys, *argv):
    """HERC3's older gene, against the mini FASTA as the index, which holds the newer gene's
    transcript too: a same-name gene that a GENCODE header does not place."""
    return _cli(tmp_path, capsys, "--background-fasta", FASTA, "--keep-duplicates", *argv,
                cfg=_cfg("HERC3", gene_id=HERC3_OLD))


def test_without_a_gtf_a_same_name_gene_may_be_a_copy(tmp_path, capsys, offline):
    rc, res, err = _herc3(tmp_path, capsys, "--transcripts-fasta", FASTA,
                          "--no-gene-background")
    assert [c["gene_id"] for c in res["background"]["same_name_copies"]] == [HERC3_NEW]
    assert "kind" not in res["background"]["same_name_copies"][0]
    assert "WARNING" in err and "may be copies" in err


def test_with_the_gtf_a_same_name_reference_gene_is_named_and_not_warned_of(
        tmp_path, capsys, offline):
    rc, res, err = _herc3(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA)
    copy, = res["background"]["same_name_copies"]
    assert (copy["gene_id"], copy["region"], copy["kind"]) == (HERC3_NEW, "4",
                                                               "reference_gene")
    assert "WARNING" not in err


def _scaffold_copy(tmp_path):
    """The mini FASTA with a LEPR transcript under a gene id the GTF does not hold, as a
    copy on a scaffold, patch or alternate locus is in GENCODE 48 and later."""
    t = sorted(_lepr_ids())[0]
    head, seq = next((h, s) for h, s in RECORDS if h.startswith(t + "."))
    f = head.split("|")
    f[0], f[1] = "ENST00000999901.1", "ENSG00000999901.1"
    return _write_fasta(tmp_path / "scaffold.fa", RECORDS + [("|".join(f), seq[:-40] + "A" * 40)])


def test_with_the_gtf_a_gene_it_does_not_hold_is_a_copy_or_another_release(
        tmp_path, capsys, offline):
    index = _scaffold_copy(tmp_path)
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--background-fasta", index)
    copy, = res["background"]["same_name_copies"]
    assert (copy["gene_id"], copy["region"], copy["kind"]) == ("ENSG00000999901", None,
                                                               "not_in_gtf")
    assert "WARNING" in err and "may be copies" not in err
    assert "that --gtf %s does not hold" % GTF in err
    assert "or genes of another release" in err


def test_a_gtf_with_the_scaffolds_places_the_copy(tmp_path, capsys, offline):
    """GENCODE's chr_patch_hapl_scaff GTF holds the genes off the reference chromosomes:
    the copy is then placed, as an Ensembl header places it."""
    index = _scaffold_copy(tmp_path)
    extra = ('KI270713.1\tHAVANA\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000999901.1"; '
             'gene_type "protein_coding"; gene_name "LEPR";\n')
    every = tmp_path / "all.gtf"
    every.write_text(gzip.decompress(pathlib.Path(GTF).read_bytes()).decode() + extra)
    rc, res, err = _cli(tmp_path, capsys, "--gtf", str(every), "--transcripts-fasta", FASTA,
                        "--background-fasta", index)
    copy, = res["background"]["same_name_copies"]
    assert (copy["region"], copy["kind"]) == ("KI270713.1", "off_reference")
    assert "on non-reference regions" in err


# --------------------------------------------------------------------------- #
# what the independent reading of this mode found
# --------------------------------------------------------------------------- #
def test_without_the_gene_background_no_gene_is_named_in_either_mode(rest):
    """The gene background is what background.gene_id names: without one, neither mode
    names a gene, so the two reports, and a rerun of either, are one."""
    by_rest = I.analyze(_cfg(), background_gene_transcripts=False)
    by_files = I.analyze(_cfg(), gtf=GTF, transcripts_fasta=FASTA,
                         background_gene_transcripts=False)
    assert by_files["background"]["gene_id"] is None
    assert _without_provenance(by_files) == _without_provenance(by_rest)


@pytest.mark.parametrize("mode", ["files", "rest"])
def test_a_rerun_from_saved_inputs_still_says_the_index_is_of_another_release(
        tmp_path, capsys, monkeypatch, mode):
    monkeypatch.setattr(urllib.request, "urlopen", Rest116())
    index, t, vid, new = _bumped_index(tmp_path)
    saved = tmp_path / "in.json"
    files = ["--gtf", GTF, "--transcripts-fasta", FASTA] if mode == "files" else []
    rc, live, err = _cli(tmp_path, capsys, *files, "--background-fasta", index,
                         "--save-inputs", str(saved))
    assert live["background"]["fasta_other_versions"] == {t: {"used": vid, "fasta": [new]}}
    rc, again, err = _cli(tmp_path, capsys, "--inputs", str(saved), "--background-fasta",
                          index)
    assert again["background"]["fasta_other_versions"] \
        == live["background"]["fasta_other_versions"]
    assert "the index was built from another release" in err


def test_a_saved_annotation_source_that_is_not_one_is_refused(tmp_path, capsys, offline):
    saved = tmp_path / "in.json"
    _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA, "--save-inputs",
         str(saved))
    doc = json.loads(saved.read_text())
    doc["annotation_source"] = {"kind": "gtf"}
    saved.write_text(json.dumps(doc))
    with pytest.raises(io.InputError, match="annotation_source"):
        io.load_inputs(str(saved))
    rc, res, err = _cli(tmp_path, capsys, "--inputs", str(saved), "--save-inputs",
                        str(tmp_path / "again.json"))
    assert rc == 1 and len(err.strip().splitlines()) == 1 and "annotation_source" in err


def test_a_truncated_gzip_index_is_one_line(tmp_path, capsys, offline):
    cut = tmp_path / "index.fa.gz"
    cut.write_bytes(pathlib.Path(FASTA).read_bytes()[:30000])
    rc, res, err = _cli(tmp_path, capsys, "--gtf", GTF, "--transcripts-fasta", FASTA,
                        "--background-fasta", str(cut))
    assert rc == 1 and len(err.strip().splitlines()) == 1, err
    assert "truncated or corrupt" in err
