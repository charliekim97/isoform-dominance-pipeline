"""extract: transcripts the config names that are not in quant.sf.

A configured transcript is absent from quant.sf when Salmon kept another copy of an
identical sequence as it built the index -- the chrY copy of a pseudoautosomal gene, in a
GENCODE-built index, is removed as a duplicate of the chrX one -- or when the config and
the index come from different releases.  Summing what is there gave a class total short
by the missing transcripts, or zero for every class, and nothing said so: a CD99 config
proposed from the chrY gene came out as 0.0000 TPM and a ratio of NA for every donor.
"""
import json
import os

import pytest

from isoform_dominance import cli, extract, index_scope

GROUPS = {"A": ["ENST00000000001", "ENST00000000002"], "B": ["ENST00000000003"]}


def _case(tmp_path, donors, groups=GROUPS):
    """``donors``: {donor: [the transcript names in its quant.sf]}."""
    qd = tmp_path / "quant"
    for d, names in donors.items():
        (qd / d).mkdir(parents=True)
        (qd / d / "quant.sf").write_text(
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
            + "".join("%s\t1000\t800\t%d\t10\n" % (n, 5 + i) for i, n in enumerate(names)))
        # one index for every donor, recorded as Salmon records it, so that the only
        # warning a test can see is the one it is about
        (qd / d / "aux_info").mkdir()
        (qd / d / "aux_info" / "meta_info.json").write_text(json.dumps(
            {"salmon_version": "1.10.3", "index_seq_hash": "58f7", "index_name_hash": "6765",
             "num_valid_targets": 642692, "keep_duplicates": False}))
    sm = tmp_path / "sm.csv"
    sm.write_text("donor,condition\n" + "".join("%s,control\n" % d for d in donors))
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "GENEX", "groups": groups,
                               "primary_comparison": ["A", "B"]}))
    return ["extract", "--config", str(cfg), "--quantdir", str(qd), "--samplemap", str(sm),
            "--cohort", "C", "--out", str(tmp_path / "pd.csv")]


def _sidecar(tmp_path):
    return json.load(open(str(tmp_path / "pd.csv") + ".index.json"))


ALL = ["ENST00000000001.4", "ENST00000000002.1", "ENST00000000003.2", "ENST00000000099.1"]
# what a GENCODE-built index holds for a config proposed from a pseudoautosomal gene's
# chrY copy: the chrX transcripts, under other ids, and none of the configured ones
NONE = ["ENST00000000011.1", "ENST00000000012.1", "ENST00000000013.1"]


def test_no_configured_transcript_in_any_donor_stops(tmp_path, capsys):
    rc = cli.main(_case(tmp_path, {"D1": NONE, "D2": NONE}))
    err = capsys.readouterr().err
    assert rc == 1
    assert "none of the 3 transcripts the config names" in err
    assert "duplicate_clusters.tsv" in err          # the usual cause, named
    assert "release" in err                         # ... and the other one
    assert not os.path.exists(tmp_path / "pd.csv")
    assert not os.path.exists(str(tmp_path / "pd.csv") + ".index.json")


def test_the_library_call_refuses_too(tmp_path):
    _case(tmp_path, {"D1": NONE})
    cfg = json.load(open(tmp_path / "cfg.json"))
    with pytest.raises(index_scope.NoConfiguredTranscripts):
        extract.run(cfg, str(tmp_path / "quant"), str(tmp_path / "sm.csv"), "C",
                    str(tmp_path / "pd.csv"))


def test_some_missing_warns_once_and_is_recorded(tmp_path, capsys):
    names = ["ENST00000000001.4", "ENST00000000003.2"]            # ...0002 is missing
    rc = cli.main(_case(tmp_path, {"D1": names, "D2": names, "D3": names}))
    err = capsys.readouterr().err
    assert rc == cli.EXIT_OK
    assert err.count("WARNING") == 1                              # once, not per donor
    assert "1 of the 3 transcripts the config names" in err
    assert "ENST00000000002" in err
    assert "duplicate_clusters.tsv" in err
    side = _sidecar(tmp_path)["missing_transcripts"]
    assert side["n_configured"] == 3
    assert side["n_missing"] == 1
    assert side["ids"] == ["ENST00000000002"]
    assert "ENST00000000002" in side["warning"]
    assert os.path.exists(tmp_path / "pd.csv")


def test_a_long_list_is_cut_in_the_warning_but_not_in_the_sidecar(tmp_path, capsys):
    groups = {"A": ["ENST%011d" % i for i in range(1, 11)], "B": ["ENST00000000100"]}
    names = ["ENST00000000001.1", "ENST00000000100.1"]           # 9 of 11 missing
    assert cli.main(_case(tmp_path, {"D1": names}, groups=groups)) == cli.EXIT_OK
    err = capsys.readouterr().err
    assert "9 of the 11 transcripts the config names" in err
    assert "and 4 more" in err                                    # the first five shown
    assert "ENST00000000010" not in err
    assert len(_sidecar(tmp_path)["missing_transcripts"]["ids"]) == 9


def test_missing_from_some_donors_only_says_how_many(tmp_path, capsys):
    rc = cli.main(_case(tmp_path, {"D1": ALL, "D2": ["ENST00000000001.4",
                                                      "ENST00000000003.2"]}))
    err = capsys.readouterr().err
    assert rc == cli.EXIT_OK
    assert "1 of the 3 transcripts the config names" in err
    assert "in 1 of 2 donors" in err
    assert _sidecar(tmp_path)["missing_transcripts"]["n_missing"] == 1


def test_one_donor_with_none_is_a_warning_when_another_has_them(tmp_path, capsys):
    """Only a cohort in which no donor has any configured transcript is refused; one
    donor without them is what a mixed-index cohort looks like, and is reported."""
    assert cli.main(_case(tmp_path, {"D1": ALL, "D2": NONE})) == cli.EXIT_OK
    err = capsys.readouterr().err
    assert "3 of the 3 transcripts the config names" in err
    assert "in 1 of 2 donors" in err


def test_all_present_is_silent(tmp_path, capsys):
    assert cli.main(_case(tmp_path, {"D1": ALL, "D2": ALL})) == cli.EXIT_OK
    assert "WARNING" not in capsys.readouterr().err
    assert _sidecar(tmp_path)["missing_transcripts"] == {
        "n_configured": 3, "n_missing": 0, "ids": [], "warning": None}


def test_full_gencode_header_names_count_as_present(tmp_path, capsys):
    names = ["ENST00000000001.4|ENSG00000000010.1|-|-|GENEX-201|GENEX|1000|protein_coding|",
             "ENST00000000002.1|ENSG00000000010.1|-|-|GENEX-202|GENEX|1000|protein_coding|",
             "ENST00000000003.2|ENSG00000000010.1|-|-|GENEX-203|GENEX|1000|protein_coding|"]
    assert cli.main(_case(tmp_path, {"D1": names})) == cli.EXIT_OK
    assert "WARNING" not in capsys.readouterr().err
    assert _sidecar(tmp_path)["missing_transcripts"]["n_missing"] == 0


# --------------------------------------------------------------------------- #
# a transcript in two groups: extract warns and sums as it always has
# --------------------------------------------------------------------------- #
def test_a_transcript_in_two_groups_is_named_and_summed_as_before(tmp_path, capsys):
    # the 2.1.1 aggregation puts its TPM in the group listed last, and that is kept --
    # changing it would change published numbers; the warning is what is new
    shared = GROUPS["A"][0]
    groups = {"A": GROUPS["A"], "B": GROUPS["B"] + [shared]}
    argv = _case(tmp_path, {"D1": ALL[:3]}, groups=groups)     # TPM 5, 6, 7
    assert cli.main(argv) == 0
    err = capsys.readouterr().err
    line = [ln for ln in err.splitlines() if shared in ln]
    assert len(line) == 1 and "WARNING" in line[0]
    assert '"A"' in line[0] and '"B"' in line[0]
    head, row = open(tmp_path / "pd.csv").read().splitlines()[:2]
    tpm = dict(zip(head.split(","), row.split(","), strict=True))
    assert (tpm["A_TPM"], tpm["B_TPM"]) == ("6.0000", "12.0000")
