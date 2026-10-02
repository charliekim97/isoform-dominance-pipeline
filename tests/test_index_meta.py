"""What `extract` makes of each donor's Salmon ``aux_info/meta_info.json``.

Two indexes built from one transcriptome, one with decoys and one without, share
``index_seq_hash``: Salmon hashes the decoy sequence separately, into
``index_decoy_seq_hash``, which is the SHA-256 of nothing (``e3b0c442...``) when there is
no decoy.  Those two indexes map different reads, so a cohort quantified against both is a
mixed cohort.  An older Salmon writes no decoy fields; then there is nothing to compare.

A meta_info.json that cannot be read is not a third refusal: the documentation promises two
(a mixed cohort, and no configured transcript in any quant.sf).  It is treated as a missing
one, with a warning that names the file.
"""
import json

import pytest

from isoform_dominance import cli, index_scope
from test_index_scope import _meta, _quantdir

NO_DECOY = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def _decoy(seq_hash, decoy_hash, n_decoys):
    return dict(_meta(seq_hash), index_decoy_seq_hash=decoy_hash,
                index_decoy_name_hash="n" + decoy_hash, num_decoy_targets=n_decoys)


# --------------------------------------------------------------------------- #
# A2: one transcriptome, two decoy sets
# --------------------------------------------------------------------------- #
def test_same_transcripts_different_decoys_is_a_mixed_cohort(tmp_path, capsys):
    argv = _quantdir(tmp_path, {"D1": _decoy("58f7", NO_DECOY, 0),
                                "D2": _decoy("58f7", "9a1c", 194)})
    assert cli.main(argv) == 1
    err = capsys.readouterr().err
    assert "different Salmon indexes" in err and "--allow-mixed-index" in err
    assert NO_DECOY in err and "9a1c" in err
    assert len(err.strip().splitlines()) == 1


def test_allow_mixed_index_combines_them_and_records_the_decoys(tmp_path, capsys):
    argv = _quantdir(tmp_path, {"D1": _decoy("58f7", NO_DECOY, 0),
                                "D2": _decoy("58f7", "9a1c", 194)})
    assert cli.main(argv + ["--allow-mixed-index"]) == cli.EXIT_OK
    assert "WARNING" in capsys.readouterr().err
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["mixed"] is True
    assert side["index_seq_hashes"] == ["58f7"]
    assert side["index_decoy_seq_hashes"] == sorted([NO_DECOY, "9a1c"])
    assert side["donors"]["D2"]["index_decoy_seq_hash"] == "9a1c"
    assert side["donors"]["D2"]["num_decoy_targets"] == 194


def test_same_transcripts_same_decoys_is_one_index(tmp_path):
    argv = _quantdir(tmp_path, {"D1": _decoy("58f7", "9a1c", 194),
                                "D2": _decoy("58f7", "9a1c", 194)})
    assert cli.main(argv) == cli.EXIT_OK


def test_a_donor_without_decoy_fields_is_not_compared_on_them(tmp_path):
    # an older Salmon wrote none; that is not evidence of another index
    argv = _quantdir(tmp_path, {"D1": _meta("58f7"), "D2": _decoy("58f7", "9a1c", 194)})
    assert cli.main(argv) == cli.EXIT_OK
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["mixed"] is False
    assert side["donors"]["D1"]["index_decoy_seq_hash"] is None


def test_different_transcripts_is_still_mixed_whatever_the_decoys(tmp_path):
    argv = _quantdir(tmp_path, {"D1": _decoy("58f7", "9a1c", 194),
                                "D2": _decoy("fdc2", "9a1c", 194)})
    assert cli.main(argv) == 1


# --------------------------------------------------------------------------- #
# A3: a meta_info.json that cannot be read is a warning, not a refusal or a traceback
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("content", [
    b"",                                               # empty
    b'{"index_seq_hash": "58f',                        # cut short
    b"[1, 2]",                                         # JSON, not an object
    b'{"index_seq_hash": "\xff\xfe58f7"}',             # not UTF-8
    json.dumps({"index_seq_hash": ["58f7"]}).encode(),  # a hash that is not a string
    json.dumps(dict(_meta("58f7"), index_decoy_seq_hash=7)).encode(),
], ids=["empty", "truncated", "list", "not-utf8", "hash-list", "decoy-hash-int"])
def test_an_unreadable_meta_info_is_one_warning_and_the_table_is_written(
        tmp_path, capsys, content):
    argv = _quantdir(tmp_path, {"D1": None, "D2": _meta("58f7")})
    bad = tmp_path / "quant" / "D1" / "aux_info" / "meta_info.json"
    bad.parent.mkdir()
    bad.write_bytes(content)
    assert cli.main(argv) == cli.EXIT_OK
    err = capsys.readouterr().err
    assert "Traceback" not in err
    lines = [ln for ln in err.strip().splitlines() if str(bad) in ln]
    assert len(lines) == 1 and "WARNING" in lines[0]
    assert (tmp_path / "pd.csv").exists()
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["donors"]["D1"] is None
    assert side["mixed"] is False


def test_the_library_reads_an_unreadable_file_as_none(tmp_path):
    q = tmp_path / "D1" / "quant.sf"
    (tmp_path / "D1" / "aux_info").mkdir(parents=True)
    q.write_text("Name\tTPM\n")
    (tmp_path / "D1" / "aux_info" / "meta_info.json").write_text("[1, 2]")
    prov = index_scope.index_provenance({"D1": str(q)})
    assert prov["donors"]["D1"] is None
    assert list(prov["unreadable_meta_info"]) == ["D1"]
