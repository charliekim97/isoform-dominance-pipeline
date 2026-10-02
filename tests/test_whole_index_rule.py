"""Which FASTA records join the system, and why the answer is never more optimistic.

A record of ``--background-fasta`` that holds a configured transcript's window found in at
most ``max_window_records`` (M) records is a column of the compatibility system; every
other record that holds a window of the system is left out, and every window it holds is
dropped from every layer.  On GENCODE 50's reference chromosomes, one column per record
that shares a window would have meant a median of 2,058 columns per config and 67,121 at
most; a few windows, in repeats, are found in very many records.

Write the exact whole-index system as every record a column.  The windows kept are the
ones no left-out record holds, so the system here is the exact one with the rows that
touch a left-out record removed.  Removing rows only shrinks the row space: what is
estimable here is estimable there.  And profiling out the left-out records adds a positive
semi-definite Schur term to the information of the rest: the standard error here is never
smaller.  These tests compare the two, with the exact system built by passing every record
as ``background_sequences``.
"""
import random

import pytest

from isoform_dominance import cli
from isoform_dominance import identifiability as I

K = 15


def _rand(r, n):
    return "".join(r.choice("ACGT") for _ in range(n))


def _mutate(r, s, n):
    s = list(s)
    for _ in range(n):
        i = r.randrange(len(s))
        s[i] = r.choice("ACGT".replace(s[i], ""))
    return "".join(s)


def _system(r, second_order=True, gene_only=True):
    """A gene, its background, and FASTA records: (config, sequences, gene, records)."""
    E = [_rand(r, r.randint(60, 140)) for _ in range(6)]
    rep = _rand(r, 40)                                   # a repeat many records hold
    T1, T2, T3 = "ENST00000000001", "ENST00000000002", "ENST00000000003"
    seqs = {T1: E[0] + E[1] + E[2] + (rep if r.random() < 0.8 else ""),
            T2: E[0] + E[3] + E[2]}
    groups = {"A": [T1], "B": [T2]}
    if r.random() < 0.4:
        seqs[T3] = E[0] + E[1] + E[4]
        groups["A"].append(T3)
    gene = {"ENST00000000011": E[3] + E[5], "ENST00000000012": E[2] + _rand(r, 80)}
    records, own = {}, []
    for i in range(r.randint(3, 22)):
        kind = r.choices(["paralogue", "repeat", "gene", "second", "copy", "random"],
                         [3, 5, 2 if gene_only else 0, 3 if second_order else 0, 1, 1])[0]
        mine = _rand(r, r.randint(30, 90))
        if kind == "paralogue":
            piece = r.choice([E[1], E[3], E[4], E[0][:50] + E[1][:50]])
            seq = _mutate(r, piece, r.randint(0, 3)) + mine
        elif kind == "repeat":
            seq = mine[:20] + rep + mine[20:]
        elif kind == "gene":
            seq = E[5][:r.randint(40, len(E[5]))] + mine
        elif kind == "second" and own:
            seq = r.choice(own) + mine                   # shares an earlier record's own part
        elif kind == "copy" and records:
            seq = r.choice(list(records.values()))
        else:
            seq = mine + _rand(r, 40)
        own.append(mine)
        records["ENST9%010d" % i] = seq
    cfg = {"gene": "GENEX", "groups": groups, "primary_comparison": ["A", "B"]}
    return cfg, seqs, gene, records


def _fasta(tmp_path, records, n):
    p = tmp_path / ("bg%d.fa" % n)
    p.write_text("".join(">%s.1\n%s\n" % (t, s) for t, s in records.items()))
    return str(p)


def _connected(seqs, gene, records):
    """The records that share a window with a column, or with a record that does, and so
    on: in the exact system the others are blocks of their own, which add to the rank and
    bear on no configured estimand."""
    windows = set().union(*(I.kmers(x, K) for x in list(seqs.values()) + list(gene.values())))
    held = {t: I.kmers(x, K) for t, x in records.items()}
    joined = set()
    while True:
        new = {t for t, w in held.items() if t not in joined and w & windows}
        if not new:
            return {t: x for t, x in records.items() if t in joined}
        joined |= new
        for t in new:
            windows |= held[t]


def _estimands(res):
    return [res["groups"]["A"], res["groups"]["B"], res["contrast"]]


def _run(cfg, seqs, gene, **kw):
    kw.setdefault("background_gene_transcripts", False)
    return I.analyze(cfg, k=K, window=K, sequences=seqs, background_sequences=gene,
                     min_log2fc=0.5, **kw)


def _compare(here, exact):
    """(a) and (b); True when here is strictly more conservative somewhere."""
    stricter = False
    for h, e in zip(_estimands(here), _estimands(exact), strict=True):
        if h["estimable"]:
            assert e["estimable"], "estimable here but not in the exact system"
        if h["estimable"] and e["estimable"]:
            for key in ("gls_relative_se", "min_resolvable_log2fc"):
                assert h[key] >= e[key] * (1 - 1e-9), (key, h[key], e[key])
                stricter |= h[key] > e[key] * (1 + 1e-9)
        stricter |= e["estimable"] and not h["estimable"]
    return stricter


def test_never_more_optimistic_than_every_record_a_column(tmp_path):
    r = random.Random(20261002)
    runs = stricter = dropped = equal_when_none_left = 0
    for n in range(70):
        cfg, seqs, gene, records = _system(r)
        fa = _fasta(tmp_path, records, n)
        exact = _run(cfg, seqs, dict(gene, **records))
        connected = _run(cfg, seqs, dict(gene, **_connected(seqs, gene, records)))
        for e, c in zip(_estimands(exact), _estimands(connected), strict=True):
            assert e["estimable"] == c["estimable"]
            assert e["gls_relative_se"] == pytest.approx(c["gls_relative_se"], rel=1e-9)
        for M in sorted({0, 1, 2, 3, 5, len(records), len(records) + 1}):
            here = _run(cfg, seqs, gene, background_fasta=fa, max_window_records=M)
            runs += 1
            stricter += _compare(here, exact)
            dropped += here["background"]["windows_dropped"]["total"] > 0
            if here["background"]["fasta_left_out"] == 0:
                # (c) nothing left out: the exact system itself, less its decoupled blocks
                for h, e in zip(_estimands(here), _estimands(connected), strict=True):
                    assert h["rank"] == e["rank"]
                    assert h["gls_relative_se"] == pytest.approx(e["gls_relative_se"],
                                                                 rel=1e-9)
                assert cli._identifiability_exit(here) == cli._identifiability_exit(exact)
                equal_when_none_left += 1
    # not vacuous: windows were dropped, the answer was stricter, and the exact case ran
    assert dropped > runs // 3 and stricter > 30 and equal_when_none_left >= 10, (
        runs, dropped, stricter, equal_when_none_left)


def test_with_m_above_every_count_and_nothing_else_left_out_it_is_the_exact_system(
        tmp_path):
    """(c): with no record holding gene-background windows alone and none sharing only
    another record's windows, M at or above the largest count leaves nothing out."""
    r = random.Random(7)
    for n in range(40):
        cfg, seqs, gene, records = _system(r, second_order=False, gene_only=False)
        fa = _fasta(tmp_path, records, n)
        exact = _run(cfg, seqs, dict(gene, **_connected(seqs, gene, records)))
        here = _run(cfg, seqs, gene, background_fasta=fa, max_window_records=len(records))
        assert here["background"]["fasta_left_out"] == 0
        for h, e in zip(_estimands(here), _estimands(exact), strict=True):
            assert (h["rank"], h["estimable"]) == (e["rank"], e["estimable"])
            assert h["gls_relative_se"] == pytest.approx(e["gls_relative_se"], rel=1e-9)
        assert cli._identifiability_exit(here) == cli._identifiability_exit(exact)


def test_a_competitor_found_once_is_always_added(tmp_path):
    r = random.Random(3)
    cfg, seqs, gene, _ = _system(r)
    one = {"ENST90000000001": list(seqs.values())[0][20:120] + _rand(r, 60)}
    fa = _fasta(tmp_path, one, 0)
    for M in (1, 2, 20):
        here = _run(cfg, seqs, gene, background_fasta=fa, max_window_records=M)
        assert list(here["background"]["fasta_competitors"]) == list(one)
        assert here["background"]["fasta_left_out"] == 0
    here = _run(cfg, seqs, gene, background_fasta=fa, max_window_records=0)
    assert here["background"]["fasta_competitors"] == {}
    assert here["background"]["fasta_left_out"] == 1


def test_a_record_left_out_takes_the_windows_it_holds_from_every_layer(tmp_path):
    """The repeat in T1 is in more than M records: its windows go from the unique count,
    the read model and the system alike, and the NOTE says how many."""
    r = random.Random(11)
    rep = _rand(r, 40)
    seqs = {"ENST00000000001": _rand(r, 199) + "A" + rep + "C" + _rand(r, 99),
            "ENST00000000002": _rand(r, 300)}
    cfg = {"gene": "GENEX", "groups": {"A": ["ENST00000000001"], "B": ["ENST00000000002"]},
           "primary_comparison": ["A", "B"]}
    # flanked by other bases than in T1, so that the repeat's windows are all they share
    records = {"ENST9%010d" % i: _rand(r, 49) + "G" + rep + "T" + _rand(r, 49)
               for i in range(4)}
    fa = _fasta(tmp_path, records, 0)
    here = _run(cfg, seqs, {}, background_fasta=fa, max_window_records=3)
    bg = here["background"]
    held = I.kmers(rep, K) & I.kmers(seqs["ENST00000000001"], K)
    assert bg["windows_dropped"]["configured"] == len(held)
    assert bg["windows_dropped"]["columns"] == {"ENST00000000001": len(held)}
    assert bg["fasta_left_out"] == 4 and bg["fasta_competitors"] == {}
    alone = _run(cfg, seqs, {})
    assert here["groups"]["A"]["n_unique_kmers"] == alone["groups"]["A"]["n_unique_kmers"] \
        - len(held)
    assert here["groups"]["A"]["expected_informative_reads"] \
        < alone["groups"]["A"]["expected_informative_reads"]
    with_four = _run(cfg, seqs, {}, background_fasta=fa, max_window_records=4)
    assert with_four["background"]["windows_dropped"]["total"] == 0
    assert list(with_four["background"]["fasta_competitors"]) == sorted(records)


def test_max_window_records_must_be_a_count():
    with pytest.raises(ValueError, match="max_window_records"):
        I.analyze({"groups": {"A": ["T1"], "B": ["T2"]}}, sequences={"T1": "A" * 50,
                                                                      "T2": "C" * 50},
                  background_gene_transcripts=False, max_window_records=-1)
    with pytest.raises(SystemExit):
        cli.main(["identifiability", "--config", "x.json", "--max-window-records", "-1"])


def test_a_column_whose_windows_are_all_dropped_is_no_precondition_failure(
        tmp_path, capsys):
    """A gene-background transcript that is all repeat loses every window to records left
    out.  It is not shorter than the window, so it is no reason for exit 2: the gene total
    leaves it out and names it, and the configured estimands are as they would be without
    it."""
    import json
    r = random.Random(5)
    rep = _rand(r, 60)
    seqs = {"ENST00000000001": _rand(r, 300), "ENST00000000002": _rand(r, 300)}
    cfg = {"gene": "GENEX", "groups": {"A": ["ENST00000000001"], "B": ["ENST00000000002"]},
           "primary_comparison": ["A", "B"]}
    gene = {"ENST00000000011": rep}
    records = {"ENST9%010d" % i: _rand(r, 40) + rep + _rand(r, 40) for i in range(3)}
    fa = _fasta(tmp_path, records, 0)
    here = _run(cfg, seqs, gene, background_fasta=fa, max_window_records=2)
    assert here["gene_total"]["transcripts_all_windows_dropped"] == ["ENST00000000011"]
    assert here["gene_total"]["transcripts_without_windows"] == []
    assert here["gene_total"]["estimable"] is True
    assert cli._identifiability_exit(here) != cli.EXIT_NOT_IDENTIFIABLE
    alone = _run(cfg, seqs, {})
    for h, a in zip(_estimands(here), _estimands(alone), strict=True):
        assert h["gls_relative_se"] == pytest.approx(a["gls_relative_se"], rel=1e-9)
    cp, sq, bs = tmp_path / "c.json", tmp_path / "s.json", tmp_path / "b.json"
    cp.write_text(json.dumps(cfg))
    sq.write_text(json.dumps(seqs))
    bs.write_text(json.dumps(gene))
    cli.main(["identifiability", "--config", str(cp), "--sequences", str(sq),
              "--background-sequences", str(bs), "--background-fasta", fa, "--window", "15",
              "--k", "15", "--max-window-records", "2"])
    assert "every window of ENST00000000011 was dropped" in capsys.readouterr().err
