"""The stored simulation outputs, and every figure the documentation quotes from them (issue #2).

CI runs no quantifier.  What it can hold is the link between the stored outputs and what
ships: (1) the stored files are the ones their READMEs name by SHA-256; (2) the package still
predicts, for each of the 49 genes, what every stored row recorded; (3) ``validation/summary.json``
and the block of ``validation/README.md`` are what ``validation/summarize.py`` computes; and (4)
every figure that the README, the CLI, the package's docstrings, the harness and
``validation/README.md`` quote from ``validation/data`` is the computed one.
"""
import copy
import importlib.util
import json
import pathlib
import re
import shutil

import pytest
from test_saved_identifiability import differences

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


summarize = _module("summarize")


@pytest.fixture(scope="module")
def analysed():
    genes, _, _ = summarize.load()
    return summarize.analyse(genes)


@pytest.fixture(scope="module")
def computed(analysed):
    return summarize.summarize(analysed=analysed)


# --------------------------------------------------------------------------- #
# (1) the stored files
# --------------------------------------------------------------------------- #
#: files example/saved/README.md names by SHA-256 that are not in the repository: the
#: GENCODE files the saved inputs were made from
GENCODE = {"gencode.v44.annotation.gtf.gz", "gencode.v44.transcripts.fa.gz",
           "gencode.v50.annotation.gtf.gz", "gencode.v50.transcripts.fa.gz"}


def stored_files(base):
    """Every file under ``base`` but its README, as paths relative to it: the files the
    README's table must name.  A hidden one (``.DS_Store``, ``._name``) is none."""
    return {p.relative_to(base).as_posix() for p in base.rglob("*")
            if p.is_file() and p.name != "README.md"
            and not any(part.startswith(".") for part in p.relative_to(base).parts)}


@pytest.mark.parametrize("folder", ["example/saved", "validation/data"])
def test_every_file_a_readme_names_is_that_file_and_none_is_unnamed(folder):
    base = ROOT / folder
    table = {k: v for k, v in summarize.sha256_table(base / "README.md").items()
             if k not in GENCODE}
    assert len(table) == 17
    present = stored_files(base)
    assert set(table) == present
    for name, sha in sorted(table.items()):
        assert summarize.file_sha256(base / name) == sha, name


@pytest.mark.parametrize("folder", ["example/saved", "validation/data"])
def test_a_hidden_file_is_not_a_stored_file(folder, tmp_path):
    """A Finder's .DS_Store, or the ._ file it writes beside a file on a foreign volume,
    is not a file the README must name."""
    copy = tmp_path / "copy"
    shutil.copytree(ROOT / folder, copy)
    for d in [copy] + [p for p in copy.iterdir() if p.is_dir()]:
        (d / ".DS_Store").write_bytes(b"\0\0\0\1Bud1")
    (copy / "._README.md").write_bytes(b"\0\5\26\7")
    assert stored_files(copy) == stored_files(ROOT / folder)
    assert len(stored_files(copy)) == 17


def test_a_changed_stored_file_stops_the_summary(tmp_path):
    data = tmp_path / "data"
    shutil.copytree(VALIDATION / "data", data)
    assert summarize.check_files(data) == []
    path = data / "sim" / "c1_d1.json"
    path.write_bytes(path.read_bytes().replace(b'"n": 30', b'"n": 31', 1))
    assert [b[0] for b in summarize.check_files(data)] == ["sim/c1_d1.json"]
    with pytest.raises(SystemExit, match="sim/c1_d1.json"):
        summarize.summarize(data)


# --------------------------------------------------------------------------- #
# (2) the stored predictions are the package's
# --------------------------------------------------------------------------- #
def test_every_stored_prediction_is_what_the_package_says_now(computed):
    assert computed["consistency"] == {"rows_compared": 720, "mismatches": []}
    assert computed["panel"]["n"] == 47
    assert sorted(computed["panel"]["left_out"]) == ["EPOR", "PFKL"]


def test_the_kallisto_panel_is_computed(computed):
    """The panel sentence said kallisto had a row for the same genes from the Salmon count;
    it is the genes with a row in all seven kallisto designs, computed."""
    p = computed["panel"]
    _, rows, _ = summarize.load()
    assert p["kallisto_genes"] == sorted(g for g in rows["kallisto"][0]
                                         if all(g in d for d in rows["kallisto"]))
    assert p["kallisto_n"] == len(p["kallisto_genes"]) == 47
    assert p["kallisto_genes"] == p["genes"]
    assert "(kallisto has a row for the same 47)" in summarize.render(computed)


def test_a_different_kallisto_panel_is_said(computed):
    s = copy.deepcopy(computed)
    s["panel"]["kallisto_genes"] = [g for g in s["panel"]["genes"] if g != "CANX"] + ["EPOR"]
    s["panel"]["kallisto_n"] = len(s["panel"]["kallisto_genes"])
    text = summarize.render(s)
    assert "the same" not in text.split("Left out:")[0]
    assert ("(kallisto has a row in all seven of its designs for 47 genes: not CANX, and "
            "also EPOR)") in text


def test_a_prediction_that_moved_is_named(analysed):
    _, rows, _ = summarize.load()
    structure = copy.deepcopy({g: a["structural"] for g, a in analysed.items()})
    assert summarize.consistency(structure, rows)[0] == []
    structure["LEPR"]["min_resolvable_log2fc"] *= 1 + 1e-7
    structure["CANX"]["verdict"] = "weakly_identifiable"
    moved = {(m[0], m[2]) for m in summarize.consistency(structure, rows)[0]}
    assert moved == {("LEPR", "min_resolvable_log2fc"), ("LEPR", "predicted_se"),
                     ("CANX", "verdict")}


# --------------------------------------------------------------------------- #
# (3) the committed summary is the computed one
# --------------------------------------------------------------------------- #
def test_summary_json_is_what_summarize_computes(computed):
    want = json.loads((VALIDATION / "summary.json").read_text(encoding="utf-8"))
    assert differences(json.loads(summarize.as_json(computed)), want, "summary") == []


def test_the_readme_block_is_the_rendered_summary(computed):
    text = (VALIDATION / "README.md").read_text(encoding="utf-8")
    block = text.split(summarize.START)[1].split(summarize.END)[0]
    assert block == "\n" + summarize.render(computed) + "\n"


def test_a_changed_stored_figure_changes_the_summary(tmp_path, analysed):
    """One stored bias changed, with its file's SHA-256 updated to match, and the summary
    no longer is the committed one."""
    data = tmp_path / "data"
    shutil.copytree(VALIDATION / "data", data)
    path = data / "sim" / "c1_d1.json"
    doc = json.loads(path.read_text())
    row = next(r for r in doc["genes"] if r["gene"] == "LEPR")
    row["bias"] = -row["bias"]
    path.write_text(json.dumps(doc, indent=1))
    readme = data / "README.md"
    old = summarize.sha256_table(VALIDATION / "data" / "README.md")["sim/c1_d1.json"]
    readme.write_text(readme.read_text().replace(old, summarize.file_sha256(path)))
    want = json.loads((VALIDATION / "summary.json").read_text(encoding="utf-8"))
    got = json.loads(summarize.as_json(summarize.summarize(data, analysed=analysed)))
    assert any(d.startswith("summary.S2.") for d in differences(got, want, "summary"))


# --------------------------------------------------------------------------- #
# (4) every quoted figure is the computed one
# --------------------------------------------------------------------------- #
def text_of(path):
    """A file's text as one line: line breaks, comment markers and the seams of implicitly
    concatenated string literals taken out, runs of white space made one space."""
    t = path.read_text(encoding="utf-8")
    if path.suffix == ".py":
        t = re.sub(r'"\s*\n\s*"', "", t)
        t = re.sub(r"\n\s*#:?", " ", t)
    return " ".join(t.split())


def _span(values, places=0, dash="–"):
    lo, hi = "%.*f" % (places, min(values)), "%.*f" % (places, max(values))
    return lo if lo == hi else lo + dash + hi


WORD = {2: "two", 3: "three"}


def quotes(s):
    """``[(file, the text it must contain)]``, each text built from the computed summary."""
    s1, s2, s7, s3, s8, v = s["S1"], s["S2"], s["S7"], s["S3"], s["S8"], s["survey"]
    five, three = summarize.FIVE, summarize.THREE
    mild = [summarize.B.index(-0.3), summarize.B.index(0.3)]
    u = summarize.UNIFORM
    rows = {(r["old_gate"], r["structure"], r["conditioning"]): r for r in s["table"]}
    good, high = rows[("pass", "estimable", "<= 10")], rows[("pass", "estimable", "> 10")]
    fail = rows[("fail", "estimable", "> 10")]
    blind = rows[("pass", "not estimable", "inf")]
    note = s["table_notes"]
    s4 = s["S4"]["summaries"]
    b2, bm2 = summarize.B.index(2.0), summarize.B.index(-2.0)

    def direction(d, i):
        return [d["five_prime_rule"][j] for j in i], [d["three_prime_rule"][j] for j in i]

    s2_five, s2_three = direction(s2, five)[0], direction(s2, three)[1]
    s7_five, s7_three = direction(s7, five)[0], direction(s7, three)[1]

    def salmon_kallisto(dash):
        return ("Salmon's error took the 5' direction for %s of the %d genes with a ratio "
                "above 1.23x and the 3' direction for %s of them, where at uniform coverage the "
                "5' direction held for %d of %d and the 3' for %d. kallisto, on the same reads, "
                "took the 5' direction for %s and the 3' for %s, but it has no bias model and "
                "leaned to the longer class at uniform coverage too (%d of %d)"
                % (_span(s2_five, 0, dash), s2["n"], _span(s2_three, 0, dash),
                   s2["five_prime_rule"][u], s2["n"], s2["three_prime_rule"][u],
                   _span(s7_five, 0, dash), _span(s7_three, 0, dash),
                   s7["three_prime_rule"][u], s7["n"]))

    readme = "README.md"
    cli = "src/isoform_dominance/cli.py"
    ident = "src/isoform_dominance/identifiability.py"
    sim = "validation/simulate_reads.py"
    vreadme = "validation/README.md"
    out = [
        # README, the skew direction (S2, S7) and LEPR (S3)
        (readme, "in the simulation its bias ran against that direction under %d of the %d "
                 "skews" % (s3["against_rule"], s3["skewed_designs"])),
        (readme, "in a %d-gene simulation (monotone positional skew) Salmon's error took the 5' "
                 "direction for %s of the %d genes above that ratio and the 3' direction for %s, "
                 "where at uniform coverage the 5' direction held for %d of %d and the 3' for %d"
                 % (v["n_genes"], _span(s2_five), s2["n"], _span(s2_three),
                    s2["five_prime_rule"][u], s2["n"], s2["three_prime_rule"][u])),
        (readme, "kallisto, on the same reads, took the 5' direction for %s and the 3' for %s, "
                 "but it has no bias model and leaned to the longer class at uniform coverage "
                 "too, %d of %d (the 5' direction %d)"
                 % (_span(s7_five), _span(s7_three), s7["three_prime_rule"][u], s7["n"],
                    s7["five_prime_rule"][u])),
        # README, spread and not accuracy (S1) and the table
        (readme, "In simulation (Salmon, %d genes, 30 replicates each; %d genes estimable with a "
                 "finite predicted SE in all seven coverage models" % (v["n_genes"], s1["n"])),
        (readme, "no gene's |bias| reaches its own predicted SE at uniform coverage: %d of %d"
                 % (s1["abs_bias_over_predicted_se"][u], s1["n"])),
        (readme, "%s of %d exceed it, depending on the direction of the skew"
                 % (_span([s1["abs_bias_over_predicted_se"][i] for i in mild]), s1["n"])),
        (readme, "stays below the predicted SE in %s of %d, as it did in %d of %d at uniform "
                 "coverage" % (_span([s1["realised_sd_under_predicted_se"][i] for i in mild]),
                               s1["n"], s1["realised_sd_under_predicted_se"][u], s1["n"])),
        (readme, "NTRK3 holds a replicate SD of %s while its bias runs from %.3f to %.2f log2"
                 % (_span(s1["NTRK3"]["realised_sd"], 3), *s1["NTRK3"]["abs_bias"])),
        (readme, "on the %d genes with a row in every design, the conditioning factor did not "
                 "sort genes by error: the %d genes the gate passes with conditioning above 10 "
                 "had a median RMSE at uniform coverage of %.3f against %.3f for the %d at or "
                 "below it, and a smaller worst bias under skew, %.2f against %.2f log2"
                 % (s["panel"]["n"], high["n"], high["uniform_rmse_median"],
                    good["uniform_rmse_median"], good["n"], high["worst_abs_bias_median"],
                    good["worst_abs_bias_median"])),
        (readme, "The gate failed %s genes the structure calls estimable (%s) and passed all %d "
                 "it calls not estimable, whose uniform RMSE (median %.3f) was no worse than the "
                 "estimable genes'"
                 % (WORD[fail["n"]], ", ".join(fail["genes"]), blind["n"],
                    blind["uniform_rmse_median"])),
        (readme, "as a group they are deeper by construction, their median gene %.1f times the "
                 "estimable genes' in fragments" % s8["median_gene_fragments"]["ratio"]),
        # the NOTE identify prints, and the --min-log2fc help (the survey)
        (cli, "in a %d-gene simulation with monotone positional skew, %s"
              % (v["n_genes"], salmon_kallisto("-"))),
        (cli, "over a %d-gene survey the median gene sat at conditioning %.0f against the "
              "default 10" % (v["n_genes"], v["median_conditioning"])),
        # identifiability.py
        (ident, "In the %d-gene simulation the direction was scored on the genes whose class "
                "effective lengths differ by more than 1.23x (log2 0.2987); 0.3 selects the "
                "same %d genes there" % (v["n_genes"], s2["n"])),
        (ident, "In a %d-gene simulation (monotone positional coverage skew) Salmon split the "
                "ambiguous mass between the classes by effective length" % v["n_genes"]),
        (ident, "for %s and %s of the %d genes whose ratio exceeds 1.23x, against %d and %d of "
                "%d at uniform coverage. kallisto, on the same reads, gave %s and %s, but it has "
                "no bias model and leaned to the longer class at uniform coverage too (%d of %d)"
                % (_span(s2_five, 0, "-"), _span(s2_three, 0, "-"), s2["n"],
                   s2["five_prime_rule"][u], s2["three_prime_rule"][u], s2["n"],
                   _span(s7_five, 0, "-"), _span(s7_three, 0, "-"),
                   s7["three_prime_rule"][u], s7["n"])),
        (ident, "ran against it under %d of %d skews with Salmon"
                % (s3["against_rule"], s3["skewed_designs"])),
        (ident, "reads %+.3f at b = 0" % s4["class total"]["spearman"][u]),
        # gls_covariance: the two medians are not set side by side, as theta sets the scale
        (ident, "Over %d estimable genes the homoskedastic factor had a median of %.1f and a "
                "range of %.0fx." % (v["n_estimable"], v["median_conditioning"],
                                     v["conditioning_range"])),
        (ident, "over the same genes a median of %.1f with theta = 1 (one expected fragment) "
                "per transcript" % v["poisson_unit_median"]),
        (ident, "at any such scale the range is %.0fx against %.0fx and the ranking barely "
                "moves (Spearman %.2f with the homoskedastic factor)"
                % (v["poisson_unit_range"], v["conditioning_range"],
                   v["poisson_unit_spearman"])),
        (ident, "and of %.2f with the design's expected fragment counts as theta"
                % v["poisson_design_median"]),
        (ident, "With the design's expected counts, which differ from one transcript to "
                "another, they are %.1fx and %.3f."
                % (v["poisson_design_range"], v["poisson_design_spearman"])),
        (ident, "over a %d-gene survey the median gene sat at %.1f against a default tau of "
                "10.0, so the default rejects %.0f%% of what it"
                % (v["n_genes"], v["median_conditioning"], 100 * v["fraction_over_tau"])),
        (ident, "a %d-gene survey put the median gene at a conditioning factor of %.1f against "
                "the default 10.0" % (v["n_genes"], v["median_conditioning"])),
        (ident, "on TPI1 it reads %.2f where the class totals themselves are not resolvable"
                % v["TPI1_counting_min_resolvable_log2fc"]),
        (ident, "at uniform coverage in that simulation the replicate SD was below it for %d of "
                "the %d estimable genes of its %d-gene panel"
                % (s1["realised_sd_under_predicted_se"][u], s1["n"], s["panel"]["n"])),
        # the harness
        (sim, "`not_identifiable` group is %.1fx deeper than the estimable one by construction "
              "(smaller class %.0f vs %.0f fragments, %.1fx; counting floor differs at "
              "Mann-Whitney p = %.4f)"
              % (s8["median_gene_fragments"]["ratio"],
                 s8["median_smaller_class_fragments"]["not_identifiable"],
                 s8["median_smaller_class_fragments"]["estimable"],
                 s8["median_smaller_class_fragments"]["ratio"],
                 s8["counting_floor_mannwhitney_p"])),
        (sim, "None of the %d transcripts in the current survey is short enough"
              % v["n_transcripts"]),
        # validation/README.md, outside the generated block
        (vreadme, "the genes with a row in all seven Salmon designs, %d of the %d"
                  % (s["panel"]["n"], v["n_genes"])),
        # said only while it is so: a kallisto panel of other genes has no such sentence
        (vreadme, "kallisto has a row for the same %d in all seven" % s["panel"]["kallisto_n"]
                  if s["panel"]["kallisto_genes"] == s["panel"]["genes"] else
                  "kallisto has a row in all seven of its designs for %d genes"
                  % s["panel"]["kallisto_n"]),
        (vreadme, "**Row 1** holds %d genes: %s." % (good["n"], ", ".join(good["genes"]))),
        (vreadme, "The %d genes with conditioning above 10 had no clearly larger error at "
                  "uniform coverage than the %d at or below it (median RMSE %.3f against %.3f, "
                  "excess %.2f against %.2f; Mann-Whitney p %.2f and %.2f), and under skew "
                  "their worst |bias| was smaller (%.2f against %.2f log2, p %.2f)"
                  % (high["n"], good["n"], high["uniform_rmse_median"],
                     good["uniform_rmse_median"], high["uniform_excess_median"],
                     good["uniform_excess_median"], s["table_comparison"]["uniform_rmse"],
                     s["table_comparison"]["uniform_excess"],
                     high["worst_abs_bias_median"], good["worst_abs_bias_median"],
                     s["table_comparison"]["worst_abs_bias"])),
        (vreadme, "**Row 3** holds %d genes, %s, whose median uniform-coverage RMSE, %.3f, is "
                  "the table's largest (their median excess, %.2f, is below the not-estimable "
                  "row's %.2f). ESR1's worst |bias|, %.2f log2, is at bias_5p -2.0, where only "
                  "%d of its 30 replicates had a defined ratio"
                  % (fail["n"], " and ".join(fail["genes"]), fail["uniform_rmse_median"],
                     fail["uniform_excess_median"], blind["uniform_excess_median"],
                     note["ESR1"]["worst_abs_bias"],
                     note["ESR1"]["replicates_with_a_ratio"][bm2])),
        (vreadme, "the %d genes the gate passes and the structure calls not estimable had no "
                  "worse RMSE at uniform coverage than the estimable genes (median %.3f against "
                  "%.3f; NTRK3 %.3f). Their large excess (%.2f)"
                  % (blind["n"], blind["uniform_rmse_median"],
                     note["estimable_uniform_rmse_median"], note["NTRK3_uniform_rmse"],
                     blind["uniform_excess_median"])),
        (vreadme, "as a group they are deeper by construction (a median gene %.1f times the "
                  "estimable genes'), so most of them have a smaller counting floor"
                  % s8["median_gene_fragments"]["ratio"]),
        (vreadme, "S8 gives a %.1f-fold larger median gene and an %.1f-fold larger smaller class"
                  % (s8["median_gene_fragments"]["ratio"],
                     s8["median_smaller_class_fragments"]["ratio"])),
    ]
    # the class_efflen_ratio table: Spearman at b = +2.0 and -2.0, and the 5' sign
    for name, x in s4.items():
        out.append((ident, "%s %+.3f %+.3f %d/%d" % (
            name, x["spearman"][b2], x["spearman"][bm2], x["five_prime_sign_right"],
            x["n_over_1_23_fold"])))
    # class_coherence: the classes under 0.05, and TPI1's two members
    low = v["classes_under_0_05_min_jaccard"]
    tpi1 = s["genes"]["TPI1"]["structural"]["transcript_lengths"][low[0]["class"]]
    out.append((ident, "Across a %d-gene survey %s classes came in under 0.05 (%s %.3f, whose two "
                       "members are %d nt and %d nt long; %s %.3f; %s %.3f), while the median "
                       "of the same figure, ``min_jaccard``, over the %d classes of its "
                       "primary comparisons that have two or more members was %.2f"
                % (v["n_genes"], WORD[len(low)], low[0]["gene"], low[0]["min_jaccard"],
                   *sorted(tpi1), low[1]["gene"], low[1]["min_jaccard"], low[2]["gene"],
                   low[2]["min_jaccard"], v["n_classes_with_pairs"], v["median_min_jaccard"])))
    # the two estimands: where they part on the panel, which the README and the gate quote
    only = v["estimable_as_abundance_only"]
    out.append((readme, "On the %d-gene panel of `validation/` the contrast is estimable as "
                        "abundance only for %d genes (%s), and as fragment totals only for %s."
                % (v["n_genes"], len(only), ", ".join(only),
                   ", ".join(v["estimable_as_fragment_only"]) or "none")))
    return out


def test_the_gate_grades_against_predictions_fixed_with_it(computed):
    """The 3.0 gate grades against gate_3_0_predictions.json, written with it on 2026-10-23:
    the bytes its SHA-256 names, the 7 genes of the primary set and the 39 of S1, and the
    predictions this code gives, which the README rounds for the primary set.  A later change
    of the code that moves a prediction fails here; the gate still grades against the file."""
    import hashlib
    import math
    path = VALIDATION / "gate_3_0_predictions.json"
    text = text_of(VALIDATION / "README.md")
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    assert "`gate_3_0_predictions.json` beside this file" in text and sha in text
    frozen = json.loads(path.read_text(encoding="utf-8"))
    assert list(frozen["primary"]) == computed["survey"]["estimable_as_abundance_only"]
    assert list(frozen["control"]) == computed["S1"]["genes"]
    assert frozen["design"] == computed["designs"]["predicted_at"]
    for part in ("primary", "control"):
        for gene, se in frozen[part].items():
            now = computed["genes"][gene]["structural"]["abundance"]["min_resolvable_log2fc"]
            assert math.isfinite(se) and se == pytest.approx(now / 1.96, rel=1e-9), gene
    assert ("For the primary set they are %s, rounded here"
            % ", ".join("%s %.3f" % (g, se) for g, se in frozen["primary"].items())
            ).replace(", VEGFA", " and VEGFA") in text
    # the sets graded are the file's, named in the text as they were on the day
    assert ("The 7 genes of the file's `primary` (%s: on 2026-10-23 the genes whose contrast is "
            "estimable as abundance and not as fragment totals)"
            % ", ".join(frozen["primary"])) in text
    assert "The 39 genes of the file's `control` (S1's on 2026-10-23)" in text


def test_the_gate_says_its_predictions_are_of_a_nearby_functional(monkeypatch):
    """The predictions are of the abundance with window positions as lengths.  With
    L - 199, near the effective length summed TPM divides by, in their place, none of the
    primary set's contrasts is estimable in the window system: the gate quotes the range of
    the residuals and the gene with the largest."""
    import numpy as np
    ident = summarize.I
    frozen = json.loads((VALIDATION / "gate_3_0_predictions.json").read_text(encoding="utf-8"))
    genes, _, _ = summarize.load()
    window = ident.DEFAULT_K
    original = ident._coefficients

    def fragment_lengths(name, c, n_positions):
        if name != "abundance":
            return original(name, c, n_positions)
        length = np.asarray(n_positions, dtype=float) + window - 1
        return original(name, c, np.maximum(1.0, length - 199.0))
    monkeypatch.setattr(ident, "_coefficients", fragment_lengths)
    residual = {}
    for gene in frozen["primary"]:
        cfg, conf, rest = summarize._split(genes[gene])
        res = ident.analyze(cfg, sequences=conf, background_sequences=rest,
                            background_gene_transcripts=False, **summarize.DESIGN)
        c = res["estimands"]["abundance"]["contrast"]
        assert not c["estimable"], gene
        residual[gene] = c["residual"]
    low, high = min(residual.values()), max(residual.values())
    assert max(residual, key=residual.get) == "VEGFA"
    text = text_of(VALIDATION / "README.md")
    assert ("none of the 7 contrasts is estimable in the window system (residual %.0e to %.0e, "
            "VEGFA's the largest)" % (low, high)).replace("e-0", "e-") in text


def test_the_gate_quotes_the_stored_errors_and_the_abundance_predictions(computed):
    """The 3.0 gate of validation/README.md quotes the stored Salmon runs and its fixed
    predictions: SMN2 and VEGFA of its stop rule, both of the primary set; the three
    predictions past the linearisation limit and how S1 graded such predictions; and the
    thresholds of the primary set and of the control, against S1."""
    import math
    _, rows, _ = summarize.load()
    uniform = rows["salmon"][summarize.UNIFORM]
    frozen = json.loads((VALIDATION / "gate_3_0_predictions.json").read_text(encoding="utf-8"))
    only = list(frozen["primary"])
    assert {"SMN2", "VEGFA"} <= set(only)
    text = text_of(VALIDATION / "README.md")

    def se(gene):
        return frozen["primary"].get(gene, frozen["control"].get(gene))
    smn2, vegfa = uniform["SMN2"], uniform["VEGFA"]
    # what makes them the stop rule: each already failed with NumReads
    assert abs(smn2["bias"]) > se("SMN2") and smn2["realised_sd"] > se("SMN2")
    assert abs(vegfa["bias"]) > se("VEGFA") and vegfa["realised_sd"] < se("VEGFA")
    assert ("SMN2 (|bias| %.3f, SD %.3f, against a predicted SE of %.3f)"
            % (abs(smn2["bias"]), smn2["realised_sd"], se("SMN2"))) in text
    assert "VEGFA (|bias| %.3f against %.3f)" % (abs(vegfa["bias"]), se("VEGFA")) in text
    # the other five of the primary set met the standard with NumReads
    for gene in set(only) - {"SMN2", "VEGFA"}:
        r = uniform[gene]
        assert abs(r["bias"]) < se(gene) and r["realised_sd"] < se(gene), gene
    assert ("5 of the 7 met the abundance standard with `NumReads` although their fragment "
            "contrast is not estimable") in text

    # past the linearisation limit: the relative SE on the natural-log scale, by which the
    # CLI stars a figure; and S1's predictions past it, which S1 graded all the same
    ident = summarize.I
    rse = {gene: se(gene) * math.log(2) for gene in only}
    for gene in only:
        assert ident.min_resolvable_log2fc(rse[gene], 1) == pytest.approx(1.96 * se(gene))
    past = sorted((g for g in only if rse[g] > ident.LINEARISATION_LIMIT), key=rse.get)
    assert ("%s, %s and %s, at a relative SE of %.3f, %.3f and %.3f on the natural-log scale "
            "against the limit of %.1f" % (*past, *(rse[g] for g in past),
                                           ident.LINEARISATION_LIMIT)) in text
    s1 = computed["S1"]
    s1_past = sum(computed["genes"][g]["structural"]["beyond_linear"] for g in s1["genes"])
    assert ("as S1 graded its %d genes, %d of whose predictions are past the limit too"
            % (s1["n"], s1_past)) in text
    assert "TNFRSF1B's predicted SE, %.3f log2" % se("TNFRSF1B") in text

    # the thresholds: of the primary set against S1's standard, and of the control
    k = summarize.UNIFORM
    assert s1["abs_bias_over_predicted_se"][k] == 0
    assert ("every one of the 7 has a |bias| below its predicted SE, at least 6 of the 7 a "
            "replicate SD below it") in text
    assert ("where no |bias| reached the predicted SE and the SD stayed below it for %d"
            % s1["realised_sd_under_predicted_se"][k]) in text
    assert ("a |bias| below the predicted SE for at least 38, and an SD below it for at least "
            "%d, the fewest any design of S1 gave with `NumReads`"
            % min(s1["realised_sd_under_predicted_se"])) in text
    met = (sum(abs(uniform[g]["bias"]) < se(g) for g in s1["genes"]),
           sum(uniform[g]["realised_sd"] < se(g) for g in s1["genes"]))
    assert ("Against the abundance predictions, `NumReads` at uniform coverage met them for "
            "%d and %d" % met) in text
    assert met[0] >= 38 and met[1] >= min(s1["realised_sd_under_predicted_se"])

def test_the_survey_gives_the_median_coherence_of_the_classes_with_pairs(computed, analysed):
    """class_coherence's docstring quoted "the median class sat near 0.7", which no
    computation gave; the median min_jaccard of the primary-comparison classes with two or
    more members is computed instead, from the survey's own coherence figures."""
    v = computed["survey"]
    js = sorted(j for a in analysed.values() for j in a["coherence"].values() if j is not None)
    assert v["n_classes_with_pairs"] == len(js) == 74
    assert v["median_min_jaccard"] == (js[36] + js[37]) / 2
    assert round(v["median_min_jaccard"], 4) == 0.4925
    assert sum(len(a["coherence"]) for a in analysed.values()) == 98


def test_the_poisson_figures_are_at_theta_one_and_at_the_design(computed):
    """gls_covariance's docstring gives the Poisson-weighted factor at theta = 1 per
    transcript (summarize.poisson_unit_factor): 16.40, 156.6x and 0.930; and with the
    design's expected counts as theta (summarize.poisson_design_factor): 353.74, 116.8x and
    0.922."""
    v = computed["survey"]
    assert (round(v["poisson_unit_median"], 2), round(v["poisson_unit_range"], 1),
            round(v["poisson_unit_spearman"], 3)) == (16.40, 156.6, 0.930)
    assert (round(v["poisson_design_median"], 2), round(v["poisson_design_range"], 1),
            round(v["poisson_design_spearman"], 3)) == (353.74, 116.8, 0.922)


def test_the_poisson_factor_takes_its_scale_from_theta(analysed):
    """What gls_covariance's docstring says of theta's scale, computed here apart from
    summarize: theta = s for every transcript gives sqrt(s) times the figure at theta = 1,
    so the same range and ranking; the design's expected fragment counts, which differ by
    transcript, give the figures summarize.poisson_design_factor gives."""
    import math
    import statistics

    import numpy as np
    ident = summarize.I

    def factor(A, c, theta):
        cov = ident.gls_covariance(A, np.asarray(theta, dtype=float))
        return math.sqrt(float(c @ cov @ c))
    genes, _, _ = summarize.load()
    have = sorted(g for g in genes if analysed[g]["poisson_unit_factor"] is not None)
    cond = [analysed[g]["structural"]["conditioning_factor"] for g in have]
    d = summarize.DESIGN
    design = {}
    for g in have:
        cfg, conf, rest = summarize._split(genes[g])
        groups = {k: [t.split(".")[0] for t in v] for k, v in cfg["groups"].items()}
        pc = cfg["primary_comparison"]
        ids = [t for v in groups.values() for t in v] + sorted(rest)
        seqs = dict(conf, **rest)
        A, _ = ident.compatibility_matrix(
            {t: ident.kmer_track(seqs[t], ident.DEFAULT_K) for t in ids}, ids)
        c = np.array([1.0 if t in groups[pc[0]] else -1.0 if t in groups[pc[1]] else 0.0
                      for t in ids])
        assert factor(A, c, [4.0] * len(ids)) == pytest.approx(
            2.0 * analysed[g]["poisson_unit_factor"], rel=1e-9)
        design[g] = factor(A, c, [ident.expected_informative_reads(
            1.0, d["tpm"], len(seqs[t]), depth=d["depth"], mean_efflen=d["mean_efflen"])
            for t in ids])
        assert design[g] == pytest.approx(analysed[g]["poisson_design_factor"], rel=1e-12)
    got = [design[g] for g in have]
    median, spread = statistics.median(got), max(got) / min(got)
    rho = summarize.H.spearman(cond, got)
    assert (round(median, 2), round(spread, 1), round(rho, 3)) == (353.74, 116.8, 0.922)


def test_every_quoted_figure_is_the_computed_one(computed):
    s = computed
    # facts the quotes rest on without printing them
    assert s["S2"]["same_genes_as_over_1_23_fold"] is True
    assert not s["S1"]["NTRK3"]["in_S1"]
    assert {r["n"] for r in s["table"] if r["old_gate"] == "fail"
            and r["structure"] == "not estimable"} == {0}
    assert s["survey"]["shortest_transcript"] >= 100        # simulate_reads.READ_LEN
    texts, missing = {}, []
    for path, quote in quotes(s):
        if path not in texts:
            texts[path] = text_of(ROOT / path)
        if quote not in texts[path]:
            missing.append("%s: %s" % (path, quote))
    assert missing == []


def test_a_changed_figure_is_a_missing_quote(computed):
    """The quotes are built from the summary: a summary one figure off finds a quote that
    no file holds."""
    s = copy.deepcopy(computed)
    texts = {p: text_of(ROOT / p) for p, _ in quotes(s)}
    assert all(q in texts[p] for p, q in quotes(s))
    for change in (lambda x: x["S2"]["five_prime_rule"].__setitem__(5, 37),
                   lambda x: x["S7"]["three_prime_rule"].__setitem__(3, 31),
                   lambda x: x["S8"]["median_smaller_class_fragments"].__setitem__(
                       "estimable", 343),
                   lambda x: x["survey"].__setitem__("median_conditioning", 66.0)):
        t = copy.deepcopy(computed)
        change(t)
        assert not all(q in texts[p] for p, q in quotes(t)), change
