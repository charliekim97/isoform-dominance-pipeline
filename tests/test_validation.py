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
        (ident, "Over %d estimable genes the homoskedastic factor had a median of %.1f and a "
                "range of %.0fx; the same genes under Poisson weighting had a median of %.1f and "
                "a range of %.0fx. The *ranking* barely moved (Spearman %.2f)"
                % (v["n_estimable"], v["median_conditioning"], v["conditioning_range"],
                   v["poisson_unit_median"], v["poisson_unit_range"],
                   v["poisson_unit_spearman"])),
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
                       "members are %d nt and %d nt long; %s %.3f; %s %.3f)"
                % (v["n_genes"], WORD[len(low)], low[0]["gene"], low[0]["min_jaccard"],
                   *sorted(tpi1), low[1]["gene"], low[1]["min_jaccard"], low[2]["gene"],
                   low[2]["min_jaccard"])))
    return out


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
