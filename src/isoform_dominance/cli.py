"""Unified command-line interface: isoform-dominance <subcommand>.

Every subcommand accepts ``--json``, which writes the full result object to stdout
instead of the human summary, so the tool composes inside a workflow manager without
anyone having to parse its printed text.

``identifiability`` does not put its structural verdict in the exit status. The verdict
moves with the annotation release the transcripts come from, so it is reported as a
field and a line of output, and the status answers only the question the user asked
with ``--min-log2fc``:

===== ==============================================================
  0   no ``--min-log2fc`` given, or both class totals and the
      contrast resolve it at the stated design
  3   ``--min-log2fc`` given and not resolved -- including when an
      estimand has no finite figure at all, or one past the
      linearisation limit (``beyond_linear``)
  2   precondition failure: the gene total itself is not estimable,
      because a transcript shorter than ``--window`` has no windows
      and an all-zero column. Not a verdict; checked first
  1   config error or network failure
===== ==============================================================
"""
import argparse
import http.client
import json
import math
import os
import sys

from . import (__version__, annotate, contamination, ensembl, extract, identifiability,
               index_scope, io, stats)

EXIT_OK = 0
EXIT_NOT_IDENTIFIABLE = 2
EXIT_EFFECT_NOT_RESOLVED = 3


def _identifiability_exit(res):
    """Exit status from the precondition and the requested effect size, never the verdict."""
    if not res["gene_total"]["estimable"]:
        return EXIT_NOT_IDENTIFIABLE
    if res["effect_resolvable"] is False:
        return EXIT_EFFECT_NOT_RESOLVED
    return EXIT_OK


def _release(value):
    """argparse type for --ensembl-release: a positive integer."""
    try:
        n = int(value)
    except ValueError:
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError("expected a positive release number, got %r" % value)
    return n


def _kv(items):
    out = {}
    for s in (items or []):
        if "=" not in s:
            raise SystemExit("argument error: expected NAME=path, got %r" % s)
        k, v = s.split("=", 1)
        if not k or not v:
            raise SystemExit("argument error: expected NAME=path, got %r" % s)
        out[k] = v
    return out


def _emit(payload):
    json.dump(payload, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")


def _net_fail(e):
    print("Ensembl request failed (%s). Check your network connection, the gene symbol, "
          "and species, or supply sequences offline." % e, file=sys.stderr)
    return 1


def _release_fail(e):
    print("Ensembl release not available: %s" % e, file=sys.stderr)
    return 1


def _bad_url(e):
    print("error: the Ensembl request could not be sent (%s); check the gene symbol, "
          "species and ids for characters a URL cannot carry" % e, file=sys.stderr)
    return 1


def cmd_annotate(a):
    try:
        cfg = annotate.run(a.gene, a.out, species=a.species, release=a.ensembl_release,
                           gene_id=a.gene_id, retries=a.retries, retry_wait=a.retry_wait)
    except ensembl.ReleaseNotServed as e:
        return _release_fail(e)
    except http.client.InvalidURL as e:
        return _bad_url(e)
    except ensembl.TRANSIENT as e:
        return _net_fail(e)
    except ValueError as e:
        print("error: %s" % e, file=sys.stderr)
        return 1
    # the notes go to stderr with --json too: they are what the JSON alone does not say
    choice = cfg.get("_gene_choice") or {}
    if choice.get("reason") and choice["reason"] != "given by --gene-id":
        print("  NOTE: %s" % choice["reason"], file=sys.stderr)
    proposal = cfg.get("_proposal", {})
    ties = proposal.get("tied_with") or []
    if ties:
        alt = cfg["primary_comparison"][0]
        # quote the recorded rule rather than paraphrase one step of it: a tie on
        # transcript count is settled by protein length and then by the acceptor
        # coordinate, and a note that stops at the protein does not say what decided this
        print("  NOTE: %s was chosen by a tie. %d other cluster(s) have as many transcripts "
              "(%d): %s. The rule is %s, so a release that adds one transcript to either "
              "can change the proposed comparison. Pick the one you mean; each is listed "
              "under _clusters."
              % (alt, len(ties), ties[0]["n_transcripts"],
                 ", ".join("%d aa at acceptor %d" % (t["rep_protein_aa"],
                                                    t["terminal_acceptor"]) for t in ties),
                 proposal.get("alternative_rule") or annotate.ALTERNATIVE_RULE),
              file=sys.stderr)
    if a.json:
        _emit(cfg)
        return EXIT_OK
    print("Proposed groups for %s -> %s" % (a.gene, a.out))
    for g, ids in cfg["groups"].items():
        print("  %s: %d transcripts" % (g, len(ids)))
    print("  primary_comparison:", cfg["primary_comparison"])
    print("  REVIEW _proposed/_clusters and rename groups before use.")
    return EXIT_OK


def _load_sequences(path, flag):
    """``{transcript_id: sequence}`` from the JSON file a sequence flag names."""
    seqs = io.load_json(path, flag + " file")
    if not (isinstance(seqs, dict)
            and all(isinstance(k, str) and isinstance(v, str) for k, v in seqs.items())):
        raise io.InputError("%s file %s is not a JSON object of transcript id to sequence"
                            % (flag, path))
    return seqs


def _fasta_note(saved, given, flag="--background-fasta",
                without="windows shared with other genes count as unique here"):
    """What to say when a rerun from saved inputs has another --background-fasta, or none.

    ``saved`` is the ``background_fasta`` record of the saved inputs (None if that run
    used none); ``given`` is this run's --background-fasta.  None when they agree.  With
    ``flag`` and ``without``, the same for another file the run read, such as --decoys.
    """
    if saved is None:
        return None
    if not given:
        return ("the saved run also used %s %s (sha256 %s...); without it, %s and the "
                "verdict can differ. Pass it again."
                % (flag, saved["path"], saved["sha256"][:12], without))
    if io.file_sha256(given) != saved["sha256"]:
        return ("%s %s is not the file the inputs were saved with (%s, sha256 %s...); the "
                "verdict can differ." % (flag, given, saved["path"], saved["sha256"][:12]))
    return None


def _analysis_note(saved, k, window, canonical, keep_duplicates=None):
    """What to say when a rerun from saved inputs is at another window or convention,
    or -- ``keep_duplicates`` is None when this run has no background, gene or FASTA --
    counts background sequences identical to a configured transcript or to one another
    differently.

    None of these is recorded in the config, so nothing else can tell a rerun that it is
    not repeating the run that was saved.  ``k`` only sets the window when ``window`` is
    None: every layer is built from the window, so a k that differs with the same window
    is the same system.
    """
    mine = {"window": window if window is not None else k, "canonical": canonical}
    notes = []
    if {key: saved[key] for key in mine} != mine:
        def _say(d):
            return ("window=%d, %s k-mers"
                    % (d["window"], "canonical" if d["canonical"] else "strand-aware"))
        notes.append("the inputs were saved at %s; this run is at %s, which is a different "
                     "compatibility system on the same sequence, so the verdict can differ."
                     % (_say(saved), _say(mine)))
    if (keep_duplicates is not None and "keep_duplicates" in saved
            and saved["keep_duplicates"] != keep_duplicates):
        notes.append("the inputs were saved %s --keep-duplicates and this run is %s it: a "
                     "background sequence identical to a configured transcript or to another "
                     "was %s then and is %s now, so the verdict can differ."
                     % (("with", "without", "counted", "not counted") if saved["keep_duplicates"]
                        else ("without", "with", "not counted", "counted")))
    return " ".join(notes) or None


def _load_saved_inputs(a, cfg):
    """The saved inputs ``--inputs`` names, split by *this* config; exits on misuse.

    A saved file holds the sequence of one grouping: the configured transcripts under
    ``sequences`` and the rest of the gene under ``background_sequences``.  Which is
    which is a property of the config, not of the file, so the two are pooled and split
    again here by the config in hand -- the rule a live run follows.  Splitting them by
    the saved grouping instead dropped a transcript the config no longer names from the
    background as well as from its class, and refused one moved the other way for "no
    sequence" while holding that sequence.
    """
    if a.sequences or a.background_sequences or a.ensembl_release is not None:
        raise ValueError("--inputs replaces --sequences, --background-sequences and "
                         "--ensembl-release; pass it alone")
    inputs = io.load_inputs(a.inputs)
    if inputs.get("gene") and cfg.get("gene") and inputs["gene"] != cfg["gene"]:
        raise ValueError("the saved inputs are for %s and the config for %s"
                         % (inputs["gene"], cfg["gene"]))
    # one symbol can name two genes: the chrY copy of a pseudoautosomal gene has its own id
    saved_id, cfg_id = inputs.get("gene_id"), cfg.get("gene_id")
    if saved_id and cfg_id and saved_id.split(".")[0] != cfg_id.split(".")[0]:
        raise ValueError("the saved inputs are for gene %s and the config for gene %s"
                         % (saved_id, cfg_id))
    # ids compared without their version, as --sequences ids are
    pool = {t.split(".")[0]: s for t, s in inputs["background_sequences"].items()}
    pool.update((t.split(".")[0], s) for t, s in inputs["sequences"].items())
    needed = {t.split(".")[0] for ids in cfg["groups"].values() for t in ids}
    absent = sorted(needed - set(pool))
    if absent:
        # never fall through to Ensembl here: the point of --inputs is a run that does not
        # depend on it, and a silent fetch would mix releases
        raise ValueError("the saved inputs have no sequence for %s, which the config names; "
                         "they were saved for another grouping" % ", ".join(absent))
    rest = {t: s for t, s in pool.items() if t not in needed}
    # a file written before 2.4.1 does not say; a background it saved is the evidence
    had = inputs["analysis"].get("gene_background", bool(inputs["background_sequences"]))
    return dict(inputs,
                sequences={t: pool[t] for t in sorted(needed)},
                background_sequences=rest, gene_background=had,
                left_out=[] if had else sorted(rest))


def cmd_identifiability(a):
    cfg = io.load_config(a.config)
    if a.save_inputs:
        # checked before the run, which can spend minutes on an archive, not after it
        if os.path.isdir(a.save_inputs):
            raise io.InputError("--save-inputs %s: is a directory; name the file to write"
                                % a.save_inputs)
        where = os.path.dirname(os.path.abspath(a.save_inputs))
        if not os.path.isdir(where):
            raise io.InputError("--save-inputs %s: directory %s does not exist"
                                % (a.save_inputs, where))
        if not os.access(where, os.W_OK):
            raise io.InputError("--save-inputs %s: directory %s is not writable"
                                % (a.save_inputs, where))
    inputs = None
    try:
        if a.inputs:
            inputs = _load_saved_inputs(a, cfg)
    except ValueError as e:
        print("config error: %s" % e, file=sys.stderr)
        return 1
    if inputs is not None:
        seqs = inputs["sequences"]
        # no gene background for a rerun of a run that had none, as for a live one
        background = (None if a.no_gene_background or not inputs["gene_background"]
                      else inputs["background_sequences"])
    else:
        seqs = _load_sequences(a.sequences, "--sequences") if a.sequences else None
        background = (_load_sequences(a.background_sequences, "--background-sequences")
                      if a.background_sequences else None)
    captured = {} if a.save_inputs else None
    try:
        res = identifiability.analyze(
            cfg, k=a.k, sequences=seqs,
            canonical=not a.strand_aware,
            window=a.window,
            background_sequences=background,
            background_fasta=a.background_fasta, decoys=a.decoys,
            background_gene_transcripts=(False if a.no_gene_background or inputs is not None
                                         else "auto"),
            read_length=a.read_length, frag_mean=a.frag_mean, frag_sd=a.frag_sd,
            paired=not a.single_end, depth=a.depth, tpm=a.tpm, n_donors=a.donors,
            conditioning_tau=a.tau, min_informative_reads=a.min_informative_reads,
            min_log2fc=a.min_log2fc, ensembl_release=a.ensembl_release,
            inputs_out=captured, keep_duplicates=a.keep_duplicates,
            retries=a.retries, retry_wait=a.retry_wait)
    except ensembl.ReleaseNotServed as e:
        return _release_fail(e)
    except http.client.InvalidURL as e:
        return _bad_url(e)
    except ensembl.TRANSIENT as e:
        return _net_fail(e)
    except ValueError as e:
        print("config error: %s" % e, file=sys.stderr)
        return 1
    except MemoryError as e:
        print("error: out of memory%s.%s" % (
            " (%s)" % e if str(e) else "",
            " Each --background-fasta record that shares a window with a configured "
            "transcript is a column of the compatibility system, so genome sequence (pass "
            "--decoys) or a repeat that many records share can make the system too large."
            if a.background_fasta else ""), file=sys.stderr)
        return 1
    if inputs is not None:
        res["annotation"]["inputs_release"] = inputs["ensembl_release"]
        # a saved sequence reaches the run as supplied; the file says which were fetched
        was = inputs.get("sequence_sources") or {}
        for t, src in (res["background"].get("identical_source") or {}).items():
            if src == "sequences" and was.get(t, "").startswith("fetched:"):
                res["background"]["identical_source"][t] = "gene"
        if background is not None and res["background"]["gene_id"] is None:
            # the gene the live run fetched its background as, which the file recorded
            res["background"]["gene_id"] = inputs.get("gene_id")
        if inputs["left_out"]:
            print("  NOTE: the saved run had no gene background, so %s, saved but not named "
                  "by this config, %s left out rather than used as background, as a live "
                  "run with --no-gene-background leaves %s out."
                  % (", ".join(inputs["left_out"]),
                     "is" if len(inputs["left_out"]) == 1 else "are",
                     "it" if len(inputs["left_out"]) == 1 else "them"), file=sys.stderr)
        for note in (_analysis_note(inputs["analysis"], a.k, a.window, not a.strand_aware,
                                    a.keep_duplicates if a.background_fasta
                                    or background is not None else None),
                     _fasta_note(inputs.get("background_fasta"), a.background_fasta),
                     _fasta_note(inputs.get("decoys"), a.decoys, "--decoys",
                                 "the genome records it names are read as competing "
                                 "sequence")):
            if note:
                print("  NOTE: " + note, file=sys.stderr)
    code = _report(a, res)
    # after the report, so a save that fails cannot take the answer with it
    if captured is not None:
        sources = captured["sequence_sources"]
        if inputs is not None:
            # a rerun supplies its saved sequence: where it came from is the file's record
            was = inputs.get("sequence_sources") or {}
            default = ("fetched:%s" % inputs["ensembl_release"]
                       if inputs["ensembl_release"] is not None else "supplied")
            sources = captured["sequence_sources"] = {t: was.get(t, default) for t in sources}
        kinds = sorted(set(sources.values()))
        # one release only when every sequence came from it
        release = (int(kinds[0].split(":", 1)[1])
                   if len(kinds) == 1 and kinds[0].startswith("fetched:") else None)
        if captured.get("gene_id") is None and inputs is not None:
            captured["gene_id"] = inputs.get("gene_id")
        io.save_inputs(a.save_inputs, captured, cfg, release, __version__,
                       background_fasta=a.background_fasta, decoys=a.decoys)
        n_sup = sum(v == "supplied" for v in sources.values())
        print("  saved this run's sequence (%s) to %s; repeat it with no network: "
              "--inputs %s" % ("Ensembl release %s" % release if release is not None
                               else "release unknown: supplied offline" if n_sup == len(sources)
                               else "%d supplied, %d fetched from %s; no one release recorded"
                               % (n_sup, len(sources) - n_sup, ", ".join(
                                   "release " + k.split(":", 1)[1] for k in kinds
                                   if k != "supplied")),
                               a.save_inputs, a.save_inputs), file=sys.stderr)
    return code


def _listed(ids, n=5):
    """The first ``n`` of ``ids``, and how many more."""
    ids = sorted(ids)
    return ", ".join(ids[:n]) + (" and %d more" % (len(ids) - n) if len(ids) > n else "")


def _report(a, res):
    """Print the report of one run -- JSON on stdout with ``--json`` -- and return the exit
    status."""
    bg = res["background"]
    same = bg.get("identical_to_configured") or {}
    twins = bg.get("identical_to_background") or {}
    where = bg.get("identical_source") or {}
    # a copy with the sequence of a configured transcript, or of a background one, takes no
    # read the configured transcripts would get: the NOTEs below, not this warning
    warning = index_scope.copy_warning(
        index_scope.without_identical(bg["same_name_copies"], dict(same, **twins)),
        "--background-fasta %s" % a.background_fasta)
    if warning:
        print("  " + warning, file=sys.stderr)

    def _pairs(d):
        ids = sorted(d)
        return ", ".join("%s (= %s)" % (r, d[r]) for r in ids[:5]) + (
            " and %d more" % (len(ids) - 5) if len(ids) > 5 else "")
    from_fasta = {r: t for r, t in same.items() if where.get(r, "fasta") == "fasta"}
    from_gene = {r: t for r, t in same.items() if r not in from_fasta}
    if from_fasta:
        print("  NOTE: --background-fasta: %d record(s) with the sequence of a configured "
              "transcript were not counted as competing sequence, because Salmon's default "
              "index keeps one of identical sequences: %s. For an index built with Salmon's "
              "--keepDuplicates, pass --keep-duplicates."
              % (len(from_fasta), _pairs(from_fasta)), file=sys.stderr)
    if from_gene:
        print("  NOTE: the gene background: %d transcript(s) with the sequence of a "
              "configured transcript were not counted as competing sequence, because "
              "Salmon's default index keeps one of identical sequences: %s. For an index "
              "built with Salmon's --keepDuplicates, pass --keep-duplicates."
              % (len(from_gene), _pairs(from_gene)), file=sys.stderr)
    if same:
        print("  NOTE: of identical sequences, Salmon's default index keeps the one that comes "
              "first in the FASTA it is built from. Where that is the copy above rather than "
              "the configured transcript, quant.sf names only the copy and `extract`, which "
              "sums the configured ids, does not count those reads. Put the copy in the same "
              "group as the transcript it equals, or build the index with Salmon's "
              "--keepDuplicates and pass --keep-duplicates here.", file=sys.stderr)
    if twins:
        print("  NOTE: %d background sequence(s) identical to another were counted once, as "
              "the first of them, because Salmon's default index keeps one of identical "
              "sequences: %s. For an index built with Salmon's --keepDuplicates, pass "
              "--keep-duplicates." % (len(twins), _pairs(twins)), file=sys.stderr)
    long = bg.get("fasta_long_records") or {}
    if long and not a.decoys:
        print("  NOTE: --background-fasta has %d record(s) longer than 1 Mb (%s), which looks "
              "like genome sequence: the gentrome of a decoy-aware index. Salmon sets aside "
              "only the reads that map better to a decoy than to any transcript, so a genome "
              "record competes with no read the transcripts explain as well, and judged "
              "against it every window inside an exon loses its uniqueness. Pass --decoys "
              "decoys.txt, or the transcript FASTA the index was built from (without the "
              "genome decoys)." % (len(long), _listed(long)), file=sys.stderr)
    elif long:
        print("  NOTE: --background-fasta has %d record(s) longer than 1 Mb that --decoys %s "
              "does not name (%s); they were read as competing sequence. If they are genome "
              "sequence, add them to it." % (len(long), a.decoys, _listed(long)),
              file=sys.stderr)
    if bg.get("decoys_absent"):
        print("  NOTE: --decoys %s names %d record(s) that --background-fasta does not hold "
              "(%s); %d decoy record(s) were left out." % (
                  a.decoys, len(bg["decoys_absent"]), _listed(bg["decoys_absent"]),
                  bg["decoys_skipped"]), file=sys.stderr)
    replaced = res["background"].get("sequence_from_fasta") or []
    if replaced:
        print("  NOTE: --background-fasta holds other sequence for %d transcript(s) of the "
              "gene background, as another release would: %s. Each is counted once, with its "
              "sequence taken from --background-fasta, which is what an index built from that "
              "FASTA holds." % (len(replaced), _listed(replaced)), file=sys.stderr)

    if a.json:
        _emit(res)
        return _identifiability_exit(res)

    # every layer is built from the window; k only sets it when --window is not given
    print("Identifiability (window=%d, %s k-mers)"
          % (res["window"], "canonical" if res["canonical"] else "strand-aware"))
    if a.window is not None and a.k != a.window:
        print("  NOTE: --k %d is not used: --window %d sets the length of every window, "
              "k-mer and unique stretch here." % (a.k, a.window), file=sys.stderr)
    rel = res["annotation"]["ensembl_release"]
    got = res["annotation"]["fetched_release"]
    saved = res["annotation"].get("inputs_release")
    print("  annotation: %s%s%s"
          % ("Ensembl release %s" % rel if rel is not None
             else "release not recorded in the config",
             "; sequence fetched from release %s" % got
             if got is not None and got != rel else "",
             "; sequence from saved inputs (%s, %s)"
             % ("release %s" % saved if saved is not None else "release not recorded",
                a.inputs) if a.inputs else ""))
    if a.inputs and rel is not None and saved != rel:
        print("  NOTE: the config was annotated against Ensembl release %s, and the saved "
              "inputs %s." % (rel, "are from release %s" % saved if saved is not None
                              else "record no one release: some or all of their sequence "
                              "was supplied, not fetched"),
              file=sys.stderr)
    if rel is None and not (a.inputs and saved is not None):
        # a rerun from inputs that record their release repeats that release's answer
        print("  NOTE: the config records no Ensembl release, and the verdict is a function "
              "of the release its transcripts came from, so this verdict is not "
              "reproducible. Re-run `annotate` to record one, or set \"ensembl_release\" in "
              "the config.", file=sys.stderr)
    elif got is not None and got != rel:
        if a.ensembl_release is None:
            print("  NOTE: the config was annotated against Ensembl release %s, but the "
                  "server now serves release %s and the sequence used here came from it. "
                  "The groups may name transcripts whose sequence has changed, and the gene "
                  "may have gained transcripts they do not name. To fetch from release %s, "
                  "pass --ensembl-release %s." % (rel, got, rel, rel), file=sys.stderr)
        else:
            print("  NOTE: sequence was fetched from Ensembl release %s, as requested, but "
                  "the config's groups were proposed from release %s. The gene may have "
                  "transcripts in release %s that the groups do not name; `annotate "
                  "--ensembl-release %s` proposes groups from it." % (got, rel, got, got),
                  file=sys.stderr)
    print("  background: %d same-gene transcript(s)%s"
          % (bg["n_background_transcripts"],
             ", FASTA %s: %d competing record(s)" % (bg["fasta"], bg["n_fasta_competitors"])
             if bg["fasta"] else ""))
    if not bg["fasta"]:
        scope = ("this gene's other transcripts"
                 if bg["n_background_transcripts"] else "the configured groups only")
        print("  NOTE: uniqueness judged against %s. A quantifier resolves fragments "
              "against the whole index, so pseudogenes, paralogues and homologous "
              "loci outside this gene are not accounted for here. Pass "
              "--background-fasta <the transcript FASTA the Salmon index was built from, "
              "without the genome decoys> for the answer that matches what the quantifier "
              "actually sees; that is the recommended way to run this command."
              % scope, file=sys.stderr)
    d = res["design"]
    print("  design: %s %dbp reads, fragments %.0f+-%.0f, depth %.0fM, TPM %.3g, n=%d"
          % ("paired" if d["paired"] else "single", d["read_length"],
             d["frag_mean"], d["frag_sd"], d["depth"] / 1e6, d["tpm"], d["n_donors"]))
    def _fc(r):
        """min resolvable |log2FC|, or why there is no usable figure."""
        v = r.get("min_resolvable_log2fc")
        if v is None or not math.isfinite(v):
            return "min |log2FC| n/a"
        return "min |log2FC| %.3f%s" % (v, "*" if r.get("beyond_linear") else "")

    incoherent = []
    for g, r in res["groups"].items():
        # unique k-mers are the class's; the stretch and the reads its best transcript's
        print("  [%s] %s: %s; class %d unique k-mers; best transcript %d bp in %d block(s), "
              "~%.0f informative reads; conditioning %.2f"
              % (r["verdict"], g, _fc(r), r["n_unique_kmers"], r["unique_length"],
                 r["n_blocks"], r["expected_informative_reads"],
                 r["conditioning_factor"]))
        coh = r.get("coherence") or {}
        if coh.get("min_jaccard") is not None and coh["min_jaccard"] < 0.05:
            incoherent.append((g, coh["min_jaccard"]))
    c = res["contrast"]
    print("  contrast %s vs %s: %s; %s, conditioning %.2f"
          % (res["primary_comparison"][0], res["primary_comparison"][1],
             _fc(c), c["verdict"], c["conditioning_factor"]))
    a_name, b_name = res["primary_comparison"][:2]
    ratio = c.get("log2_efflen_ratio")
    if ratio is not None:
        eff = c["class_mean_efflen"]
        print("  effective length: %s %.0f bp vs %s %.0f bp (mean per transcript), "
              "log2 ratio %.2f" % (a_name, eff[a_name], b_name, eff[b_name], ratio))

        def _pos(g):
            p = c["distinguishing_window_position"][g]
            if not p["n"]:
                return "%s none" % g
            return ("%s median %.2f [%.2f-%.2f], %d positions over %d transcript(s)"
                    % (g, p["median"], p["q1"], p["q3"], p["n"],
                       res["groups"][g]["n_transcripts"]))
        print("  distinguishing windows (0 = 5' end, 1 = 3' end of each transcript): "
              "%s; %s" % (_pos(a_name), _pos(b_name)))
        if c["efflen_direction_in_band"]:
            shorter, longer = (a_name, b_name) if ratio < 0 else (b_name, a_name)
            print("  NOTE: the classes differ %.1f-fold in mean effective length. Under "
                  "positional coverage skew a 5'-skewed library tends to inflate the "
                  "shorter class (%s) and a 3'-skewed one the longer (%s): in a 49-gene "
                  "simulation (Salmon, one quantifier, monotone positional skew) the 5' "
                  "direction held for 35-36 of the 39 genes with a ratio above 1.23x, and "
                  "the 3' direction for 30-32 of them, against 20 of 39 at uniform "
                  "coverage. No other quantifier, gene panel or form of skew was tested. "
                  "Measure your libraries' gene-body coverage (RSeQC geneBody_coverage.py) "
                  "to know which direction applies. This tendency assumes the classes are "
                  "ambiguous where the skew concentrates reads. When the shorter class is "
                  "itself distinguished at that end -- see its distinguishing-window "
                  "position above -- the direction can reverse; the interaction was not "
                  "measured." % (2 ** abs(ratio), shorter, longer),
                  file=sys.stderr)
    if any(r.get("beyond_linear") for r in list(res["groups"].values()) + [c]):
        print("  * first-order figure past the linearisation limit: read it as "
              "'not resolvable at this design', not as a value.", file=sys.stderr)
    for g, j in incoherent:
        print("  NOTE: class %s is not coherent -- its least similar pair of transcripts "
              "shares %.1f%% of windows. `annotate` groups by the 3' terminal-exon "
              "acceptor alone, so a class can hold transcripts that have nothing else in "
              "common; a precision figure for such a class describes a quantity nobody "
              "asked for. Review the grouping before using any number above."
              % (g, 100.0 * j), file=sys.stderr)
    noise = res["counting_noise"]
    print("  unique-read counting floor on log2 ratio: SE %.3f per donor, "
          "min resolvable |log2FC| %.3f at n=%d"
          % (noise["log2_ratio_se"], noise["min_resolvable_log2fc"], d["n_donors"]))
    print("    (precision of log2(n_a/n_b), the informative-read counts of each class's "
          "best single transcript -- that ratio equals the class ratio only when the two "
          "classes have the same informative fraction, so it can read better than the "
          "class comparison actually is. The per-class figures above are for the class "
          "totals themselves.)", file=sys.stderr)
    gt = res["gene_total"]
    if not gt["estimable"]:
        print("  PRECONDITION FAILED: the gene total is not estimable. %s shorter than "
              "window=%d, so %s no windows and an all-zero column in the compatibility "
              "system; no class total or contrast built on it is well posed. This does "
              "not depend on the grouping, the design or any threshold. Exit status 2."
              % (", ".join(gt["transcripts_without_windows"]) or "A transcript is",
                 res["window"],
                 "it has" if len(gt["transcripts_without_windows"]) == 1 else "they have"),
              file=sys.stderr)
    if d["min_log2fc"] is not None:
        past = [g for g in res["primary_comparison"][:2] if res["groups"][g]["beyond_linear"]]
        past += ["the contrast"] if c["beyond_linear"] else []
        print("  EFFECT SIZE: |log2FC| %.3f %s at this design%s"
              % (d["min_log2fc"],
                 "resolved by both class totals and the contrast"
                 if res["effect_resolvable"] else "NOT resolved",
                 "" if res["effect_resolvable"] or not past else
                 " (%s beyond the linearisation limit, relative SE > %.1f)"
                 % (", ".join(past), identifiability.LINEARISATION_LIMIT)))
    print("  VERDICT:", res["verdict"])
    for reason in res["reasons"]:
        print("    - %s" % reason, file=sys.stderr)
    if d["min_log2fc"] is None:
        print("  The verdict is a screening flag, not the exit status. It says whether each "
              "class total and the contrast lie in the row space of a sequence-derived "
              "compatibility surrogate at window=%d, for exactly the transcripts in this "
              "config. It does not say whether a quantifier will recover the comparison, "
              "and it changes with the annotation release those transcripts came from. "
              "Pass --min-log2fc <the smallest effect you need> for an exit status on the "
              "resolvable effect size." % res["window"], file=sys.stderr)
    if res["verdict"] == "not_identifiable":
        failed = [g for g in res["groups"] if not res["groups"][g].get("estimable")]
        if not res["contrast"].get("estimable"):
            failed.append("the class contrast")
        print("  Outside the row space of the compatibility surrogate at this window "
              "length: %s. That surrogate is built from sequence, not from the "
              "observation model of a sequencing run, so this is a screening verdict "
              "rather than a statement about the data: a longer --window, a different "
              "grouping, or long reads may change it."
              % (", ".join(failed) or "an estimand"), file=sys.stderr)
    return _identifiability_exit(res)


def cmd_extract(a):
    notes = []
    n = extract.run(io.load_config(a.config), a.quantdir, a.samplemap, a.cohort, a.out,
                    allow_mixed_index=a.allow_mixed_index, notes=notes)
    for note in notes:
        print("  " + note, file=sys.stderr)
    if a.json:
        _emit({"out": a.out, "n_donors": n, "cohort": a.cohort})
    else:
        print("wrote %s (n=%d donors)" % (a.out, n))
    return EXIT_OK


def cmd_stats(a):
    res = stats.run(io.load_config(a.config), a.condition, _kv(a.perdonor), a.out,
                    n_boot=a.n_boot, seed=a.seed)
    if a.json:
        _emit(res)
        return EXIT_OK
    for det in res["detail"]:
        line = ("  %-12s n=%d  %d/%d  fold=%.1fx [%.1f-%.1f]  P=%.4g (%s)"
                % (det["cohort"], det["n"], det["n_greater"], det["n"],
                   det["median_fold"], det["fold_ci"][0], det["fold_ci"][1], det["p"],
                   det["wilcoxon_method"] or "no test"))
        if det["underpowered"]:
            line += "  (floor %.4g: cannot reach 0.05)" % det["p_floor"]
        print(line)
    cn, cgt, cp, cfold = res["combined"]
    print("  %-12s n=%d  %d/%d  fold=%.1fx  P=%.4g (%s)"
          % ("POOLED", cn, cgt, cn, cfold, cp,
             res["pooled"]["wilcoxon_method"] or "no test"))
    for key in ("stouffer", "stratified_signed_rank"):
        c = res["combination"][key]
        print("  %-12s k=%d cohorts  P=%.4g" % (key.upper(), c["k"], c["p"]))
    print("  headline combination: %s. POOLED above is reported for continuity only: "
          "it pools donors across independent cohorts and so ignores the cohort factor."
          % res["headline_combination"])
    print("wrote %s.{png,pdf,svg} + %s_stats.csv" % (a.out, a.out))
    return EXIT_OK


def cmd_qc(a):
    rows = contamination.run(io.load_config(a.config, need_groups=False), _kv(a.markers),
                             _kv(a.target), a.out)
    if a.json:
        _emit([{"cohort": r[0], "n": r[1], "rho": r[2], "p": r[3], "ratio": r[4]}
               for r in rows])
        return EXIT_OK
    for name, n, rho, p, ratio in rows:
        print("  %-12s n=%d  rho=%+.3f  P=%.3f  contam/tissue=%.3f" % (name, n, rho, p, ratio))
    print("wrote %s.{png,pdf,svg} + %s_scores.csv" % (a.out, a.out))
    return EXIT_OK


def cmd_selftest(a):
    from . import _selftest
    if a.json:
        res = _selftest.result()
        _emit(res)
        return EXIT_OK if res["ok"] else 1
    return _selftest.main()


def build_parser():
    p = argparse.ArgumentParser(prog="isoform-dominance",
                                description="Isoform-usage quantification and discrimination from bulk RNA-seq.")
    p.add_argument("--version", action="version", version="isoform-dominance %s" % __version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def _json(sp):
        sp.add_argument("--json", action="store_true",
                        help="emit the full result as JSON on stdout")
        return sp

    def _net(sp):
        sp.add_argument("--retries", type=int, default=ensembl.DEFAULT_RETRIES,
                        help="retries per Ensembl request after the first attempt, on HTTP "
                             "429/500/502/503/504, connection errors and read timeouts "
                             "(default: %(default)s)")
        sp.add_argument("--retry-wait", type=float, default=ensembl.DEFAULT_RETRY_WAIT,
                        help="seconds before the first retry, doubling each time; an HTTP "
                             "429 waits for its Retry-After instead (default: %(default)s)")
        sp.add_argument("--ensembl-release", type=_release, default=None, metavar="N",
                        help="Ensembl release to fetch from (default: the current one at "
                             "rest.ensembl.org). Earlier releases are read from Ensembl's "
                             "REST archive -- release 110 is GENCODE 44 -- and the release "
                             "the server reports is checked against N before anything is "
                             "fetched")
        return sp

    s = _net(_json(sub.add_parser("annotate",
                                  help="gene symbol -> proposed isoform groups (Ensembl)")))
    s.add_argument("--gene", required=True); s.add_argument("--species", default="homo_sapiens")
    s.add_argument("--gene-id", metavar="ENSG",
                   help="the Ensembl gene to propose groups from, when the symbol names more "
                        "than one on the reference chromosomes")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_annotate)

    s = _net(_json(sub.add_parser(
        "identifiability", aliases=["identify"],
        help="are the classes measurable by short reads, and how precisely?")))
    s.add_argument("--config", required=True)
    s.add_argument("--k", type=int, default=identifiability.DEFAULT_K,
                   help="k-mer length; sets the default window when --window is not given, "
                        "and nothing else (default: %(default)s)")
    s.add_argument("--window", type=int, default=None,
                   help="window length for the compatibility system (default: k). A "
                        "different window gives a different system, not a uniformly "
                        "sharper one: the rank is not monotone in it.")
    s.add_argument("--save-inputs", metavar="FILE",
                   help="write the sequence this run used -- the configured transcripts, "
                        "the gene background and the Ensembl release they came from -- to "
                        "FILE, so the run can be repeated with --inputs and no network after "
                        "the release's REST archive is gone. A --background-fasta is recorded "
                        "by path and SHA-256, not copied")
    s.add_argument("--inputs", metavar="FILE",
                   help="run on sequence saved by --save-inputs, with no request at all; "
                        "replaces --sequences, --background-sequences and --ensembl-release")
    s.add_argument("--sequences", help="optional JSON {transcript_id: cdna} (offline)")
    s.add_argument("--background-sequences", help="optional JSON {transcript_id: cdna} of background transcripts")
    s.add_argument("--background-fasta",
                   help="FASTA (optionally gzipped) to judge uniqueness against: the "
                        "transcript FASTA the Salmon index was built from, without the "
                        "genome decoys (or pass --decoys). Each record that shares a window "
                        "with a configured transcript is a column of the compatibility "
                        "system")
    s.add_argument("--decoys", metavar="FILE",
                   help="Salmon's decoys.txt, one record name per line, when "
                        "--background-fasta is a decoy-aware index's gentrome: those records "
                        "are skipped, as Salmon sets aside only reads that map better to a "
                        "decoy than to any transcript")
    s.add_argument("--keep-duplicates", action="store_true",
                   help="the Salmon index was built with --keepDuplicates: count a "
                        "--background-fasta record with a configured transcript's sequence as "
                        "competing sequence (by default it is not, as Salmon's default index "
                        "keeps one of identical sequences)")
    s.add_argument("--no-gene-background", action="store_true",
                   help="do not fetch the gene's other transcripts as background (v2.1 behaviour)")
    s.add_argument("--strand-aware", action="store_true",
                   help="do not fold k-mers to their canonical form (v2.1 behaviour)")
    s.add_argument("--read-length", type=int, default=identifiability.DEFAULT_READ_LENGTH)
    s.add_argument("--frag-mean", type=float, default=identifiability.DEFAULT_FRAG_MEAN)
    s.add_argument("--frag-sd", type=float, default=identifiability.DEFAULT_FRAG_SD)
    s.add_argument("--single-end", action="store_true")
    s.add_argument("--depth", type=float, default=identifiability.DEFAULT_DEPTH,
                   help="mapped fragments per library")
    s.add_argument("--tpm", type=float, default=identifiability.DEFAULT_TPM,
                   help="TPM given to every transcript of the gene, background included, "
                        "for the read model and the Poisson weights (default: %(default)s)")
    s.add_argument("--donors", type=int, default=1)
    s.add_argument("--min-log2fc", type=float, default=None,
                   help="smallest |log2 fold change| you need both class totals and the "
                        "contrast to resolve; when given it replaces --tau in the verdict and sets the "
                        "exit status (3 when not resolved), which is the recommended way "
                        "to run this (--tau has no calibrated value: over a 49-gene survey "
                        "the median gene sat at conditioning 65 against the default 10)")
    s.add_argument("--tau", type=float, default=identifiability.DEFAULT_CONDITIONING_TAU,
                   help="structural conditioning factor above which a class is only weakly "
                        "identifiable; a diagnostic that does not affect the exit status")
    s.add_argument("--min-informative-reads", type=float,
                   default=identifiability.DEFAULT_MIN_INFORMATIVE_READS)
    s.set_defaults(func=cmd_identifiability)

    s = _json(sub.add_parser("extract", help="quant.sf -> per-donor isoform-group TPM"))
    for x in ("config", "quantdir", "samplemap", "cohort", "out"):
        s.add_argument("--" + x, required=True)
    s.add_argument("--allow-mixed-index", action="store_true",
                   help="combine donors quantified against different Salmon indexes "
                        "(different index_seq_hash in aux_info/meta_info.json) instead of "
                        "stopping; the sidecar <out>.index.json records which is which")
    s.set_defaults(func=cmd_extract)

    s = _json(sub.add_parser("stats", help="paired Wilcoxon, cohort combination + figure"))
    s.add_argument("--config", required=True); s.add_argument("--condition", default="control")
    s.add_argument("--perdonor", action="append", required=True, help="NAME=perdonor.csv (repeatable)")
    s.add_argument("--n-boot", type=int, default=stats.DEFAULT_N_BOOT,
                   help="bootstrap replicates for the fold-change interval (0 disables)")
    s.add_argument("--seed", type=int, default=stats.DEFAULT_SEED)
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_stats)

    s = _json(sub.add_parser("qc", help="contamination control"))
    s.add_argument("--config", required=True)
    s.add_argument("--markers", action="append", required=True, help="NAME=marker_tpm.csv")
    s.add_argument("--target", action="append", required=True, help="NAME=perdonor.csv")
    s.add_argument("--out", required=True); s.set_defaults(func=cmd_qc)

    s = _json(sub.add_parser("selftest", help="run the download-free reproducibility test"))
    s.set_defaults(func=cmd_selftest)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args) or 0
    except io.InputError as e:              # a file named on the command line is unusable
        print("error: %s" % e, file=sys.stderr)
        return 1
    except OSError as e:                    # ... or cannot be opened: one line, not a trace
        if e.filename is None:
            raise
        print("error: %s: %s" % (e.strerror or e, e.filename), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
