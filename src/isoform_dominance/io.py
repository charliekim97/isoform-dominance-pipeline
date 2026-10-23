"""Shared IO helpers: config, sample-map and saved-inputs loading."""
import collections.abc
import csv
import datetime
import gzip
import hashlib
import json
import math
import numbers
import os

#: The ``format`` field of a file written by :func:`save_inputs`.
INPUTS_FORMAT = "isoform-dominance/inputs/1"


class InputError(ValueError):
    """A file named on the command line that cannot be used as it is."""


def load_json(path, what):
    """Parse the JSON file ``path``; :class:`InputError` naming ``what`` if it is not JSON
    in UTF-8, the only encoding JSON may be exchanged in (RFC 8259)."""
    with open(path, encoding="utf-8") as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            raise InputError("%s %s is not valid JSON (%s)" % (what, path, e)) from e
        except UnicodeDecodeError as e:
            raise InputError("%s %s is not UTF-8 text (byte 0x%02x at offset %d)"
                             % (what, path, e.object[e.start], e.start)) from e


def finite_json(obj):
    """``obj`` with every number that is not finite -- NaN, inf, -inf -- replaced by None,
    through dicts, lists and tuples (a tuple becomes a list, as JSON has it).

    RFC 8259 has no such values.  Python's ``json`` writes them as ``NaN`` and
    ``Infinity`` and reads them back, but a standard parser refuses the whole document,
    so everything this package writes goes through here and :func:`dump_json` or
    ``allow_nan=False``.  Where a null needs a reason, the report carries a field beside
    it that says it (``estimable``, ``finite_se``, ``defined``, ``fold_ci_method``)."""
    if isinstance(obj, dict):
        return {k: finite_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [finite_json(v) for v in obj]
    if isinstance(obj, numbers.Real) and not isinstance(obj, numbers.Integral):
        return float(obj) if math.isfinite(obj) else None
    return obj


def dump_json(obj, path, **kw):
    """Write ``obj`` to ``path`` as standard JSON: a number that is not finite raises
    ``ValueError`` rather than being written as a token no standard parser reads."""
    with open(path, "w") as f:
        json.dump(obj, f, allow_nan=False, **kw)


def is_gzip(path):
    """Is ``path`` gzipped?  Told by its first two bytes (``1f 8b``), not by its name: a
    plain file called ``.gz`` and a gzipped one called ``.GZ`` or nothing at all are both
    common."""
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def open_text(path):
    """Open a possibly gzipped text file for reading (:func:`is_gzip`)."""
    return gzip.open(path, "rt") if is_gzip(path) else open(path)


def open_bytes(path):
    """Open a possibly gzipped file for reading bytes (:func:`is_gzip`)."""
    return gzip.open(path, "rb") if is_gzip(path) else open(path, "rb")


def load_config(path, need_groups=True):
    """A config written by ``annotate`` or by hand: a JSON object, with a ``groups`` object
    mapping group names to transcript ids unless ``need_groups`` is false (``qc`` reads
    only ``contamination_qc``).  With ``need_groups`` the groups and ``primary_comparison``
    are held to :func:`groups_shape_error`, and a config that fails it is an
    :class:`InputError` naming the file and the place."""
    cfg = load_json(path, "config")
    if not isinstance(cfg, dict):
        raise InputError("config %s is not a JSON object" % path)
    if need_groups and not isinstance(cfg.get("groups"), dict):
        raise InputError("config %s has no \"groups\" object mapping group names to "
                         "transcript ids; `annotate` writes one" % path)
    if need_groups:
        wrong = groups_shape_error(cfg)
        if wrong:
            raise InputError("config %s: %s" % (path, wrong))
    return cfg


def _shown(value):
    """``value`` as the config would write it, cut short when long."""
    try:
        text = json.dumps(value, ensure_ascii=False)
    except (TypeError, ValueError):
        text = repr(value)
    return text if len(text) <= 60 else text[:57] + "..."


def _a_list(value):
    """Whether ``value`` reads as a list: something with a length that can be iterated, and
    is not a string, bytes or a mapping.  A JSON array; from Python, a tuple, a set or an
    array too, as 2.8.0 took them."""
    if isinstance(value, (str, bytes, bytearray, collections.abc.Mapping)):
        return False
    try:
        len(value)
        iter(value)
    except TypeError:                       # a number, None, an iterator, a 0-d array
        return False
    return True


#: what is wrong with a ``primary_comparison`` given, from Python, as a set
_A_SET = ("primary_comparison is %s, a set: it has no first and second group; give them in "
          "a list")


def groups_shape_error(config):
    """What is wrong with the shape of ``config``'s classes and comparison, as a sentence
    that starts with the place (``groups``, ``group "<name>"`` or ``primary_comparison``),
    or None.

    ``groups`` maps one group name or more to a non-empty list of transcript ids, each a
    non-empty string with no white space; ``primary_comparison``, when present, is a list of
    group names, each a string.  From Python a list may be any collection that is not a
    string or a mapping -- a tuple, a set, an array -- except that a comparison must have an
    order, and so is no set.  Through 2.8.0 nothing checked this: a class given as a string
    was read one character at a time as transcript ids, an empty class gave a verdict, an id
    with a space was read as it stood (asked of Ensembl, and found in no ``quant.sf``, unless
    a version followed it), and a comparison given as a string was read as its letters.
    That the comparison names two groups, both the config's and not one twice, is checked
    where the comparison is used: by :func:`primary_pair` (``stats``) and by
    ``identifiability``; ``extract`` writes no fraction column then."""
    groups = config.get("groups")
    if not isinstance(groups, collections.abc.Mapping):
        return "groups is %s, not an object mapping group names to transcript ids" % (
            _shown(groups))
    if not groups:
        return ("groups is empty: a config names one group or more, each with a list of "
                "transcript ids")
    for name, ids in groups.items():
        place = "group %s" % _shown(name)
        if not _a_list(ids):
            return "%s is %s, not a list of transcript ids" % (place, _shown(ids))
        if len(ids) == 0:                   # an array has no truth value
            return "%s lists no transcript; a group lists one or more transcript ids" % place
        for t in ids:
            if not (isinstance(t, str) and t and not any(c.isspace() for c in t)):
                return ("%s lists %s, which is not a transcript id: a transcript id is a "
                        "non-empty string with no white space" % (place, _shown(t)))
    if "primary_comparison" not in config:
        return None
    pc = config["primary_comparison"]
    if isinstance(pc, collections.abc.Set):
        return _A_SET % _shown(pc)
    if not _a_list(pc):
        return "primary_comparison is %s, not a list of group names" % _shown(pc)
    for g in pc:
        if not isinstance(g, str):
            return "primary_comparison lists %s, which is not a group name" % _shown(g)
    return None


def file_sha256(path, chunk=1 << 20):
    """Hex SHA-256 of the file's bytes, read in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def read_decoys(path):
    """The record names in Salmon's ``decoys.txt``: one per line, its first word, any
    leading ``>`` dropped; blank lines skipped.  A file with none is an :class:`InputError`.
    """
    names = []
    with open_text(path) as f:
        for line in f:
            word = line.strip().lstrip(">").split()
            if word:
                names.append(word[0])
    if not names:
        raise InputError("--decoys %s names no record" % path)
    return names


def transcript_to_group(groups):
    """{group: [ENST,...]} -> {ENST(no version): group}."""
    m = {}
    for g, txs in groups.items():
        for t in txs:
            m[t.split(".")[0]] = g
    return m


def shared_transcripts(groups):
    """``{transcript (no version): [group, ...]}`` for each transcript more than one group
    names, groups in config order.  :func:`transcript_to_group` gives such a transcript to
    the group listed last."""
    seen = {}
    for g, txs in groups.items():
        for t in txs:
            names = seen.setdefault(t.split(".")[0], [])
            if g not in names:
                names.append(g)
    return {t: gs for t, gs in seen.items() if len(gs) > 1}


def load_sample_map(path):
    """CSV with columns donor,condition[,SRR] -> {donor: condition}."""
    cond = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        if "donor" not in (reader.fieldnames or []):
            raise InputError("sample map %s has no \"donor\" column (its columns: %s); "
                             "expected donor,condition[,SRR]"
                             % (path, ", ".join(reader.fieldnames or []) or "none"))
        for r in reader:
            cond[r["donor"]] = r.get("condition", "NA")
    return cond


def _comparison(config):
    """``config``'s ``primary_comparison`` as a list -- without one, the config's first two
    groups -- for the code that reads it in order.  ``ValueError`` for a set, from Python:
    it has no first and second group (2.8.0 raised ``TypeError`` there, and a ``list()``
    of it would take Python's order of the moment)."""
    pc = config.get("primary_comparison", list(config["groups"])[:2])
    if isinstance(pc, collections.abc.Set):
        raise ValueError(_A_SET % _shown(pc))
    return list(pc)


def primary_pair(config):
    """The two groups named in ``primary_comparison`` (for the paired comparison); without
    one, the config's first two groups.  ``ValueError`` when the comparison does not name
    exactly two groups, names a group the config does not have, or names one group twice,
    as ``identifiability`` refuses such a comparison, and for a set; ``stats`` says it in
    one line naming the config."""
    groups = config["groups"]
    pc = _comparison(config)
    if len(pc) < 2:
        raise ValueError(
            "primary_comparison must name two isoform groups for a paired "
            "comparison; got %r. `annotate` proposes a pair only when an "
            "alternative terminal-exon cluster exists; edit the config to define "
            "two groups." % (pc,))
    if len(pc) > 2:
        # through 2.8.0 the first two were compared, and the others said nothing
        raise ValueError("primary_comparison must name two isoform groups for a paired "
                         "comparison; got %r, %d groups" % (pc, len(pc)))
    missing = list(dict.fromkeys(g for g in pc if g not in groups))
    if missing:
        # through 2.8.0 the per-donor table refused it, for a column the config never named
        raise ValueError("primary_comparison names group(s) not in config['groups']: %r"
                         % (missing,))
    if pc[0] == pc[1]:
        # each donor's class against itself: every ratio 1, and P NaN
        raise ValueError("primary_comparison compares group %r with itself; name two "
                         "different groups" % (pc[0],))
    return pc[0], pc[1]


def save_inputs(path, captured, config, release, version, background_fasta=None,
                decoys=None, max_window_records=None, annotation_source=None):
    """Write the sequence an ``identifiability`` run used, so it can be repeated offline.

    ``captured`` is what :func:`isoform_dominance.identifiability.analyze` put in
    ``inputs_out``; ``release`` is the Ensembl release that sequence came from, or None
    when it was supplied from an unknown one.  Ensembl keeps REST archives for recent
    releases only, and stopped publishing new releases on the REST API after 116, so
    this file, not the archive, is what keeps a verdict reproducible.

    A ``background_fasta`` is not copied -- a transcriptome FASTA runs to hundreds of
    megabytes -- but its path, size and SHA-256 are recorded, so that a rerun from this
    file can say when it is given no FASTA, or another one; and so are a ``decoys`` file's.

    ``analysis`` records the k, window and k-mer convention the run used.  None of the
    three is in the config, so without them nothing can tell a rerun from this file that
    it is building a different compatibility system on the same sequence.  It records
    ``keep_duplicates`` too, when the run had a background, gene or FASTA: whether a
    background sequence identical to a configured transcript or to another was counted;
    and, with a FASTA, ``max_window_records``, which decides the FASTA records that join
    the system; and ``gene_background``: whether the run judged
    uniqueness against the gene's other transcripts.  A rerun under another grouping moves a
    transcript the config no longer names into that background only if there was one.

    ``sequence_sources`` gives each sequence's origin, ``"supplied"`` or
    ``"fetched:<release>"``.  ``release`` is recorded only when every sequence was fetched
    from it: a run on partly supplied sequence did not use one release, and saying it did
    would date the supplied sequence to a release it may not come from.

    ``gene_id`` is the Ensembl gene the gene background was fetched as -- or, when none was
    fetched, the config's ``gene_id`` -- so that a rerun can be refused a config of another
    gene of the same name.

    A run on local files (``identifiability --gtf --transcripts-fasta``) records them as
    ``annotation_source`` -- names, sizes, SHA-256 and what the GTF's header says -- and
    ``release`` is theirs; each sequence's source is ``"file:<release>"``.  ``versions``
    gives each sequence's versioned id where the run knew it, so that a rerun can still
    say when ``--background-fasta`` holds a configured transcript at another version.
    """
    def _file(p):
        return p and {"path": str(p), "bytes": os.path.getsize(p), "sha256": file_sha256(p)}
    fasta = _file(background_fasta) or None
    doc = {"format": INPUTS_FORMAT, "package_version": version,
           "saved": datetime.date.today().isoformat(),
           "gene": config.get("gene"), "species": config.get("species", "homo_sapiens"),
           "gene_id": captured.get("gene_id") or config.get("gene_id"),
           "ensembl_release": release,
           "config_ensembl_release": config.get("ensembl_release"),
           "analysis": dict({key: captured[key] for key in ("k", "window", "canonical")},
                            gene_background=bool(captured.get("gene_background",
                                                              captured["background_sequences"]))),
           "background_fasta": fasta,
           "decoys": _file(decoys) or None,
           "annotation_source": annotation_source,
           # which version of each sequence, where known: a rerun can still tell an index
           # of another release
           "versions": captured.get("versions") or {},
           "sequence_sources": captured.get("sequence_sources"),
           "sequences": captured["sequences"],
           "background_sequences": captured["background_sequences"]}
    if background_fasta and max_window_records is not None:
        doc["analysis"]["max_window_records"] = int(max_window_records)
    if background_fasta or doc["analysis"]["gene_background"]:
        # without a background there was no copy to count or not, and a rerun with one
        # must not be told that the setting changed
        doc["analysis"]["keep_duplicates"] = bool(captured.get("keep_duplicates"))
    dump_json(doc, path)
    return doc


def load_inputs(path):
    """Read a file written by :func:`save_inputs`; :class:`InputError` if it is not one.

    Every field the caller goes on to read is checked here, so that a hand-edited or
    truncated file is one line naming the file rather than a ``KeyError`` or a
    ``TypeError`` from somewhere inside the analysis.
    """
    doc = load_json(path, "saved inputs")
    fmt = doc.get("format") if isinstance(doc, dict) else None
    if fmt != INPUTS_FORMAT:
        raise InputError("%s is not a saved-inputs file (format %r, expected %r); write one "
                         "with `identifiability --save-inputs`" % (path, fmt, INPUTS_FORMAT))
    for key in ("sequences", "background_sequences"):
        seqs = doc.get(key)
        if not (isinstance(seqs, dict)
                and all(isinstance(k, str) and isinstance(v, str) for k, v in seqs.items())):
            raise InputError("saved inputs %s: %s is not an object of transcript id to "
                             "sequence" % (path, key))
    if "ensembl_release" not in doc:
        raise InputError("saved inputs %s records no ensembl_release; write the file with "
                         "`identifiability --save-inputs`" % path)
    if not isinstance(doc["ensembl_release"], (int, type(None))) \
            or isinstance(doc["ensembl_release"], bool):
        raise InputError("saved inputs %s: ensembl_release is %r, expected a release number "
                         "or null" % (path, doc["ensembl_release"]))
    analysis = doc.get("analysis")
    if analysis is None:
        raise InputError("saved inputs %s records no analysis (k, window, canonical); it "
                         "was written by a version before 2.4, and which compatibility "
                         "system it was saved from cannot be recovered from it" % path)
    if not (isinstance(analysis, dict)
            and all(isinstance(analysis.get(key), int) and not isinstance(analysis.get(key), bool)
                    for key in ("k", "window"))
            and isinstance(analysis.get("canonical"), bool)
            and isinstance(analysis.get("keep_duplicates", False), bool)
            and isinstance(analysis.get("max_window_records", 0), int)
            and not isinstance(analysis.get("max_window_records", 0), bool)
            and isinstance(analysis.get("gene_background", False), bool)):
        raise InputError("saved inputs %s: analysis is %r, expected k and window as "
                         "integers, and canonical, keep_duplicates and gene_background as "
                         "booleans (and max_window_records an integer)" % (path, analysis))
    if not isinstance(doc.get("gene_id"), (str, type(None))):
        raise InputError("saved inputs %s: gene_id is %r, expected an Ensembl gene id or null"
                         % (path, doc["gene_id"]))
    sources = doc.get("sequence_sources")
    if sources is not None and not (isinstance(sources, dict) and all(
            isinstance(k, str) and isinstance(v, str) for k, v in sources.items())):
        raise InputError("saved inputs %s: sequence_sources is %r, expected null or an object "
                         "of transcript id to \"supplied\" or \"fetched:<release>\""
                         % (path, sources))
    src = doc.get("annotation_source")
    if src is not None and not (
            isinstance(src, dict) and isinstance(src.get("kind"), str)
            and isinstance(src.get("transcripts_fasta"), dict)
            and isinstance(src["transcripts_fasta"].get("file"), str)
            and isinstance(src.get("gtf"), (dict, type(None)))):
        raise InputError("saved inputs %s: annotation_source is %r, expected null or the "
                         "files a run read (kind, gtf, transcripts_fasta)" % (path, src))
    versions = doc.get("versions")
    if versions is not None and not (isinstance(versions, dict) and all(
            isinstance(k, str) and isinstance(v, str) for k, v in versions.items())):
        raise InputError("saved inputs %s: versions is %r, expected null or an object of "
                         "transcript id to versioned id" % (path, versions))
    for key in ("background_fasta", "decoys"):
        rec = doc.get(key)
        if rec is not None and not (isinstance(rec, dict)
                                    and isinstance(rec.get("path"), str)
                                    and isinstance(rec.get("sha256"), str)):
            raise InputError("saved inputs %s: %s is %r, expected null or an object with a "
                             "path and a sha256" % (path, key, rec))
    return doc
