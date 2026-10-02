"""Shared IO helpers: config, sample-map and saved-inputs loading."""
import csv
import datetime
import gzip
import hashlib
import json
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
    only ``contamination_qc``)."""
    cfg = load_json(path, "config")
    if not isinstance(cfg, dict):
        raise InputError("config %s is not a JSON object" % path)
    if need_groups and not isinstance(cfg.get("groups"), dict):
        raise InputError("config %s has no \"groups\" object mapping group names to "
                         "transcript ids; `annotate` writes one" % path)
    return cfg


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


def primary_pair(config):
    """The two groups named in ``primary_comparison`` (for the paired comparison)."""
    groups = config["groups"]
    pc = config.get("primary_comparison", list(groups)[:2])
    if len(pc) < 2:
        raise ValueError(
            "primary_comparison must name two isoform groups for a paired "
            "comparison; got %r. `annotate` proposes a pair only when an "
            "alternative terminal-exon cluster exists; edit the config to define "
            "two groups." % (pc,))
    return pc[0], pc[1]


def save_inputs(path, captured, config, release, version, background_fasta=None,
                decoys=None, max_window_records=None):
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
           "sequence_sources": captured.get("sequence_sources"),
           "sequences": captured["sequences"],
           "background_sequences": captured["background_sequences"]}
    if background_fasta and max_window_records is not None:
        doc["analysis"]["max_window_records"] = int(max_window_records)
    if background_fasta or doc["analysis"]["gene_background"]:
        # without a background there was no copy to count or not, and a rerun with one
        # must not be told that the setting changed
        doc["analysis"]["keep_duplicates"] = bool(captured.get("keep_duplicates"))
    with open(path, "w") as f:
        json.dump(doc, f)
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
    for key in ("background_fasta", "decoys"):
        rec = doc.get(key)
        if rec is not None and not (isinstance(rec, dict)
                                    and isinstance(rec.get("path"), str)
                                    and isinstance(rec.get("sha256"), str)):
            raise InputError("saved inputs %s: %s is %r, expected null or an object with a "
                             "path and a sha256" % (path, key, rec))
    return doc
