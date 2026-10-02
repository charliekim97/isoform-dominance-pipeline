"""Index scope: same-name copies on non-reference regions, and which index quantified whom.

From GENCODE release 48, ``gencode.vX.transcripts.fa.gz`` holds the transcripts on
scaffolds, patches and alternate loci as well as on the reference chromosomes; releases
44-47 hold the reference chromosomes only.  The change is in the release-48 entry of the
changelog in GENCODE's FTP ``_README.TXT``, whose file description still calls the file
reference-chromosome only.  Ensembl's ``cdna.all`` holds both too (checked for release 116).

A gene with a copy on such a region is in the FASTA twice, under its own name and another
gene id.  An index built from that FASTA gives the quantifier a second target for the
gene's reads: Salmon folds an identical copy into the reference transcript but splits
reads with a copy that differs, and kallisto splits them either way.  ``extract`` sums only
the transcript ids in the config, so the copy's share is missing from the class total.

A same-name gene on a reference chromosome is not such a copy: an index built from
reference chromosomes keeps it.  An Ensembl header names the region, and such a gene is
not reported.  A GENCODE header does not, so there it is reported too, as a possible copy:
the chrY copy of a pseudoautosomal gene such as CD99 or SHOX, which GENCODE 44, 48 and 50
give its own gene id, or a distinct gene sharing the name, such as HERC3.

The copies are found from the FASTA headers, which describe exactly the file the index was
built from, and not from Ensembl REST: ``xrefs/symbol/homo_sapiens/SMN1`` also returns
SMN2's gene ids.
"""
import json
import os

from .io import InputError, open_text

#: The ``format`` field of the ``<out>.index.json`` sidecar :mod:`extract` writes.
INDEX_FORMAT = "isoform-dominance/index/1"

#: What is kept of a Salmon ``aux_info/meta_info.json``.  The decoy fields matter because
#: two indexes of one transcriptome, one with decoys and one without, share
#: ``index_seq_hash``; Salmon hashes the decoys separately (the SHA-256 of nothing,
#: ``e3b0c442...``, when there are none).  An older Salmon writes no decoy fields.
META_KEYS = ("salmon_version", "index_seq_hash", "index_name_hash", "num_valid_targets",
             "keep_duplicates", "index_decoy_seq_hash", "num_decoy_targets")

#: The fields that say which index quantified a donor; each must be a string when present.
_HASH_KEYS = ("index_seq_hash", "index_decoy_seq_hash")

_ENSEMBL_REGION_KINDS = ("chromosome:", "scaffold:", "primary_assembly:")

#: Ensembl's names for the GRCh38 reference chromosomes, as a cDNA header gives them
#: (``chromosome:GRCh38:<name>:...``): the regions the README's filter keeps.
REFERENCE_REGIONS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])


class MixedIndexError(InputError):
    """Donors in one cohort were quantified against different indexes."""


class NoConfiguredTranscripts(InputError):
    """No donor's quant.sf has any transcript the config names."""


def _unversioned(ident):
    """``ENSG00000182378.15_PAR_Y`` -> ``ENSG00000182378``: the version, and the PAR_Y
    suffix older GENCODE releases put after it, name the same gene."""
    return ident.split(".")[0]


def parse_header(head):
    """A FASTA header (without ``>``), or a quant.sf name, as a record; None if the header
    is in neither the GENCODE nor the Ensembl cDNA format.

    GENCODE: ``transcript|gene|havana gene|havana transcript|transcript name|gene name|...``.
    Ensembl: ``transcript cdna <kind>:GRCh38:<region>:... gene:<id> ... gene_symbol:<name>``.
    The record is ``{"transcript", "gene_id", "gene_name", "region"}``, ids unversioned;
    ``region`` is None for GENCODE, whose header does not carry it.
    """
    head = head.strip()
    fields = head.split("|")
    if len(fields) >= 6 and fields[0] and fields[1] and fields[5]:
        return {"transcript": _unversioned(fields[0]), "gene_id": _unversioned(fields[1]),
                "gene_name": fields[5], "region": None}
    tokens = head.split()
    if not tokens:
        return None
    gene = name = region = None
    for tok in tokens[1:]:
        if tok.startswith("gene:"):
            gene = tok[len("gene:"):]
        elif tok.startswith("gene_symbol:"):
            name = tok[len("gene_symbol:"):]
        elif region is None and tok.startswith(_ENSEMBL_REGION_KINDS):
            parts = tok.split(":")
            region = parts[2] if len(parts) > 2 else None
    if not (gene and name):
        return None
    return {"transcript": _unversioned(tokens[0]), "gene_id": _unversioned(gene),
            "gene_name": name, "region": region}


def same_name_copies(records, target_ids, gene_names=()):
    """Genes that share a name with the target transcripts' gene but not its id.

    ``records`` are :func:`parse_header` results (None is skipped); ``target_ids`` are the
    configured transcript ids.  The target gene is whatever gene id and name the records
    give those transcripts, so a config that names the gene differently, or not at all,
    does not change the answer; ``gene_names`` adds names to match.  Nothing is reported
    when none of the target transcripts is among the records: the reference gene cannot
    then be told from its copies.

    A record on a reference chromosome (:data:`REFERENCE_REGIONS`; only an Ensembl header
    says) is not a copy.  A GENCODE header names no region, so there every same-name gene
    id is returned, with ``region`` None: each may be a copy, or a same-name gene on a
    reference chromosome.

    Returns ``[{"gene_id", "gene_name", "region", "transcripts"}]`` sorted by gene id.
    """
    records = [r for r in records if r]
    targets = {_unversioned(t) for t in target_ids}
    ref = [r for r in records if r["transcript"] in targets]
    if not ref:
        return []
    ref_genes = {r["gene_id"] for r in ref}
    names = {r["gene_name"] for r in ref} | set(gene_names or ())
    copies = {}
    for r in records:
        if (r["gene_name"] in names and r["gene_id"] not in ref_genes
                and r["region"] not in REFERENCE_REGIONS):
            c = copies.setdefault(r["gene_id"], {"gene_id": r["gene_id"],
                                                 "gene_name": r["gene_name"],
                                                 "region": r["region"], "transcripts": set()})
            c["transcripts"].add(r["transcript"])
    return [dict(c, transcripts=sorted(c["transcripts"]))
            for _, c in sorted(copies.items())]


def _relevant(heads, targets, names, versions=None):
    """The parsed ``heads`` that are a target transcript or carry one of ``names``; and in
    ``versions``, each target's versioned ids, the header's first ``|``-field or word."""
    out = []
    for head in heads:
        if versions is not None:
            first = head.strip().split("|")[0].split()
            if first and _unversioned(first[0]) in targets:
                versions.setdefault(_unversioned(first[0]), set()).add(first[0])
        r = parse_header(head)
        if r and (r["transcript"] in targets or r["gene_name"] in names):
            out.append(r)
    return out


def copies_in(read_heads, target_ids, gene_names=(), versions=None):
    """:func:`same_name_copies` over the headers ``read_heads()`` yields, holding only the
    records that can matter rather than every header of a whole-transcriptome FASTA.

    ``read_heads`` is called again, once, when the name the headers give the target
    transcripts is not among ``gene_names``: a copy can come before its reference gene.
    ``versions``, a dict, receives the versioned ids the headers give each target
    transcript (``{id: {versioned id, ...}}``).
    """
    targets = {_unversioned(t) for t in target_ids}
    names = {n for n in gene_names or () if n}
    recs = _relevant(read_heads(), targets, names, versions)
    found = {r["gene_name"] for r in recs if r["transcript"] in targets}
    if found - names:
        names |= found
        recs = _relevant(read_heads(), targets, names)
    return same_name_copies(recs, targets, names)


def fasta_copies(path, target_ids, gene_names=(), versions=None):
    """:func:`copies_in` for a plain or gzipped FASTA (told by its magic bytes)."""
    def heads():
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    yield line[1:]
    return copies_in(heads, target_ids, gene_names, versions)


def without_identical(copies, identical):
    """:func:`same_name_copies`' result less the records ``identical`` names.

    ``identical`` maps a FASTA record id to the configured transcript whose sequence it
    is.  Salmon's default index keeps one of identical sequences, so such a record takes no
    read from the configured transcript and is no reason to warn that reads are split; a
    copy left with no other record is dropped.
    """
    identical = {_unversioned(r) for r in (identical or ())}
    out = []
    for c in copies:
        rest = [t for t in c["transcripts"] if t not in identical]
        if rest:
            out.append(dict(c, transcripts=rest))
    return out


def copy_warning(copies, source):
    """The one-paragraph warning for :func:`same_name_copies`' result, or None.

    A copy with a region (an Ensembl header) is on a non-reference region.  One without (a
    GENCODE header) may instead be a same-name gene on a reference chromosome, and the
    warning says so rather than calling it a copy.
    """
    if not copies:
        return None
    placed = [c for c in copies if c["region"]]
    unplaced = [c for c in copies if not c["region"]]
    names = "/".join(sorted({c["gene_name"] for c in copies}))

    def listed(cs):           # the transcript ids are in the JSON; here, how many
        return "; ".join("%s%s, %d transcript(s)"
                         % (c["gene_id"], " on %s" % c["region"] if c["region"] else "",
                            len(c["transcripts"])) for c in cs)

    def n_tx(cs):
        return sum(len(c["transcripts"]) for c in cs)

    parts = []
    if placed:
        parts.append("%s has %d transcript(s) named %s on non-reference regions, under %d "
                     "gene id(s) other than the configured transcripts' (%s)."
                     % (source, n_tx(placed), names, len(placed), listed(placed)))
    if unplaced:
        parts.append("%s has %d transcript(s) named %s under %d gene id(s) other than the "
                     "configured transcripts' (%s), which may be copies on scaffolds, "
                     "patches or alternate loci."
                     % (source, n_tx(unplaced), names, len(unplaced), listed(unplaced)))
    parts.append("If the Salmon index was built from it, reads are split between the "
                 "configured transcripts and such copies, which `extract` does not count, "
                 "and the class totals are biased. Build the index from "
                 "reference-chromosome transcripts only (see the Full workflow, step 0).")
    if unplaced:
        parts.append("These headers do not say where a gene lies (GENCODE's never do), so "
                     "a gene id listed without a region may instead be on a reference "
                     "chromosome -- the chrY copy of a pseudoautosomal gene, which recent "
                     "GENCODE releases give its own gene id, or a different gene that shares "
                     "the name. That is not an index-scope problem: an index rebuilt from "
                     "reference chromosomes keeps it. An id that is in the "
                     "reference-chromosome GTF (gencode.vX.annotation.gtf.gz) is one of "
                     "these.")
    return "WARNING: " + " ".join(parts)


class UnreadableMetaInfo(Exception):
    """A meta_info.json that is there but cannot be read as Salmon's; the message names it."""


def read_meta_info(quant_path):
    """The index fields of the Salmon ``aux_info/meta_info.json`` beside a quant.sf, or
    None when there is none (kallisto, or an older Salmon).

    Raises :class:`UnreadableMetaInfo` when the file is there but is not a JSON object in
    UTF-8 whose hashes are strings: :func:`index_provenance` treats that donor as having
    none, and says so, rather than refusing the cohort or failing with a traceback.
    """
    path = os.path.join(os.path.dirname(quant_path), "aux_info", "meta_info.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as f:
            meta = json.load(f)
    except (ValueError, OSError) as e:          # JSONDecodeError and UnicodeDecodeError
        raise UnreadableMetaInfo("%s: %s" % (path, e)) from e
    if not isinstance(meta, dict):
        raise UnreadableMetaInfo("%s: a JSON %s, not an object" % (path, type(meta).__name__))
    for key in _HASH_KEYS:
        if meta.get(key) is not None and not isinstance(meta[key], str):
            raise UnreadableMetaInfo("%s: %s is %r, not a string" % (path, key, meta[key]))
    return {key: meta.get(key) for key in META_KEYS}


def _index_key(meta):
    """``(index_seq_hash, index_decoy_seq_hash)``, the second None when not recorded."""
    return meta["index_seq_hash"], meta.get("index_decoy_seq_hash")


def index_provenance(donor_quants):
    """``{donor: quant.sf path}`` -> the sidecar document: each donor's index fields, the
    distinct ``index_seq_hash`` and ``index_decoy_seq_hash`` values, whether there is more
    than one index, the donors with no meta_info.json, and those whose meta_info.json could
    not be read (``{donor: why}``; their entry in ``donors`` is None too).

    Two donors are on different indexes when their ``index_seq_hash`` differs, or when both
    record an ``index_decoy_seq_hash`` and it differs.  A donor whose Salmon recorded no
    decoy hash is compared on the transcripts alone.
    """
    donors, unreadable = {}, {}
    for d, q in sorted(donor_quants.items()):
        try:
            donors[d] = read_meta_info(q)
        except UnreadableMetaInfo as e:
            donors[d], unreadable[d] = None, str(e)
    known = [m for m in donors.values() if m and m.get("index_seq_hash")]
    hashes = sorted({m["index_seq_hash"] for m in known})
    decoys = sorted({m["index_decoy_seq_hash"] for m in known
                     if m.get("index_decoy_seq_hash") is not None})
    return {"format": INDEX_FORMAT, "donors": donors, "index_seq_hashes": hashes,
            "index_decoy_seq_hashes": decoys,
            "mixed": len(hashes) > 1 or len(decoys) > 1,
            "missing_meta_info": [d for d, m in donors.items()
                                  if m is None and d not in unreadable],
            "unreadable_meta_info": unreadable}


def mixed_index_message(prov):
    """Which donors went with which index, for the error and for the warning."""
    by = {}
    for d, m in prov["donors"].items():
        if m and m.get("index_seq_hash"):
            by.setdefault(_index_key(m), []).append(d)

    def _one(key, ds):
        seq, decoy = key
        first = prov["donors"][ds[0]]
        if decoy is None:
            return ("index_seq_hash %s (%s targets): %s"
                    % (seq, first.get("num_valid_targets"), ", ".join(ds)))
        return ("index_seq_hash %s, index_decoy_seq_hash %s (%s targets, %s decoys): %s"
                % (seq, decoy, first.get("num_valid_targets"), first.get("num_decoy_targets"),
                   ", ".join(ds)))
    return "; ".join(_one(key, ds)
                     for key, ds in sorted(by.items(), key=lambda kv: (kv[0][0], kv[0][1] or "")))
