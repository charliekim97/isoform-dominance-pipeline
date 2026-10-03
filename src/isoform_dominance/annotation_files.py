"""Local annotation files: a GENCODE (or Ensembl) GTF and a transcript FASTA, read offline.

Ensembl 116 (June 2026) is the last release of the legacy Ensembl platform, and the REST
API that :mod:`isoform_dominance.annotate` and :mod:`isoform_dominance.identifiability`
read is kept "for e116 for long term use" only; the new platform has no REST API.  An
index built from GENCODE 51 or later can only be matched by its own files.  This module
reads them, with no network, and gives the GTF's genes as records shaped like the REST
``lookup/...?expand=1`` answer -- exactly the fields the package reads -- so that
:func:`isoform_dominance.annotate.transcripts_of` and everything after it run unchanged
on either source.

Three rules make the two sources agree; each was measured on a 109-gene survey against
REST at releases 110, 114 and 116 (GENCODE 44, 48 and 50):

1. The protein length is ``(sum of CDS nt + (3 - frame of the 5'-most CDS) % 3) // 3``
   (:func:`protein_length`).  GENCODE's CDS features leave out the stop codon, and a CDS
   that starts mid-codon (``cds_start_NF``) is padded at its 5' end as Ensembl pads it.
2. The canonical transcript carries the tag ``Ensembl_canonical`` exactly; not
   ``Ensembl_canonical_extended`` (new in GENCODE 50) and not MANE.
3. Clusters are listed in content order, in both modes
   (:func:`isoform_dominance.annotate.cluster_by_terminal_exon`): REST and the GTF list a
   gene's transcripts in different orders.

Ids lose their version.  GENCODE 25-43 write the chrY copy of a pseudoautosomal gene
under the chrX ids plus ``_PAR_Y``; such records are kept apart (``_par_y``) and are
neither a candidate gene nor background.  A region is the seqname without ``chr``, and
``chrM`` is ``MT``, as Ensembl names them.

Only the standard library is used.
"""
import contextlib
import gzip
import os
import re
import zlib

from . import io

#: Bytes read from the GTF at a time.  Only a block that holds a needle is split into
#: lines; parsing every line of GENCODE 50's comprehensive GTF took over 40 s.
BLOCK = 1 << 24

_ATTR = re.compile(r'\s*([^\s;"]+)\s+(?:"([^"]*)"|([^;\s]+))\s*;')
#: GENCODE's ``##description``: ``... version 50 (Ensembl 116)``, ``version M37 (...)``.
_DESC = re.compile(r"version (M?\d+) \(Ensembl (\d+)\)")
_FEATURES = ("gene", "transcript", "exon", "CDS")
_BASIC = b'tag "basic";'
#: A file of fewer transcripts than this is not called basic however its lines are tagged:
#: an extract of a few genes from the comprehensive GTF can hold only basic transcripts
#: (every single-transcript gene's is), and GENCODE's basic GTF holds tens of thousands.
MIN_BASIC_TRANSCRIPTS = 100
PAR_Y = "_PAR_Y"


class AnnotationFileError(io.InputError):
    """A GTF or FASTA that cannot give what was asked of it."""


class BasicGTF(AnnotationFileError):
    """A GENCODE ``basic`` GTF, which leaves out transcripts the index holds."""


@contextlib.contextmanager
def _reading(path):
    """A truncated or corrupt gzip file is one :class:`AnnotationFileError`, not a trace:
    a half-downloaded GTF is the likeliest bad input."""
    try:
        yield
    except (EOFError, zlib.error, gzip.BadGzipFile) as e:
        raise AnnotationFileError("%s is not a complete gzip file (%s); download it again"
                                  % (path, e)) from e


def provenance(path):
    """``{file, bytes, sha256}`` of ``path``."""
    return {"file": os.path.basename(str(path)), "bytes": os.path.getsize(path),
            "sha256": io.file_sha256(path)}


def header(path):
    """What the GTF's header says: ``provider``, ``description``, ``date``,
    ``gencode_release`` and ``ensembl_release``, None for what it does not say.

    The releases come from GENCODE's ``##description`` (``... version 50 (Ensembl
    116)``); ``gencode_release`` is an integer, or a string such as ``"M37"`` for mouse.
    An Ensembl GTF (``#!genome-build ...``) names no release, nor may a GENCODE release
    after 50, whose header has not been seen; a header with none is accepted.
    """
    out = {}
    with _reading(path), io.open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##") and ":" in line:
                k, v = line[2:].split(":", 1)
                out[k.strip()] = v.strip()
            elif line.startswith("#!") and " " in line:
                k, v = line[2:].split(" ", 1)
                out[k.strip()] = v.strip()
    m = _DESC.search(out.get("description", ""))
    gencode = None
    if m:
        gencode = int(m.group(1)) if m.group(1).isdigit() else m.group(1)
    return {"provider": out.get("provider"), "description": out.get("description"),
            "date": out.get("date") or out.get("genebuild-last-updated"),
            "gencode_release": gencode,
            "ensembl_release": int(m.group(2)) if m else None}


def region(seqname):
    """GENCODE ``chr1``/``chrM`` -> Ensembl ``1``/``MT``."""
    c = seqname[3:] if seqname.startswith("chr") else seqname
    return "MT" if c == "M" else c


def _unv(ident):
    """``ENSG00000182378.15_PAR_Y`` -> ``ENSG00000182378_PAR_Y``; others lose the version."""
    base = ident.split(".")[0]
    if ident.endswith(PAR_Y) and not base.endswith(PAR_Y):
        base += PAR_Y
    return base


def _blocks(path, block=BLOCK):
    """The file's bytes in pieces of about ``block``, each ending at a line's end: the part
    of a block after its last newline is carried into the next."""
    with _reading(path), io.open_bytes(path) as fh:
        rest = b""
        while True:
            blk = fh.read(block)
            if not blk:
                break
            blk = rest + blk
            cut = blk.rfind(b"\n") + 1
            body, rest = blk[:cut], blk[cut:]
            if body:
                yield body
        if rest:
            yield rest + b"\n"


def _matching_lines(path, found, block, stats):
    """Decoded lines of ``path`` for which ``found(bytes)`` is true, and in ``stats`` the
    count of feature lines below the gene -- a gene line is never tagged ``basic`` -- and of
    those tagged ``basic``.  ``found`` is first asked of a whole block, and a block it does
    not match is not split into lines."""
    for body in _blocks(path, block):
        n = body.count(b"\n") - body.count(b"\n#") - (1 if body.startswith(b"#") else 0)
        stats["lines"] += n - body.count(b"\tgene\t")
        stats["transcripts"] += body.count(b"\ttranscript\t")
        stats["basic"] += body.count(_BASIC)
        if found(body):
            for line in body.split(b"\n"):
                if line and found(line):
                    yield line.decode()


def _attrs(field):
    a = {}
    for k, v1, v2 in _ATTR.findall(field):
        v = v1 if (v1 or not v2) else v2
        if k == "tag":
            a.setdefault("tag", []).append(v)
        else:
            a.setdefault(k, v)
    return a


def scan(path, symbol=None, gene_id=None, transcript_ids=(), block=BLOCK, info=None):
    """Every gene of the GTF named ``symbol``, with ``gene_id``, or holding one of
    ``transcript_ids``, as REST-shaped records (see the module docstring).

    ``symbol`` is matched exactly first, and only if no gene has that exact name is the
    file read again for a case-insensitive match (REST's ``lookup/symbol`` ignores case);
    ``info["case_insensitive"]`` says so.  A gene found through ``transcript_ids`` alone
    holds only the lines of those transcripts.  A GENCODE ``basic`` GTF -- every line
    below a gene tagged ``basic``, in a file of at least :data:`MIN_BASIC_TRANSCRIPTS`
    transcripts -- raises :class:`BasicGTF`: it leaves out transcripts the index holds, and
    at release 110 it changed 51 of 109 proposals.
    """
    info = {} if info is None else info
    needles, keeps = [], []
    if gene_id:
        base = gene_id.split(".")[0]
        needles += [b'gene_id "' + base.encode() + b'.', b'gene_id "' + base.encode() + b'"']
        keeps.append(lambda a: a["gene_id"].split(".")[0] == base)
    if symbol:
        needles.append(('gene_name "%s";' % symbol).encode())
        keeps.append(lambda a: a.get("gene_name") == symbol)
    tids = {t.split(".")[0] for t in transcript_ids}
    for t in sorted(tids):
        needles += [b'transcript_id "' + t.encode() + b'.', b'transcript_id "' + t.encode()
                    + b'"']
    if tids:
        keeps.append(lambda a: a.get("transcript_id", "").split(".")[0] in tids)
    if not needles:
        raise ValueError("scan needs a symbol, a gene id or transcript ids")

    def found(text):
        return any(n in text for n in needles)

    stats = {"lines": 0, "transcripts": 0, "basic": 0}
    genes = _records(_matching_lines(path, found, block, stats),
                     lambda a: any(k(a) for k in keeps))
    _refuse_basic(path, stats)
    info["case_insensitive"] = False
    if symbol and not gene_id and not any(g["display_name"] == symbol for g in genes):
        # bytes.lower() and a plain search: re.IGNORECASE over every block took three times
        # as long as the first read on a file of GENCODE 50's size
        low = ('gene_name "%s";' % symbol).lower().encode()
        loose = _records(_matching_lines(path, lambda text: low in text.lower(), block,
                                         {"lines": 0, "transcripts": 0, "basic": 0}),
                         lambda a: (a.get("gene_name") or "").upper() == symbol.upper())
        if loose:
            info["case_insensitive"] = True
            have = {g["id"] for g in genes}
            genes += [g for g in loose if g["id"] not in have]
    return genes


def _refuse_basic(path, stats):
    if stats["transcripts"] >= MIN_BASIC_TRANSCRIPTS and stats["basic"] >= stats["lines"]:
        raise BasicGTF(
            "%s is a GENCODE basic GTF: every one of its %d transcript, exon, CDS and UTR "
            "lines is tagged \"basic\". It leaves out transcripts a Salmon index of the full "
            "transcriptome holds (at release 110 it changed 51 of 109 proposals); pass the "
            "comprehensive gencode.vN.annotation.gtf.gz" % (path, stats["lines"]))


def _records(lines, keep):
    genes = {}
    for line in lines:
        if not line or line[0] == "#":
            continue
        p = line.rstrip("\r\n").split("\t")
        if len(p) < 9 or p[2] not in _FEATURES:
            continue
        a = _attrs(p[8])
        if "gene_id" not in a or not keep(a):
            continue                    # the needle matched some other attribute
        gfull = a["gene_id"]
        if "." not in gfull and "gene_version" in a:
            gfull = "%s.%s" % (gfull, a["gene_version"])
        gid = _unv(gfull)
        g = genes.get(gid)
        if g is None:
            g = genes[gid] = {"id": gid.replace(PAR_Y, ""), "version": None,
                              "display_name": a.get("gene_name"),
                              "seq_region_name": region(p[0]),
                              "strand": 1 if p[6] == "+" else -1,
                              "start": None, "end": None, "biotype": None,
                              "canonical_transcript": None, "Transcript": [],
                              "_tx": {}, "_par_y": gid.endswith(PAR_Y)}
        if p[2] == "gene":
            ver = gfull.split(".")[1].split("_")[0] if "." in gfull else None
            g.update(display_name=a.get("gene_name"), start=int(p[3]), end=int(p[4]),
                     biotype=a.get("gene_type") or a.get("gene_biotype"),
                     version=int(ver) if ver and ver.isdigit() else None)
            continue
        if "transcript_id" not in a:
            continue
        tfull = a["transcript_id"]
        if "." not in tfull and "transcript_version" in a:
            tfull = "%s.%s" % (tfull, a["transcript_version"])
        tid = _unv(tfull)
        t = g["_tx"].get(tid)
        if t is None:
            t = g["_tx"][tid] = {"id": tid.replace(PAR_Y, ""), "version": None,
                                 "biotype": None, "is_canonical": 0, "Translation": None,
                                 "Exon": [], "_cds": [], "_versioned": tfull}
            g["Transcript"].append(t)
        if p[2] == "transcript":
            ver = tfull.split(".")[1].split("_")[0] if "." in tfull else None
            t.update(version=int(ver) if ver and ver.isdigit() else None,
                     biotype=a.get("transcript_type") or a.get("transcript_biotype"),
                     is_canonical=int("Ensembl_canonical" in a.get("tag", ())))
            if t["is_canonical"]:
                g["canonical_transcript"] = tfull.replace(PAR_Y, "")
        elif p[2] == "exon":
            t["Exon"].append({"start": int(p[3]), "end": int(p[4])})
        else:
            t["_cds"].append((int(p[3]), int(p[4]), p[7]))
    for g in genes.values():
        for t in g["Transcript"]:
            length = protein_length(t.pop("_cds"), g["strand"])
            t["Translation"] = {"length": length} if length else None
        del g["_tx"]
    return list(genes.values())


def protein_length(cds, strand):
    """Ensembl's ``Translation.length`` from a transcript's GTF CDS features
    ``[(start, end, frame)]``, which leave out the stop codon.

    Ensembl pads the 5' end of a CDS that starts mid-codon (``cds_start_NF``) with
    ``(3 - frame) % 3`` N's, the GTF frame being that of the 5'-most CDS -- on the minus
    strand the one that ends furthest along -- and drops a trailing partial codon:
    ``(nt + pad) // 3``.  This equals REST's ``Translation.length`` for 1292/1292 (release
    110 against GENCODE 44), 1350/1350 (114, 48) and 3680/3680 (116, 50) translated
    transcripts of a 109-gene survey; ``nt // 3`` is wrong for 80, 77 and 77 of them, and
    for 7,412 CDS-bearing transcripts of GENCODE 50.
    """
    if not cds:
        return None
    nt = sum(e - s + 1 for s, e, _ in cds)
    first = min(cds, key=lambda c: c[0]) if strand == 1 else max(cds, key=lambda c: c[1])
    frame = int(first[2]) if first[2] in ("0", "1", "2") else 0
    return (nt + (3 - frame) % 3) // 3


def fasta_id(head):
    """The versioned id of a FASTA header (without ``>``): the first ``|``-field of a
    GENCODE header, the first word of an Ensembl one; None when it has none."""
    words = head.strip().split("|")[0].split()
    return words[0] if words else None


def read_fasta(path, ids):
    """``{unversioned id: (versioned id, sequence)}`` for the records of ``ids``.

    ``_PAR_Y`` records are skipped.  Two records of one of ``ids`` raise
    :class:`AnnotationFileError`: which of them the index holds cannot be told.
    """
    want = {i.split(".")[0] for i in ids}
    out, cur, first, chunks = {}, None, None, []

    def flush():
        if cur is not None:
            if cur in out:
                raise AnnotationFileError("%s holds %s twice (%s and %s)"
                                          % (path, cur, out[cur][0], first))
            out[cur] = (first, "".join(chunks).upper())

    with _reading(path), io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                first = fasta_id(line[1:])
                cur = first.split(".")[0] if first else None
                if cur not in want or first.endswith(PAR_Y):
                    cur = None
                chunks = []
            elif cur is not None:
                chunks.append(line.strip())
        flush()
    return out


def sequences_for(gene, fasta, versions=None):
    """Every transcript of ``gene`` (a :func:`scan` record) from ``fasta``: ``{id: seq}``.

    A transcript the FASTA lacks, or holds at another version, is an
    :class:`AnnotationFileError`: a background silently short of transcripts, or sequence
    of another release, changes the answer.  ``versions``, a dict, receives each id's
    versioned id as the FASTA gives it.
    """
    want = {t["id"]: t for t in gene["Transcript"]}
    got = read_fasta(fasta, want)
    missing = sorted(set(want) - set(got))
    if missing:
        raise AnnotationFileError(
            "%s has no record for %d of the %d transcripts the GTF gives %s (%s%s). A subset "
            "FASTA -- GENCODE's basic or pc_transcripts, or Ensembl's cdna.all, which has no "
            "ncRNA -- does not hold them all; pass gencode.vN.transcripts.fa.gz of the GTF's "
            "release" % (fasta, len(missing), len(want), gene["id"], ", ".join(missing[:5]),
                         " and %d more" % (len(missing) - 5) if len(missing) > 5 else ""))
    wrong = sorted(t for t, (v, _) in got.items() if version(v) is not None
                   and want[t]["version"] is not None and version(v) != want[t]["version"])
    if wrong:
        t = wrong[0]
        raise AnnotationFileError(
            "%s and the GTF disagree on the version of %d transcript(s) of %s (%s is %s in "
            "the FASTA and %s.%s in the GTF): they are of different releases"
            % (fasta, len(wrong), gene["id"], t, got[t][0], t, want[t]["version"]))
    if versions is not None:
        versions.update((t, v) for t, (v, _) in got.items())
    return {t: s for t, (_, s) in got.items()}


def version(versioned):
    """``ENST00000349533.11`` -> 11, ``ENST00000381192.10_PAR_Y`` -> 10; None without a
    version."""
    if "." not in versioned:
        return None
    v = versioned.split(".")[1].split("_")[0]
    return int(v) if v.isdigit() else None
