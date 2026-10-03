"""Gene symbol -> proposed isoform groups via the Ensembl REST API, or a local GTF.

The hard part of isoform analysis is deciding which transcripts form a functional
group. This module clusters a gene's protein-coding transcripts by their 3' terminal-exon
splice acceptor — the alternative last exon that distinguishes functional isoform classes
(e.g. a long signalling form vs a short truncated form) — and proposes a two-group
comparison (canonical-isoform cluster vs the largest alternative cluster) that the user
reviews and renames before use.

A symbol can name more than one gene on the reference chromosomes, and Ensembl's
``lookup/symbol`` returns one of them without saying so: for the pseudoautosomal genes
CD99, CRLF2, CSF2RA, IL3RA and SHOX the chrY copy, for HERC3 and DUSP13B the newer of two
genes of that name.  :func:`build_config` looks for the others and chooses by
:data:`GENE_RULE`, or stops.
"""
import json
import urllib.parse

from . import annotation_files, ensembl
from .ensembl import DEFAULT_RETRIES, DEFAULT_RETRY_WAIT
from .index_scope import REFERENCE_REGIONS

ENSEMBL = ensembl.SERVER

#: How :func:`build_config` chooses among genes that share the symbol; recorded in a
#: config's ``_gene_choice`` whenever there was a choice to make.
GENE_RULE = ("the gene on the reference chromosomes (1-22, X, Y, MT) whose name is the "
             "symbol; of a chrX/chrY pair -- a pseudoautosomal gene -- the chrX one, the copy "
             "a Salmon index built from GENCODE keeps; any other choice is --gene-id's")


class AmbiguousGene(ValueError):
    """A symbol names more than one gene on the reference chromosomes, and none is chosen."""


class NotOnReference(ValueError):
    """A human symbol whose every gene lies off the reference chromosomes (an alternate
    locus, a patch or a scaffold), so the recommended index does not contain it."""


#: The species :data:`GENE_RULE`'s reference chromosomes and chrX/chrY pairs are about.
REFERENCE_SPECIES = "homo_sapiens"


def _off_reference(g, species):
    """Is ``g`` known to lie off the reference chromosomes?  Only for
    :data:`REFERENCE_SPECIES`, and only when the record names its region."""
    region = g.get("seq_region_name")
    return species == REFERENCE_SPECIES and region is not None \
        and region not in REFERENCE_REGIONS


def _not_on_reference(gene, records):
    lines = "; ".join("%s on %s" % (g["id"], g.get("seq_region_name")) for g in records)
    return NotOnReference(
        "%s has no gene on the reference chromosomes (1-22, X, Y, MT), only %s. The "
        "recommended reference-chromosome index does not contain it; to propose groups "
        "from one of these anyway, pass it with --gene-id" % (gene, lines))


def _q(text):
    """``text`` as one URL path segment: a symbol can hold a space or a slash."""
    return urllib.parse.quote(str(text), safe="")


def _get(path, timeout=ensembl.DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES,
         retry_wait=DEFAULT_RETRY_WAIT, server=None):
    return ensembl.get_json(path, server=server, timeout=timeout, retries=retries,
                            retry_wait=retry_wait)


def _post(path, body, timeout=ensembl.DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES,
          retry_wait=DEFAULT_RETRY_WAIT, server=None):
    return ensembl.request_json(path, body, server=server, timeout=timeout, retries=retries,
                                retry_wait=retry_wait)


def fetch_transcripts(gene, species="homo_sapiens", retries=DEFAULT_RETRIES,
                      retry_wait=DEFAULT_RETRY_WAIT, server=None):
    """Return {gene, gene_id, species, strand, transcripts:[{id, protein_aa, terminal_acceptor, is_canonical}]}.

    For the gene ``lookup/symbol`` gives; :func:`build_config` also checks for other genes
    of the name.  The lookup is retried as described in :mod:`isoform_dominance.ensembl`.
    ``server`` is a base URL from :func:`isoform_dominance.ensembl.resolve_server`; None
    means the current release.
    """
    g = _get("/lookup/symbol/%s/%s?expand=1" % (_q(species), _q(gene)),
             retries=retries, retry_wait=retry_wait, server=server)
    return transcripts_of(g, gene, species)


def transcripts_of(g, gene, species="homo_sapiens"):
    """:func:`fetch_transcripts`' result for an expanded ``lookup`` record ``g``."""
    strand = g["strand"]
    canonical = (g.get("canonical_transcript") or "").split(".")[0]
    out = []
    for t in g.get("Transcript", []):
        if t.get("biotype") != "protein_coding":
            continue
        tl = t.get("Translation") or {}
        plen = tl.get("length")
        if not plen:
            continue
        exons = t["Exon"]
        if strand == 1:
            term = max(exons, key=lambda e: e["end"]); acc = term["start"]
        else:
            term = min(exons, key=lambda e: e["start"]); acc = term["end"]
        tid = t["id"].split(".")[0]
        out.append({"id": tid, "protein_aa": plen, "terminal_acceptor": acc,
                    "is_canonical": (t.get("is_canonical", 0) == 1) or (tid == canonical)})
    if not out:
        raise ValueError("No protein-coding transcripts with a translation found for %s" % gene)
    return {"gene": gene, "gene_id": g.get("id"), "species": species, "strand": strand,
            "transcripts": out}


def _where(g):
    """A candidate as ``_gene_choice`` records it."""
    return {"gene_id": g["id"],
            "location": "%s:%s-%s" % (g.get("seq_region_name"), g.get("start"), g.get("end")),
            "n_transcripts": len(g.get("Transcript") or [])}


def choose_gene(gene, species, looked_up, **net):
    """The gene ``gene`` means, and the record of a choice when there was one to make.

    ``looked_up`` is the expanded record ``lookup/symbol`` gave.  The other genes of the
    name come from ``xrefs/symbol``, which also lists genes that carry the symbol only as
    a synonym (SMN2 for SMN1) and copies on alternate loci; one ``lookup/id`` for those
    keeps the genes whose display name is the symbol and, for human, that lie on a
    reference chromosome (:data:`isoform_dominance.index_scope.REFERENCE_REGIONS`).  With
    one such gene nothing is recorded, and when xrefs/symbol lists no gene besides
    ``looked_up`` no ``lookup/id`` is made.  Of a human chrX/chrY pair the chrX gene is
    taken (:data:`GENE_RULE`).  Any other set of two or more raises :class:`AmbiguousGene`
    listing them.  A human symbol none of whose genes is on a reference chromosome raises
    :class:`NotOnReference`; ``--gene-id`` takes one of them anyway.

    The region rule and the chrX/chrY rule are about GRCh38's chromosome names, so for any
    other species every gene of the name is a candidate: zebrafish chromosomes run to 25,
    and fly and worm ones have other names altogether.

    Returns ``(expanded gene record, choice or None)``; the choice is
    :func:`choose_among`'s.
    """
    xr = _get("/xrefs/symbol/%s/%s?object_type=gene" % (_q(species), _q(gene)), **net)
    if not isinstance(xr, list):
        raise ValueError("Ensembl xrefs/symbol for %s answered a %s, not a list"
                         % (gene, type(xr).__name__))
    ids = sorted({x["id"] for x in xr if isinstance(x, dict) and x.get("type") == "gene"
                  and str(x.get("id", "")).startswith("ENS")})
    if not ids:
        if _off_reference(looked_up, species):
            raise _not_on_reference(gene, [looked_up])
        return looked_up, {"rule": GENE_RULE, "chosen": looked_up.get("id"),
                           "candidates": [],
                           "reason": "xrefs/symbol listed no gene for %s, so other genes of "
                                     "that name were not looked for" % gene}
    others = [i for i in ids if i != looked_up.get("id")]
    if not others:
        if _off_reference(looked_up, species):
            raise _not_on_reference(gene, [looked_up])
        return looked_up, None
    got = _post("/lookup/id", {"ids": others, "expand": 1}, **net)
    name = (looked_up.get("display_name")
            if (looked_up.get("display_name") or "").upper() == gene.upper() else gene)
    records = [looked_up] + [got[i] for i in others if isinstance(got, dict) and got.get(i)]
    g, choice = choose_among(gene, [g for g in records if g.get("display_name") == name],
                             species)
    return (g or looked_up), choice


def choose_among(gene, named, species=REFERENCE_SPECIES):
    """The gene ``gene`` means among ``named``, the expanded records of every gene whose
    display name is the symbol, and the record of the choice when there was one to make.

    The rule both sources follow, REST (:func:`choose_gene`) and a GTF
    (:func:`build_config_from_gtf`).  For human, only a gene on a reference chromosome is a
    candidate, and when every one of ``named`` is known to lie off them
    :class:`NotOnReference` is raised.  With one candidate nothing is recorded; of a human
    chrX/chrY pair the chrX gene is taken (:data:`GENE_RULE`); any other set of two or more
    raises :class:`AmbiguousGene` listing them.  Returns ``(record or None, choice or
    None)``, None when no candidate is left.
    """
    human = species == REFERENCE_SPECIES
    named = sorted(named, key=lambda g: g["id"])
    genes = [g for g in named if g.get("seq_region_name") in REFERENCE_REGIONS] if human \
        else named
    if not genes and named and all(_off_reference(g, species) for g in named):
        raise _not_on_reference(gene, named)
    if len(genes) <= 1:
        return (genes[0] if genes else None), None
    choice = {"rule": GENE_RULE, "candidates": [_where(g) for g in genes]}
    by_region = {g["seq_region_name"]: g for g in genes}
    if human and len(genes) == 2 and set(by_region) == {"X", "Y"}:
        x, y = by_region["X"], by_region["Y"]
        choice.update(chosen=x["id"], reason=(
            "%s is a pseudoautosomal gene: %s on chrX and %s on chrY, one sequence. Salmon "
            "keeps only the first of identical sequences -- the chrX copy in a GENCODE FASTA, "
            "the chrY copy in Ensembl's cDNA FASTA, which lists chrY first -- so the chrX "
            "gene is used; for an index built from Ensembl cDNA, use --gene-id %s"
            % (gene, x["id"], y["id"], y["id"])))
        return x, choice
    raise AmbiguousGene(
        "%s names %d genes%s: %s. Pick one with --gene-id"
        % (gene, len(genes), " on the reference chromosomes" if human else "", "; ".join(
            "%s at %s, %d transcript%s" % (c["gene_id"], c["location"], c["n_transcripts"],
                                           "" if c["n_transcripts"] == 1 else "s")
            for c in choice["candidates"])))


def cluster_by_terminal_exon(info):
    """Group transcripts by 3' terminal-exon acceptor coordinate.

    Clusters are listed by content -- most transcripts first, then the longer
    representative protein, then the lower acceptor coordinate -- and so are a config's
    ``_clusters`` and ``_proposal.tied_with``.  Through 2.5 clusters with as many
    transcripts kept the order in which the source listed the transcripts, and REST at
    release 116 and the GENCODE 50 GTF list them in different orders for 104 of 109 survey
    genes."""
    clusters = {}
    for t in info["transcripts"]:
        clusters.setdefault(t["terminal_acceptor"], []).append(t)
    out = []
    for acc, txs in clusters.items():
        lens = sorted(x["protein_aa"] for x in txs)
        out.append({"acceptor": acc, "rep_aa": lens[len(lens) // 2], "n": len(txs),
                    "canonical": any(x["is_canonical"] for x in txs),
                    "ids": sorted(x["id"] for x in txs)})
    return sorted(out, key=lambda c: (-c["n"], -c["rep_aa"], c["acceptor"]))


#: How :func:`propose_groups` picks the alternative class; recorded in every config.
ALTERNATIVE_RULE = ("the non-canonical cluster with the most transcripts; a tie goes to the "
                    "longer representative protein, and a tie on both to the lower terminal-"
                    "acceptor coordinate")


def alternative_ties(clusters, groups, primary):
    """The clusters the alternative class was chosen over on protein length alone.

    :func:`propose_groups` takes the non-canonical cluster with the most transcripts.  When
    another cluster has as many, the choice fell to representative protein length -- a
    tie-break, not a biological criterion -- and a later release that adds one transcript
    to either cluster changes the proposed comparison.  Across six Ensembl releases (110,
    112-116), 26 to 39 of the 84 to 100 two-class proposals in a 109-gene survey were ties
    of this kind.  Of the 84 genes proposed a pair at every one of them, 30 were proposed a
    different pair at some release; 12 of the 30 had been a tie at the release before the
    first change, and in 8 -- LEPR and NTRK3 among them -- the new alternative was a
    cluster the old one had been tied with.  In all 30 it was the alternative class that
    changed, never the canonical one.
    """
    if len(primary) < 2:
        return []
    alt_ids, canon_ids = set(groups[primary[0]]), set(groups[primary[1]])
    alt = next(c for c in clusters if set(c["ids"]) == alt_ids)
    return [c for c in clusters
            if set(c["ids"]) not in (alt_ids, canon_ids) and c["n"] == alt["n"]]


def propose_groups(info):
    """The canonical cluster, the alternative (:data:`ALTERNATIVE_RULE`) and all clusters.

    Every tie is broken by content, down to the acceptor coordinate, so the proposal is a
    function of the annotation alone.  Before 2.4 a tie on both transcript count and
    protein length fell to the order in which the server listed the transcripts, which
    nothing guarantees: 249 of the 14,054 two-class genes of GENCODE 50 (1.8%) are such
    ties, and for FOXO1 and STK11 the REST order and the GTF order pick differently.
    """
    clusters = cluster_by_terminal_exon(info)
    canon = (next((c for c in clusters if c["canonical"]), None)
             or max(clusters, key=lambda c: (c["rep_aa"], -c["acceptor"])))
    others = sorted([c for c in clusters if c is not canon],
                    key=lambda c: (-c["n"], -c["rep_aa"], c["acceptor"]))
    alt = others[0] if others else None
    lbl_canon = "iso_%daa" % canon["rep_aa"]
    groups = {lbl_canon: canon["ids"]}
    primary = [lbl_canon]
    if alt:
        lbl_alt = "iso_%daa" % alt["rep_aa"]
        if lbl_alt == lbl_canon:
            lbl_alt += "_alt"
        groups[lbl_alt] = alt["ids"]
        primary = [lbl_alt, lbl_canon]  # alternative (often shorter) first
    return groups, primary, clusters


def _given(g, gene, gene_id, species):
    """The ``_gene_choice`` of a gene named by ``--gene-id``, which must be a gene of the
    name ``gene``."""
    if (g.get("display_name") or "").upper() != gene.upper():
        raise ValueError("%s is %s, not %s" % (gene_id, g.get("display_name"), gene))
    choice = {"rule": GENE_RULE, "chosen": g["id"], "candidates": [_where(g)],
              "reason": "given by --gene-id"}
    if _off_reference(g, species):
        choice["reason"] += (
            "; %s is on %s, not a reference chromosome: the recommended "
            "reference-chromosome index does not contain this gene"
            % (g["id"], g.get("seq_region_name")))
    return choice


def build_config(gene, species="homo_sapiens", release=None, gene_id=None, **retry):
    """Build a reviewable config.json dict for `gene` from Ensembl annotation.

    ``release`` is the Ensembl release to propose the groups from; None means the one
    ``rest.ensembl.org`` currently serves, and an earlier one is read from Ensembl's
    REST archive (see :func:`isoform_dominance.ensembl.resolve_server`).
    ``ensembl_release`` records the release the groups were proposed from, as the server
    reported it.  The identifiability verdict is a function of that release, so a
    config without it cannot be re-run to the same answer.

    ``gene_id`` names the Ensembl gene outright, for a symbol that names several (see
    :func:`choose_gene`); it must be a gene of that name.  ``gene_id`` in the config is the
    gene the groups were proposed from, and ``_gene_choice`` records how it was chosen
    when there was a choice.
    """
    server = ensembl.resolve_server(release, **retry)
    net = dict(retry, server=server)
    if gene_id:
        g = _get("/lookup/id/%s?expand=1" % _q(gene_id.split(".")[0]), **net)
        choice = _given(g, gene, gene_id, species)
    else:
        g = _get("/lookup/symbol/%s/%s?expand=1" % (_q(species), _q(gene)), **net)
        g, choice = choose_gene(gene, species, g, **net)
    info = transcripts_of(g, gene, species)
    release = ensembl.release_number(_get("/info/data", **net))
    return _config(gene, species, release, info, choice)


def _config(gene, species, release, info, choice):
    """The config both sources write, from :func:`transcripts_of`'s ``info``."""
    groups, primary, clusters = propose_groups(info)
    ties = alternative_ties(clusters, groups, primary)
    cfg = {
        "gene": gene, "gene_id": info["gene_id"], "species": species,
        "ensembl_release": release,
        "groups": groups, "primary_comparison": primary,
        "_proposal": {"alternative_rule": ALTERNATIVE_RULE,
                      "tied_with": [{"terminal_acceptor": c["acceptor"],
                                     "rep_protein_aa": c["rep_aa"],
                                     "n_transcripts": c["n"],
                                     "transcripts": c["ids"]} for c in ties]},
        "_proposed": ("Auto-proposed by `isoform-dominance annotate`. Groups = protein-coding "
                      "transcripts sharing a 3' terminal-exon splice acceptor (isoform-defining "
                      "alternative last exon). REVIEW and rename to functional names "
                      "(e.g. short/long) before use; smaller clusters are listed under _clusters."),
        "_clusters": [{"terminal_acceptor": c["acceptor"], "rep_protein_aa": c["rep_aa"],
                       "n_transcripts": c["n"], "contains_canonical": c["canonical"],
                       "transcripts": c["ids"]} for c in clusters],
    }
    if choice:
        cfg["_gene_choice"] = choice
    return cfg


def build_config_from_gtf(gene, gtf, species="homo_sapiens", gene_id=None, release=None,
                          notes=None, block=annotation_files.BLOCK):
    """:func:`build_config` from a local GTF, with no network.

    ``gtf`` is GENCODE's comprehensive ``gencode.vN.annotation.gtf.gz`` (or an Ensembl
    GTF); :mod:`isoform_dominance.annotation_files` reads it as records shaped like REST's,
    and the gene is chosen by the rule REST follows (:func:`choose_among`), so the config
    is the one :func:`build_config` writes from the release the file is of, but for
    ``annotation_source``: ``{kind: "gtf", file, bytes, sha256, provider,
    gencode_release, ensembl_release, date, description, n_transcripts}``, the file the
    groups were proposed from, what its header says, and how many transcripts of every
    biotype the GTF gives the gene -- the gene background ``identifiability --gtf`` uses.

    ``ensembl_release`` is the one the header names (``##description: ... version 50
    (Ensembl 116)``); ``release``, when given, must be the same.  A header that names none
    -- an Ensembl GTF, or a GENCODE release whose header has not been seen -- is accepted,
    and ``release`` is recorded as the file's.  What a REST run cannot see is appended to
    ``notes``: a header with no release, a symbol found only ignoring case, and a gene none
    of whose transcripts the GTF tags ``Ensembl_canonical`` (REST always names one).
    """
    af = annotation_files
    notes = [] if notes is None else notes
    hdr = af.header(gtf)
    if hdr["ensembl_release"] is None:
        notes.append("%s names no Ensembl release in its header; %s" % (
            gtf, "the release is recorded as %d, from --ensembl-release" % release
            if release is not None else "no release is recorded. Pass --ensembl-release N "
            "to record the one the file is of"))
    elif release is not None and int(release) != hdr["ensembl_release"]:
        raise ValueError("%s is Ensembl release %d (\"%s\"); --ensembl-release %d contradicts "
                         "it" % (gtf, hdr["ensembl_release"], hdr["description"], release))
    if gene_id:
        records = [g for g in af.scan(gtf, gene_id=gene_id, block=block) if not g["_par_y"]]
        if not records:
            raise ValueError("%s has no gene %s" % (gtf, gene_id.split(".")[0]))
        g = records[0]
        choice = _given(g, gene, gene_id, species)
    else:
        found = {}
        records = [g for g in af.scan(gtf, symbol=gene, info=found, block=block)
                   if not g["_par_y"]]
        named = [g for g in records if g["display_name"] == gene] \
            or [g for g in records if (g["display_name"] or "").upper() == gene.upper()]
        if not named:
            raise ValueError("%s has no gene named %s (gene_name, compared ignoring case too); "
                             "use the approved symbol, or --gene-id" % (gtf, gene))
        if found.get("case_insensitive"):
            notes.append("%s has no gene named exactly %s; %s was found ignoring case, as "
                         "Ensembl's REST lookup finds it" % (
                             gtf, gene, ", ".join(sorted({g["display_name"] for g in named}))))
        # every GTF record names its region, so a gene is chosen or NotOnReference raised
        g, choice = choose_among(gene, named, species)
    if not any(t["is_canonical"] for t in g["Transcript"]):
        notes.append("%s tags none of %s's transcripts Ensembl_canonical, so the canonical "
                     "class is the cluster with the longest representative protein; REST, "
                     "which names a canonical transcript, can propose another. Review the "
                     "groups." % (gtf, g.get("display_name") or gene))
    info = transcripts_of(g, gene, species)
    cfg = _config(gene, species, hdr["ensembl_release"] if hdr["ensembl_release"] is not None
                  else release, info, choice)
    cfg["annotation_source"] = dict(
        kind="gtf", **af.provenance(gtf),
        **{k: hdr[k] for k in ("provider", "gencode_release", "ensembl_release", "date",
                               "description")},
        n_transcripts=len(g["Transcript"]))
    return cfg


def run(gene, out, species="homo_sapiens", release=None, gene_id=None, gtf=None, notes=None,
        **retry):
    if gtf:
        cfg = build_config_from_gtf(gene, gtf, species, gene_id=gene_id, release=release,
                                    notes=notes)
    else:
        cfg = build_config(gene, species, release=release, gene_id=gene_id, **retry)
    with open(out, "w") as f:
        json.dump(cfg, f, indent=2)
    return cfg
