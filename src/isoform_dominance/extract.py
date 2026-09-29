"""Extract per-donor isoform-group TPM from Salmon quant.sf (stdlib only)."""
import csv, os, glob, json
from . import index_scope
from .io import transcript_to_group, load_sample_map


def quant_paths(quantdir):
    """``{donor: quantdir/<donor>/quant.sf}``; FileNotFoundError if there is none."""
    quants = sorted(glob.glob(os.path.join(quantdir, "*", "quant.sf")))
    if not quants:
        raise FileNotFoundError(
            "No Salmon output found at %s. Expected one quant.sf per donor at "
            "%s/<donor>/quant.sf." % (os.path.join(quantdir, "*", "quant.sf"), quantdir))
    return {os.path.basename(os.path.dirname(q)): q for q in quants}


def extract(config, quantdir, samplemap, cohort):
    """Return [(donor, condition, {group: tpm}), ...]. Reads quantdir/<donor>/quant.sf."""
    groups = config["groups"]
    tx2grp = transcript_to_group(groups)
    cond = load_sample_map(samplemap)
    rows = []
    for donor, q in sorted(quant_paths(quantdir).items()):
        gt = dict.fromkeys(groups, 0.0)
        with open(q) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or "Name" not in reader.fieldnames or "TPM" not in reader.fieldnames:
                raise ValueError(
                    "%s is not a valid Salmon quant.sf (missing 'Name'/'TPM' columns); "
                    "found columns: %s" % (q, reader.fieldnames))
            for row in reader:
                g = tx2grp.get(row["Name"].split(".")[0])
                if g:
                    gt[g] += float(row["TPM"])
        rows.append((donor, cond.get(donor, "NA"), gt))
    return rows


def write_perdonor(config, rows, cohort, out):
    groups = config["groups"]
    pc = config.get("primary_comparison", list(groups)[:2])
    has_pair = len(pc) == 2
    a, b = (pc[0], pc[1]) if has_pair else (None, None)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        head = ["cohort", "donor", "condition"] + ["%s_TPM" % g for g in groups]
        if has_pair:
            head.append("%s_fraction" % a)
        w.writerow(head)
        for donor, c, gt in rows:
            line = [cohort, donor, c] + ["%.4f" % gt[g] for g in groups]
            if has_pair:
                s = gt[a] + gt[b]
                line.append("%.4f" % (gt[a] / s) if s > 0 else "NA")
            w.writerow(line)
    return out


def _quant_names(path):
    """A function yielding the Name column of one quant.sf, for :func:`index_scope.copies_in`."""
    def names():
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                yield row["Name"]
    return names


def _first_row_name(path):
    with open(path) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"), None)
    return row["Name"] if row else ""


def quant_copies(config, prov, quants):
    """Same-name copies among the targets of each index the cohort was quantified against.

    Only an index built without Salmon's ``--gencode`` keeps GENCODE's whole header as the
    target name, and only then does quant.sf say which gene each target belongs to; plain
    names give nothing to check, and nothing is reported.
    """
    one_per_index = {}
    for donor, q in sorted(quants.items()):
        meta = prov["donors"].get(donor) or {}
        one_per_index.setdefault(meta.get("index_seq_hash"), q)
    targets = list(transcript_to_group(config["groups"]))
    found = {}
    for q in one_per_index.values():
        if "|" not in _first_row_name(q):
            continue
        for c in index_scope.copies_in(_quant_names(q), targets,
                                       [config["gene"]] if config.get("gene") else ()):
            found.setdefault(c["gene_id"], c)
    return [found[g] for g in sorted(found)]


def run(config, quantdir, samplemap, cohort, out, allow_mixed_index=False, notes=None):
    """Write the per-donor CSV ``out`` and, beside it, ``<out>.index.json``: which Salmon
    index quantified each donor, from its ``aux_info/meta_info.json``.

    Donors quantified against different indexes (different ``index_seq_hash``) are
    refused with :class:`index_scope.MixedIndexError` before anything is written, unless
    ``allow_mixed_index``.  Warnings -- a combined mixed cohort, donors with no
    meta_info.json, same-name copies among the index's targets -- are appended to
    ``notes`` when a list is given.
    """
    notes = [] if notes is None else notes
    quants = quant_paths(quantdir)
    prov = index_scope.index_provenance(quants)
    if prov["mixed"]:
        which = index_scope.mixed_index_message(prov)
        if not allow_mixed_index:
            raise index_scope.MixedIndexError(
                "the donors of cohort %s were quantified against different Salmon indexes "
                "(%s). Their TPMs are not comparable: quantify them against one index, or "
                "pass --allow-mixed-index to combine them anyway" % (cohort, which))
        notes.append("WARNING: the donors of cohort %s were quantified against different "
                     "Salmon indexes (%s); combined because of --allow-mixed-index"
                     % (cohort, which))
    if prov["missing_meta_info"]:
        notes.append("WARNING: no aux_info/meta_info.json for donor(s) %s, so which index "
                     "quantified them is not recorded, and a mix of indexes cannot be "
                     "detected (kallisto, or Salmon output without aux_info)"
                     % ", ".join(prov["missing_meta_info"]))
    rows = extract(config, quantdir, samplemap, cohort)
    copies = quant_copies(config, prov, quants)
    warning = index_scope.copy_warning(copies, "the index behind %s" % quantdir)
    if warning:
        notes.append(warning)
    write_perdonor(config, rows, cohort, out)
    with open(out + ".index.json", "w") as f:
        json.dump(dict(prov, cohort=cohort, same_name_copies=copies), f, indent=1)
    return len(rows)
