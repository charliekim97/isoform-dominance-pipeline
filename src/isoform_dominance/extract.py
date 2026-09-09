"""Extract per-donor isoform-group TPM from Salmon quant.sf (stdlib only)."""
import csv, os, glob
from .io import transcript_to_group, load_sample_map


def extract(config, quantdir, samplemap, cohort):
    """Return [(donor, condition, {group: tpm}), ...]. Reads quantdir/<donor>/quant.sf."""
    groups = config["groups"]
    tx2grp = transcript_to_group(groups)
    cond = load_sample_map(samplemap)
    quants = sorted(glob.glob(os.path.join(quantdir, "*", "quant.sf")))
    if not quants:
        raise FileNotFoundError(
            "No Salmon output found at %s. Expected one quant.sf per donor at "
            "%s/<donor>/quant.sf." % (os.path.join(quantdir, "*", "quant.sf"), quantdir))
    rows = []
    for q in quants:
        donor = os.path.basename(os.path.dirname(q))
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


def run(config, quantdir, samplemap, cohort, out):
    rows = extract(config, quantdir, samplemap, cohort)
    write_perdonor(config, rows, cohort, out)
    return len(rows)
