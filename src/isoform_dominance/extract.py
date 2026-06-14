"""Extract per-donor isoform-group TPM from Salmon quant.sf (stdlib only)."""
import csv, os, glob
from .io import transcript_to_group, load_sample_map, primary_pair


def extract(config, quantdir, samplemap, cohort):
    """Return [(donor, condition, {group: tpm}), ...]. Reads quantdir/<donor>/quant.sf."""
    groups = config["groups"]
    tx2grp = transcript_to_group(groups)
    cond = load_sample_map(samplemap)
    rows = []
    for q in sorted(glob.glob(os.path.join(quantdir, "*", "quant.sf"))):
        donor = os.path.basename(os.path.dirname(q))
        gt = {g: 0.0 for g in groups}
        with open(q) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = tx2grp.get(row["Name"].split(".")[0])
                if g:
                    gt[g] += float(row["TPM"])
        rows.append((donor, cond.get(donor, "NA"), gt))
    return rows


def write_perdonor(config, rows, cohort, out):
    groups = config["groups"]
    a, b = primary_pair(config)
    has_pair = len(config.get("primary_comparison", list(groups)[:2])) == 2
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
