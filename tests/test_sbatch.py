"""`scripts/01_salmon_quant.sbatch`: where it downloads from and where it writes.

ENA files a run under ``vol1/fastq/<first six>/<a second directory>/<run>``, and the second
directory depends on how many digits the run number has: none for six, ``00`` and the last
digit for seven, ``0`` and the last two for eight, the last three for nine.  The script
always made ``0`` and the last two, which is right for eight digits only: ERR2060213, of
seven, is at ``ERR206/003``, not ``ERR206/013``.

The function is run under ``bash`` as written in the script, so it must stay bash-3.2
syntax (macOS's ``/bin/bash``).
"""
import os
import re
import shutil
import subprocess

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), os.pardir, "scripts",
                      "01_salmon_quant.sbatch")
BASH = shutil.which("bash")
pytestmark = pytest.mark.skipif(BASH is None or os.name == "nt", reason="needs bash")
ROOT = "https://ftp.sra.ebi.ac.uk/vol1/fastq/"


def _function(name):
    text = open(SCRIPT).read()
    m = re.search(r"^%s\(\) *\{.*?^\}$" % name, text, re.S | re.M)
    assert m, "no multi-line %s() { ... } in the script" % name
    return m.group(0)


def _ena(acc):
    r = subprocess.run([BASH, "-c", _function("ena") + '\nena "$1"', "_", acc],
                       capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr


@pytest.mark.parametrize("acc,where", [
    ("SRR123456", "SRR123/SRR123456"),                # six digits: no second directory
    ("ERR2060204", "ERR206/004/ERR2060204"),          # seven: 00 + the last digit
    ("ERR2060213", "ERR206/003/ERR2060213"),          # ... observed in ENA's fastq_ftp
    ("SRR24003390", "SRR240/090/SRR24003390"),        # eight: 0 + the last two
    ("SRR10134643", "SRR101/043/SRR10134643"),
    ("SRR123456789", "SRR123/789/SRR123456789"),      # nine: the last three
])
def test_the_ena_path_follows_the_run_number_length(acc, where):
    code, out, err = _ena(acc)
    assert code == 0, err
    assert out == ROOT + where


def test_a_run_number_of_another_length_is_an_error(tmp_path):
    code, out, err = _ena("SRR12345")
    assert code != 0 and "SRR12345" in err


def test_the_output_directory_is_per_sample_map(tmp_path):
    # two cohorts whose sample maps reuse donor names (ctrl1..ctrl5) must not share
    # quant/<donor>: the script skips a donor whose quant.sf exists, so the second cohort
    # silently kept the first one's
    text = open(SCRIPT).read()
    assert re.search(r'^OUTDIR=\$\{OUTDIR:-\$PWD/quant/\$\(basename "\$SAMPLE_MAP" \.csv\)\}$',
                     text, re.M)
    r = subprocess.run([BASH, "-c", 'SAMPLE_MAP=example/sample_map_GSE228458.csv; '
                        + re.search(r"^OUTDIR=.*$", text, re.M).group(0) + '; echo "$OUTDIR"'],
                       capture_output=True, text=True, cwd=str(tmp_path))
    assert r.stdout.strip() == os.path.join(os.path.realpath(str(tmp_path)), "quant",
                                            "sample_map_GSE228458")
