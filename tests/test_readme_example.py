"""The README's worked example is what the code prints, from inputs in the repository.

The `identifiability` block runs the saved LEPR inputs of GENCODE 50 (Ensembl 116) in
`example/saved/`: every line it shows is the command's stdout, run here with every socket
refused, and it exits 0; with `--min-log2fc 0.1` the same run exits 3.  The command is read
from the README and run from the repository root with the paths as written, because the
header line repeats the `--inputs` path as given.  The `annotate` block's transcript counts
are the saved config's.
"""
import json
import pathlib
import re
import shlex
import socket

import pytest

from isoform_dominance import cli

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
CONFIG = ROOT / "example" / "saved" / "LEPR.gencode50.config.json"


def _blocks():
    return re.findall(r"```bash\n(.*?)```", README.read_text(encoding="utf-8"), re.S)


def _identifiability_block():
    """``(argv, output lines, the closing note)`` of the README's identifiability example."""
    block = next(b for b in _blocks()
                 if b.startswith("isoform-dominance identifiability --config example/saved/"))
    command, rest = block.split("\n#", 1)
    lines = ("#" + rest).rstrip("\n").split("\n")
    assert all(x.startswith("# ") for x in lines), lines
    argv = shlex.split(command.replace("\\\n", " "))
    assert argv[0] == "isoform-dominance"
    return argv[1:], [x[2:] for x in lines[:-1]], lines[-1][2:]


@pytest.fixture
def offline(monkeypatch):
    """No connection can be made: the run must use nothing but the saved inputs."""
    def refuse(*a, **k):
        raise AssertionError("the README example tried to reach the network")
    monkeypatch.setattr(socket, "socket", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket, "getaddrinfo", refuse)


def test_every_line_of_the_identifiability_block_is_what_the_command_prints(
        offline, monkeypatch, capsys):
    argv, shown, note = _identifiability_block()
    assert "--inputs" in argv and "--min-log2fc" in argv
    monkeypatch.chdir(ROOT)
    code = cli.main(argv)
    assert capsys.readouterr().out.rstrip("\n").split("\n") == shown
    assert code == 0
    assert note.startswith("exit status 0;")
    # the note's second half: a smaller effect is not resolved, exit 3
    m = re.search(r"with --min-log2fc ([0-9.]+) the same run reports NOT resolved and exits 3",
                  note)
    assert m
    argv[argv.index("--min-log2fc") + 1] = m.group(1)
    assert cli.main(argv) == 3
    assert "EFFECT SIZE: |log2FC| %.3f NOT resolved" % float(m.group(1)) in capsys.readouterr().out


def test_the_annotate_block_shows_the_saved_config():
    block = next(b for b in _blocks() if b.startswith("isoform-dominance annotate --gene LEPR"))
    shown = dict(re.findall(r"^#\s+(iso_\w+)\s*:\s*(\d+) transcripts", block, re.M))
    cfg = json.loads(CONFIG.read_text())
    assert {g: int(n) for g, n in shown.items()} == {g: len(v) for g, v in cfg["groups"].items()}
    canonical = re.search(r"^#\s+(iso_\w+)\s*:.*canonical\)", block, re.M).group(1)
    assert [c["transcripts"] for c in cfg["_clusters"] if c["contains_canonical"]] == [
        cfg["groups"][canonical]]
