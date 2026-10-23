# Contributing

Contributions, bug reports, and feature requests are welcome.

## Reporting bugs or asking questions

Open an issue at
<https://github.com/charliekim97/isoform-dominance-pipeline/issues>. For a bug,
please include:

- the command you ran and the full error message,
- your operating system and Python version (`python --version`),
- the package version (`isoform-dominance --version`),
- a minimal example or input file if possible.

For general questions about usage, open an issue with the `question` label.

## Getting support

The GitHub issue tracker is the primary support channel. There is no separate
mailing list or chat; please open an issue (using the `question` label for usage
questions) and the maintainer will respond there.

## Development setup

```bash
git clone https://github.com/charliekim97/isoform-dominance-pipeline
cd isoform-dominance-pipeline
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q          # all tests are offline and must pass
```

`docs/index.md`, the documentation site's front page, is made from the README: after
editing the README, run `python scripts/make_docs_index.py` (`tests/test_docs_index.py`
fails until the two agree). To build the site as CI does, `pip install "mkdocs==1.6.1"` and
`mkdocs build --strict`.

## Submitting changes

1. Fork the repository and create a branch for your change.
2. Add or update tests; the suite must stay green and offline (no network in tests).
3. Keep the public CLI behaviour backward-compatible, or note the break clearly.
4. Open a pull request describing the change and the motivation. CI must pass: tests
   on Python 3.10–3.14 (Linux) and 3.12 (macOS, Windows), tests at the declared dependency
   floors, an install from the built wheel, the source distribution's own tests, `ruff`,
   the committed docs example, and a check that no commit credits an AI assistant; and,
   when the documentation or the paper changes, the strict build of the documentation site
   and the JOSS draft PDF.

## Scope

This tool is intentionally focused on the single-gene isoform-dominance question.
Genome-wide differential transcript usage is well served by existing packages
(`DEXSeq`, `DRIMSeq`, `satuRn`, `IsoformSwitchAnalyzeR`); proposals that broaden
scope toward those tools are unlikely to be merged. Improvements to grouping,
identifiability, statistics, documentation, and robustness are very welcome.

## Code of conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
