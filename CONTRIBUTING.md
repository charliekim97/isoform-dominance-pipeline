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

## Submitting changes

1. Fork the repository and create a branch for your change.
2. Add or update tests; the suite must stay green and offline (no network in tests).
3. Keep the public CLI behaviour backward-compatible, or note the break clearly.
4. Open a pull request describing the change and the motivation. CI
   (Python 3.10–3.12) must pass.

## Scope

This tool is intentionally focused on the single-gene isoform-dominance question.
Genome-wide differential transcript usage is well served by existing packages
(`DEXSeq`, `DRIMSeq`, `satuRn`, `IsoformSwitchAnalyzeR`); proposals that broaden
scope toward those tools are unlikely to be merged. Improvements to grouping,
identifiability, statistics, documentation, and robustness are very welcome.

## Code of conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
