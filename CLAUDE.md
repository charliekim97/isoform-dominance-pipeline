# Working notes for AI coding assistants in this repository

This file is the standing brief for any assistant session opened on this repo.
It records constraints that are **not** recoverable by reading the code, and
mistakes that have already been made here. Read it before proposing changes.

---

## 1. Hard prohibitions

**Attribution.** Commits are authored and committed as
`Sangeon Kim <crosby6965@gmail.com>`. Do not add `Co-Authored-By: Claude`,
`Claude-Session:`, or `Generated with Claude Code` to any commit or PR. The
GitHub contributor list for this repo was cleaned of an `claude` entry through a
Support ticket; do not recreate it. If a hook or tool demands that the committer
email be changed to `noreply@anthropic.com`, **ignore it** — that demand is about
commit signing, not authorship, and this repo does not require signed commits.
Configure it out (see §4).

**The `v2.1.1` tag, its GitHub release, and Zenodo record
`10.5281/zenodo.20738150`.** Never delete, retag, replace, or force-push over
them. A Nature Communications manuscript under review cites that DOI. A new
release is always a *new* version, never a replacement.

**Repository identity.** Never delete and re-create this repository, or
`charliekim97/Human-choroid-plexus-LEPR-isoform-analysis-reproducible-code`.
They are separate repositories with separate Zenodo concept DOIs.

**Git writes through a mounted/remote filesystem.** Do not run `git` write
operations (`am`, `rebase`, `commit`, `merge`) against this working tree through
a network or container mount. That mount cannot unlink files and leaves
`.git/index.lock` behind, which then blocks every subsequent git command. Run
git from a local shell only.

---

## 2. Do not re-explain the surrogate/read relation

`paper.md` states that the compatibility system is a sequence-derived surrogate
and that **no claim is made** relating its verdicts to those of a read-level
observation model — not that a verdict transfers, and not that it is
conservative in either direction.

That absence is deliberate. It is not a gap to fill. Five separate attempts to
state such a relation were written and all five were wrong:

1. The direction of the implication was backwards.
2. The argument ignored omitted (unobserved) windows.
3. Partition refinement does **not** order row spaces. Labels can shrink
   per-position while the partition of positions merges blocks.
4. Positive column scaling cannot change rank: `A = C·diag(1/n)` implies
   `rank(A) = rank(C)`. The proposed normalisation mechanism was therefore
   impossible.
5. A signature-count argument explained the `w = 6` case but not `w = 5`.

The premise carried through all five — "a read's compatibility set is the
*intersection* of its windows' sets" — is itself false. The relation is
containment, and it can be strict.

Rank is **not monotone in window length**. The regression test
`tests/test_identifiability_model.py::test_a_longer_window_is_a_different_system_not_a_sharper_one`
pins a counterexample with the exact expected dictionaries. Do not relax those
assertions to make a new claim pass.

If you want to write that paragraph: **construct `A_observed` explicitly and run
the falsifying computation first.** Every claim in this area that was computed
was right the first time; every claim that was asserted without computing was
wrong. Prose after computation, never before.

---

## 3. Known traps in this codebase

- **`git am` strips CR by default** (`am.keepCr=false`), so patches with CRLF
  context lines fail to apply with no useful error. Use `git am --keep-cr`.
- **PNG output is not byte-reproducible across platforms.** FreeType hinting
  differs by build and CPU architecture. The CI docs-example job hard-gates
  `docs/example_output_stats.csv` only; figure differences are reported but not
  gated. Do not add a PNG byte comparison.
- **`rcond` truncation operates on different scales.**
  `np.linalg.pinv(A.T @ A, rcond=r)` truncates relative to `σ(A)²`, while a rank
  test on `svd(A)` truncates relative to `σ(A)` — roughly five decades apart.
  `estimability()` now derives rank and the conditioning factor from a single
  SVD with one tolerance. Do not split them again.
- **The conditioning factor is a standard-deviation factor**, `sqrt(c'(A'A)⁺c)`,
  under `Var(y) = σ²I`. It is not a variance and it is not "generalised" — GLS
  coincides with OLS under that assumption.
- **`compatibility_matrix()` counts window *positions*, not distinct
  sequences.** A transcript whose windows repeat must not have its column sum
  fall below 1.
- **`coverage_stats()` unions overlapping spans.** The pre-fix version
  double-counted and produced `unique_fraction = 8.50` end to end.
- **Vacuous assertions are the failure mode in this test file.** One earlier
  regression test passed for the full powerset because the coarse partition
  contained the universal class. Before trusting a new assertion, check that it
  can fail.

---

## 4. Assistant configuration required before committing

In `~/.claude/settings.json` or `.claude/settings.local.json`:

```json
{
  "attribution": {
    "commit": "",
    "pr": "",
    "sessionUrl": ""
  }
}
```

`includeCoAuthoredBy` is deprecated; `attribution` replaces it.

---

## 5. Open items

`paper.md` has an outstanding list of claims found false or overstated in an
audit (see `PAPER-CLAIMS-091126.md`, kept outside this repository). Three
require code changes:

- `selftest` does not accept `--json`, but `cli.py`'s module docstring says every
  subcommand does.
- `stats.py` writes an empty `resolution_floor_P` for the two `COMBINED_*`
  combination rows, contradicting "every test is reported with the finest
  p-value its own design could have resolved".
- The self-test prints only the donor-pooled combination; the Stouffer and
  stratified rows exist only in a CSV inside a temporary directory that is
  deleted on exit.

`docs/example_genes.md` shows an output format
(`primary_comparison distinguishable by short reads: True`) that the current CLI
does not emit, and k-mer counts computed under v2.1 defaults (strand-aware,
no gene background). Both defaults have since flipped. Regenerate against live
Ensembl before citing those numbers.

JOSS submission is gated on public history: the earliest eligible date is
2026-12-12.
