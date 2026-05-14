# Contributing to Anupravaha-I

Thank you for taking an interest. Anupravaha-I is primarily an educational solver, so contributions that **clarify code, fix bugs, improve tutorials, or extend documentation** are especially welcome.

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md) — by participating you agree to abide by its terms. For **security-sensitive issues** see [`SECURITY.md`](SECURITY.md).

## Quick start

```bash
# Fork the project on GitHub (top-right "Fork" button), then:
git clone https://github.com/<your-username>/anupravaha-1-dev.git
cd anupravaha-1-dev
mkdir build && cd build
cmake .. && make -j4
```

Branch off `main`, commit your change (with `-s` for DCO sign-off), push to your fork, and open a Pull Request against `anupravaha-one/anupravaha-1-dev`.

## Verified build environment

The current `main` builds clean (zero errors, GCC's expected 22 buffer-size warnings) with the toolchain below on a stock Ubuntu 22.04 LTS / WSL2 setup. Newer versions almost certainly work; older usually do — open a PR or issue to extend this matrix.

| Tool | Tested version | Required | Notes |
|---|---|---|---|
| Ubuntu | 22.04.5 LTS (jammy) | any modern Linux | Native or via WSL2 on Windows 10/11 |
| GCC | 11.4.0 | yes — C, C++, Fortran | `gcc`, `g++`, `gfortran` |
| CMake | 3.22.1 | ≥ 2.6 (per `CMakeLists.txt`) | |
| GNU Make | 4.3 | any recent | |
| GNU ld (binutils) | 2.38 | any recent | |
| `libhdf5-dev` | 1.10.7 | yes | Headers used by the vendored CGNS sources |
| `libcgns-dev` | 3.4.0-3 | optional | CGNS is also vendored under `Cgns/` — system package not strictly required |
| Git | 2.34.1 | any recent | For cloning |

Last verified clean build: **2026-05-13** — `cmake ..` + `make -j4` produces `build/bin/Anupravaha_1` (the executable) and `build/bin/libAPK.so` (the user-relink library). Build wall-time ~5–10 min on a fast disk, longer on `/mnt/c` or `/mnt/d` (WSL 9p mount).

## What kinds of changes are welcome

- Bug fixes (please link to a reproducer or describe how you triggered the bug)
- Tutorial additions and corrections (the `docs/tutorials/` tree on the website)
- New physics modules under `Solver/`, with a matching theory page and tutorial
- Performance improvements that don't reduce readability
- Documentation in any form — README clarifications, code comments, theory pages

## What to avoid

- Wholesale reformatting of existing files (large whitespace-only diffs make review hard)
- Changes that drop support for the canonical tutorials without a clear reason
- Adding heavy external dependencies — Anupravaha-I deliberately stays close to plain C + CMake + CGNS

## Tests

We aim to follow the standard for mature CFD codes: **every bug fix gets a regression test that fails before the fix and passes after**; **every new feature ships with at least one demonstration case**. See [`tests/README.md`](tests/README.md) for the layout and how to add a case.

For bug-fix PRs, please include:

- A `.cgns` mesh (or a pointer to an existing tutorial mesh) that reproduces the bug
- A recorded `.inp` input file (see [`tests/README.md`](tests/README.md))
- A short note in the PR description explaining what the failure mode was *before* your fix and what it is *after*

CI builds the code on every PR. Once at least one regression case is committed, CI additionally runs it as a smoke test. Submissions that don't build will be rejected automatically.

If you genuinely can't write a test for your change (e.g., a typo fix in a comment, or a build-system tweak), say so explicitly in the PR description and a reviewer may waive the requirement.

## Sign your work (Developer Certificate of Origin)

Anupravaha-I uses the [Developer Certificate of Origin](https://developercertificate.org/) to keep the chain of contribution legally clean. **Every commit must carry a `Signed-off-by:` line:**

```bash
git commit -s -m "Fix off-by-one in APKFindMaxDiff"
```

This is *not* a copyright assignment — it's a one-line statement that you have the right to contribute the code under the project's GPL-3.0 license (the full DCO text is at the link above; please read it once).

If you forget, amend the last commit with `git commit --amend -s`, or sign off a batch with `git rebase --signoff <upstream>`. The trailer ends up looking like:

```
Fix off-by-one in APKFindMaxDiff

Signed-off-by: Jane Doe <jane@example.com>
```

CI will reject PRs whose commits aren't all signed off.

## Style

The existing C code uses K&R-style braces, 4-space indentation, and `APK`-prefixed function names. Please match the surrounding code in any file you touch rather than imposing a new style. We deliberately do **not** auto-format the repository — large mechanical reformatting PRs make review impossible and hide real changes. Reformat only the lines you're already touching, never adjacent code.

## Pull requests

A pull request template is auto-loaded when you open a new PR — please fill it in. Summary of expectations:

- One topic per PR — easier to review and revert if needed
- **All commits signed off** (`git commit -s`) — see [DCO](#sign-your-work-developer-certificate-of-origin)
- A test added or updated, **or** an explicit note explaining why none is needed
- Reference any related issue with `Closes #N`
- Describe what you changed, why, and what you tested
- Expect **1 maintainer approval** before merge; we may ask for changes

## Reporting issues

Use the GitHub Issues tab — bug and feature templates are auto-loaded when you click *New issue*. Helpful bug reports include:

- Anupravaha version / commit
- OS, compiler, CGNS version
- The `<case>.cgns` mesh (or a minimal reproducer) if the issue is solver-side
- The `Settings/` directory contents for the run that failed
- What you expected vs. what happened

## Maintainer checklist (for project administrators)

Before flipping the project to public visibility on GitHub:

- [ ] **Settings → General → Danger Zone → Change visibility**: *Public*
- [ ] **Settings → Branches → Add branch protection rule** for `main`: require PR reviews before merging, require status checks (Build) to pass, restrict direct pushes, disallow force pushes
- [ ] **Settings → Tags → New rule** protecting `v*`: only Maintainers can create or push these
- [ ] Branch protection rule → *Require approvals*: **≥ 1**, *Require review from Code Owners*, *Dismiss stale pull request approvals when new commits are pushed*
- [ ] **Settings → Password and authentication** → enable two-factor authentication on every Maintainer account
- [ ] **Settings → Code security and analysis** → enable Dependabot alerts, Secret scanning, Private vulnerability reporting
- [ ] `.github/CODEOWNERS` is already shipped — review it for completeness and update the GitHub username if it differs from `@vatsalyasharma`
- [ ] CI workflow `.github/workflows/build.yml` is already shipped — verify it passes on a sample PR
- [ ] Optional: configure off-platform backup mirroring via a separate GitHub Actions workflow if desired
