# Contributing to Anupravaha-I

Thank you for taking an interest. Anupravaha-I is primarily an educational solver, so contributions that **clarify code, fix bugs, improve tutorials, or extend documentation** are especially welcome.

## Quick start

```bash
# Fork the project on GitLab (top-right "Fork" button), then:
git clone https://gitlab.com/<your-username>/Anupravaha_1.git
cd Anupravaha_1/src/AnuPravahaStructured
mkdir build && cd build
cmake .. && make -j
```

Branch off `main`, commit your change, push to your fork, and open a Merge Request against `vatsalyasharma/Anupravaha_1`.

## Verified build environment

The current `main` builds clean (zero errors, GCC's expected 22 buffer-size warnings) with the toolchain below on a stock Ubuntu 22.04 LTS / WSL2 setup. Newer versions almost certainly work; older usually do — open an MR or issue to extend this matrix.

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

## Style

The existing C code uses K&R-style braces, 4-space indentation, and `APK`-prefixed function names. Please match the surrounding code in any file you touch rather than imposing a new style.

## Merge requests

- One topic per MR — easier to review and revert if needed
- Reference any related issue with `Closes #N`
- Describe what you changed, why, and what you tested
- Expect review before merge; we may ask for changes

## Reporting issues

Use the GitLab Issues tab. Helpful issue reports include:

- Anupravaha version / commit
- OS, compiler, CGNS version
- The `<case>.cgns` mesh (or a minimal reproducer) if the issue is solver-side
- The `Settings/` directory contents for the run that failed
- What you expected vs. what happened

## Maintainer checklist (for project administrators)

Before flipping the project to public visibility on GitLab:

- [ ] Settings → General → Visibility: **Public**
- [ ] Settings → Repository → Protected branches: protect `main`, "Allowed to push" = **No one**, "Allowed to merge" = Maintainers
- [ ] Settings → Repository → Protected tags: protect `v*`
- [ ] Settings → Merge requests → Approvals required: **≥ 1**, "Prevent approval by author" enabled
- [ ] Settings → Account → Two-Factor Authentication enabled on every Maintainer account
- [ ] Settings → Repository → Repository mirroring: configure a backup mirror
- [ ] Optional: add a `CODEOWNERS` file to require specific reviewers for `src/Solver/` etc.
- [ ] Optional: add a CI job that builds the C code on every MR, so broken changes can't merge
