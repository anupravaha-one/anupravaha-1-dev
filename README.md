# Anupravaha-I

[![CI](https://github.com/anupravaha-one/anupravaha-1-dev/actions/workflows/build.yml/badge.svg)](https://github.com/anupravaha-one/anupravaha-1-dev/actions/workflows/build.yml)
[![contributors](https://img.shields.io/github/contributors/anupravaha-one/anupravaha-1-dev)](https://github.com/anupravaha-one/anupravaha-1-dev/graphs/contributors)
[![license: GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](LICENSE)

Educational, structured multi-block, finite-volume Navier–Stokes solver developed at the Indian Institute of Technology Kanpur. Anupravaha-I covers incompressible and weakly-compressible flow, conjugate heat transfer, turbulence (k–ε, RNG k–ε, k–ω, SST), natural and mixed convection, and liquid-metal MHD.

The 2018 release in this repository is a frozen reference implementation, used in graduate-level CFD coursework and as a teaching companion to the printed manual.

**Full documentation, theory, and 20+ worked tutorials live on the project website:** https://vatsalyasharma.gitlab.io/Anupravaha_1

## Layout

| Path | Contents |
|---|---|
| `APKMain.{c,h}` | Top-level driver |
| `Bc/` | Boundary-condition modules |
| `Cgns/` | CGNS I/O bindings and helpers |
| `Common/` | Shared utilities and data structures |
| `Geometry/` | Mesh / geometry handling |
| `Initialize/` | Field and parameter initialization |
| `Io/` | Input parsing, output writers, post-processing helpers |
| `Solver/` | Core algorithms and physics modules — Combustion, LMMHD, MPFV, Radiation, Solidification, Turbulence, … |
| `User/` | User-defined functions, boundary profiles, source terms |
| `CMakeLists.txt` | Top-level CMake build script |

## Building

Requirements:

- CMake ≥ 2.6
- A C compiler (GCC, Clang, MSVC)
- The CGNS library — https://cgns.github.io

Out-of-source build:

```bash
mkdir build && cd build
cmake ..
make -j
```

The executable lands under `build/bin/`. See the [Tutorials](https://vatsalyasharma.gitlab.io/Anupravaha_1/tutorials/) for end-to-end runs.

## Running

Anupravaha expects a project directory containing a `<case>.cgns` mesh. On first launch it creates `Settings/` and `SOLUTIONS/` subdirectories beside the executable that hold run-time parameters and converged fields.

## License

Released under the GNU General Public License v3.0 — see [`LICENSE`](LICENSE).

## Citation

If you use Anupravaha in academic work, please cite using the metadata in [`CITATION.cff`](CITATION.cff). GitHub renders this as a "Cite this repository" button on the repo home page.

## Contributing

Issues and pull requests are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md).
