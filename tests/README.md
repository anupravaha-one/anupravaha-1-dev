# Tests

Regression tests for Anupravaha-I. Each case is a self-contained directory that the CI pipeline runs to verify the solver still produces correct output for a known scenario.

## Layout

```
tests/
├── README.md                    # this file
└── <case-name>/
    ├── <case>.cgns              # input mesh (or a symlink into ../../docs/assets/tutorial_meshes/ if shared)
    ├── <case>.inp               # recorded interactive responses (see "Recording an .inp" below)
    ├── reference/
    │   ├── residuals.dat        # expected residual history
    │   └── midplane.dat         # expected key flow quantity (varies per case)
    └── tolerance.yaml           # numerical tolerances for comparison
```

## Why this layout

Anupravaha asks for run parameters interactively (convection scheme, time stepping, pressure solver, etc.). A test case is **fully reproducible** when its `.inp` file deterministically answers every prompt — the solver then runs without any human in the loop. The `reference/` directory holds the known-good output, and `tolerance.yaml` declares how much numerical drift is acceptable.

## Recording an `.inp` for a new case

The `.inp` file is just a verbatim transcript of the answers you'd type interactively. The simplest recipe:

```bash
# 1. Run the solver manually, answering each prompt
./Anupravaha_1 mesh.cgns sol.cgns

# 2. The exact answers you typed end up in Settings/coefficients.dat in the order they were asked
cp Settings/coefficients.dat tests/<case-name>/<case-name>.inp

# 3. From now on, replay with stdin redirection:
./Anupravaha_1 mesh.cgns sol.cgns < tests/<case-name>/<case-name>.inp
```

## Adding a new regression test

1. Pick a tutorial case from the website's `docs/tutorials/`.
2. Create `tests/<case-name>/` and copy the mesh in.
3. Record the `.inp` as above.
4. Run once to convergence, copy the residual history and a key line probe into `reference/`.
5. Write `tolerance.yaml` — relative tolerance for residuals (typically `1e-6`), absolute tolerance for flow quantities (case-specific).
6. Update `.gitlab-ci.yml` to add a job for the new case, *or* leave it for the nightly suite if it takes more than a few minutes.

A `tolerance.yaml` template:

```yaml
residuals:
  rel_tol: 1.0e-6
flow_quantities:
  midplane_u:
    abs_tol: 1.0e-4
```

## Running tests locally

```bash
# Build once
cd ../  &&  mkdir -p build  &&  cd build  &&  cmake ..  &&  make -j4

# Then for each case
cd ../tests/<case-name>
../../build/bin/Anupravaha_1 <case>.cgns sol.cgns < <case>.inp
# Comparator: TBD (a small Python driver will land in tests/compare.py)
```

## CI

The pipeline in [`.gitlab-ci.yml`](../.gitlab-ci.yml) currently runs a **build** job on every push to a MR branch. Once at least one regression case is committed here, a **smoke-test** job will be added that runs the fastest case (target: under 60 s). The full regression suite will run nightly on `main` (planned).

## Coming soon

- `tests/compare.py` — small driver that runs a case, parses CGNS output, applies `tolerance.yaml`, exits non-zero if any tolerance is violated.
- First real case: `tests/lid-driven-cavity/` — once a maintainer records the `.inp`.
