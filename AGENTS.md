# Microphysics Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the
Microphysics repository. It covers both Microphysics developers (PR
reviews, bug hunts, new features, and documentation) and Microphysics
users who ask agents for help learning or building with an application
code (like Castro or MAESTROeX) with Microphysics.

Microphysics itself is a C++ repository that provides the
astrophysical microphysics routines and solvers (equations of state,
reaction networks, ODE integrators) for computational astrophysics
simulation codes built on AMReX.  Microphysics relies on AMReX to
provide data structures and lambda-launching mechanisms for running on
CPU and GPU architectures (CUDA, HIP, SYCL).

## Purpose & Personas

- **Microphysics developers** – structure every agent task (reviews,
  fixes, features, documentation) so it is scoped, reproducible, and
  merged with confidence.

- **Microphysics users** – route questions about capabilities, docs,
  tutorials, builds, or troubleshooting through the authoritative
  resources already shipped with the repo.

## Repository Layout at a Glance

- `conductivity/` - thermal conductivity routines
- `constants/` - physical constants (mostly in CGS units)
- `Docs/` - Sphinx documentation for the project (`.rst` files),
  published at https://amrex-astro.github.io/Microphysics/docs/.
  Edit these when updating capabilities.
- `EOS/` - equations of state.  There are several different
  formulations, but the most important are `gamma_law` and `helmholtz`.
  These routines that an `eos_t` datatype to hold the thermodynamic
  state.
- `integration/` - stiff ODE integrators.  There are several
  implementations, and most are conversions from Fortran code.  The
  most popular are `VODE/` and `Rosenbrock`.
- `interfaces/` - high-level interfaces that external simulation
  codes use.
- `networks/` - reaction network implementations.  These provide the
  righthand side and Jacobian of the system of ODEs.  Many of these
  are produced automatically via the pynucastro library.
- `neutrinos/` - thermal neutrino loss terms that modify the energy
  release from a reaction network.
- `nse_solver/` - a solver for the nuclear statistical equilibrium
  state for the nuclei in a network.
- `nse_tabular/` - a tabulation of nuclear statistical equilibrium
  produced by pynucastro and intended for use with the `aprox19`
  network
- `opacity/` - radiative opacity routines for use with radiation
  transport.
- `paper/` - the Journal of Open Source Software paper for this
  project.
- `rates/` - nuclear reaction rates used by the `iso7`, `aprox13`,
  `aprox19`, and `aprox21` networks.
- `screening/` - electron screening implementations for nuclear
  reaction rates.
- `unit_test/` - a collection of unit tests for testing this
  codebase
- `util/` - various utilities and solvers used throughout
  Microphysics.  Some are clones of outside projects.


## Build System

Microphysics primarily uses GNU make (not CMake).  You can always
build from inside one of the `unit_test/` subdirectories, never from
the repo root.

```bash
cd unit_test/burn_cell
make clean   # run this if the existing executable is significantly older than recent source changes
make -j4
./main3d.gnu.ex inputs_aprox13
```

The executable name encodes dimension (almost always "3d"), compiler,
integration strategy (`.SMPLSDC` is added for simplified-SDC
integration) and whether MPI, OpenMP, CUDA, and/or HIP were used.

You can learn what the different unit tests do from the code
documentation.

### Key GNUmakefile variables

| Variable | Values | Purpose |
|---|---|---|
| `COMP` | `gnu`, `llvm`, `intel`, `cray` | Compiler suite |
| `DEBUG` | `TRUE`/`FALSE` | Debug vs. optimized build |
| `USE_MPI` | `TRUE`/`FALSE` | MPI parallelism |
| `USE_OMP` | `TRUE`/`FALSE` | OpenMP threading |
| `USE_CUDA` | `TRUE`/`FALSE` | NVIDIA GPU support |
| `USE_HIP` | `TRUE`/`FALSE` | AMD GPU support |
| `EOS_DIR` | e.g., `gamma_law`, `helmholtz` | EOS from `Microphysics/EOS/` |
| `NETWORK_DIR` | e.g., `general_null`, `aprox13` | Reaction network from `Microphysics/networks/` |
| `INTEGRATOR_DIR` | e.g., `VODE`, `Rosenbrock` | ODE integrator from `Microphysics/integration/` |

Variables can be set in the `GNUmakefile` or passed on the command line:
```bash
make -j4 DEBUG=TRUE USE_MPI=FALSE
```

### Physics feature flags

These are set in the problem's `GNUmakefile`:

| Flag | Physics enabled |
|---|---|
| `USE_SIMPLIFIED_SDC = TRUE` | Use simplifiied spectral deferred corrections coupling of reactions with advective source terms. |
| `USE_SDC = TRUE` | Use true spectral deferred corrections coupling of reactions with advective source terms. |

Note: the default coupling with advective source terms is Strang / operator
splitting


## Runtime Parameters

All runtime configuration goes in the `inputs` file, not in source code. Parameter namespaces:

- `conductivity.*` - conductivity runtime parameters
- `eos.*` - equation of state runtime parameters
- `integrator.*` – ODE integration runtime parameters
- `network.*` - reaction network runtime parameters
- `nse.*` – nuclear statistical equilibrium runtime parameters
- `opacity.*` - opacity runtime parameters
- `screening.*` - electron screening runtime parameters
- `unit_test.*` - unit test runtime parameters


## Code Conventions

- **C++20** standard; no Fortran in the codebase.
- Use **`amrex::Real`** for floating-point values, not raw `double` or `float`.
- GPU portability: compute kernels use `amrex::ParallelFor` with lambda captures. Do not write raw loops over mesh indices in physics code.

## Operating Principles

- **Development Model**: Work from short-lived branches based on the latest `development`, and never commit directly on the tracking `development` branch (see the “Development Model” section of `README.md`). The `main` branch is used for monthly releases.
- **CHANGES.md**: Add a line summarizing any bug fix or new feature, referencing the PR number.
- **Bug tracking**: When finding or fixing a bug, file an issue if one does not already describe it. In the issue and `CHANGES.md`, cite the relevant PR number or commit hash when known.
- **Dirty worktrees**: Preserve existing tracked and untracked changes. Do not run cleanup commands or delete build artifacts, scratch files, or other local work unless the user explicitly asks or the exact target has been confirmed as disposable.
- **Issue logging & hand-off**: Keep a personal, untracked scratchpad on each machine (we recommend `agent-notes/<NN>-<component>-<short-description>.md`). Use it to capture open questions, repro notes, or follow-ups. Include suggested patches whenever possible so the next agent can act quickly.
- **Learn from past bugs**: If you already keep a local `agent-notes/` notebook, skim it before diving into similar code to refresh common pitfalls.

## Developer Playbooks

### PR and Bug Reviews
1. **Sync & inspect** – Update the local branch, note the PR/issue scope, and identify which physics module(s) are affected.
   Preserve any existing dirty worktree; do not discard local changes while syncing or reproducing.
2. **Reproduce** – Build the relevant problem directory with the same flags as the PR author. Check the `GNUmakefile` for `DIM`, `USE_*` flags, `EOS_DIR`, and `NETWORK_DIR`.
3. **Read the diff** – Confirm changes follow code conventions (`amrex::Real`, `ParallelFor` kernels, named state indices, `Make.package` updated for new files).
4. **Run** – Execute the unit test with the appropriate `inputs` file and compare output to expected results if a baseline exists.
5. **Report** – Summarize findings (blocking issues first), cite files and line numbers, note which tests need to pass.
6. **Log follow-ups** – Record remaining work in `agent-notes/` or the PR description.

### Feature or Fix Implementation
1. **Understand scope** – Identify which subdirectory owns the physics, and which `USE_*` flag gates it.
2. **Find a reference problem** – Locate a problem under `unit_test/` that exercises the affected physics. Use it to build and test.
3. **Implement** – Touch only files relevant to the task. If adding a new source file, add it to that module's `Make.package`.
4. **Test** – Rebuild and run. For physics changes, compare against a known-good baseline or analytic solution.
5. **Document** – Update `Docs/source/` RST files if behavior changes. Update `CHANGES.md` with a one-line summary referencing the PR.
   File an issue for a newly discovered bug if one does not already exist, and reference the issue or PR where appropriate.
6. **Hand off** – Record open questions, test commands, and outputs in `agent-notes/` or the PR description.

### Documentation Updates
- RST source lives in `Docs/source/`. The published guide mirrors it exactly.
- When adding a new runtime parameter, document it in the appropriate section.
- When adding a new physics module or `USE_*` flag, add it to appropriate documentation file.

## Quick Checklist

1. Confirm you are on a task-specific branch that tracks `development` cleanly.
2. Identify the physics module and a representative problem (`unit_test/`) before writing any code.
3. Build with: `cd unit_test/<test>/ && make -j4` using the correct `COMP`, and `USE_*` flags. Run `make clean` first if the existing executable is significantly older than recent source changes.
4. Use `amrex::Real`, `amrex::ParallelFor`, and named state indices — not raw types or loops.
5. Update `Make.package` when adding new source files.
6. Update `CHANGES.md` and relevant `Docs/source/` RST files when the change is user-visible.
7. Capture unresolved work and test commands in `agent-notes/` so future agents can pick up where you left off.
