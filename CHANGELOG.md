# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Local development instructions in the installation docs: the compilers and libraries `packaging/build_ipopt.sh` needs on Linux, macOS and Windows, and the Docker/cibuildwheel route for building a wheel without installing any of them.

### Changed

* One package again: the IPOPT binding is built into `compas_cra` itself as `compas_cra._native`, and `compas_cra` ships as platform wheels for CPython 3.9-3.13 on Windows, macOS (Apple Silicon and Intel) and manylinux (x86_64 and aarch64). `pip install compas_cra` installs exactly one distribution, solver included.
* The build backend is scikit-build-core, which drives the CMake project at the repository root. It reads no `requirements.txt`, so the dependencies now live in `pyproject.toml`.
* The cibuildwheel settings live in `[tool.cibuildwheel]` in `pyproject.toml` rather than in workflow environment variables, so `cibuildwheel` reproduces the CI wheel locally with no environment set. Only the runner-dependent parts (workspace paths, the macOS deployment target, the MSYS2 location) stay in `pipeline.yml`.
* One CI pipeline for every push and every release: `build.yml` (pushes and PRs to main) and `release.yml` (`v*` tags) both call the reusable `pipeline.yml`, which differs only by its `publish` input. A release therefore runs the exact chain that was already green on main — IPOPT from source, all five platforms, the smoke tests, the test suite against a built wheel and the docs — and then publishes it. `wheels.yml` and `docs.yml` are gone, folded into that pipeline.
* `packaging/check_release.py` now checks for the `compas_cra._native._core` extension rather than a vendored `ipopt` executable, which the statically linked build no longer ships, and additionally requires every supported CPython on every platform. Verified against the artifacts of the PR build: it passes the complete set and rejects a partial one.
* `invoke release` is defined in `tasks.py` instead of coming from `compas_invocations2`, without its local `python -m build` step: building a wheel here needs a staged IPOPT tree, and the published artifacts are built by CI on all five platforms regardless.
* `packaging/build_ipopt.sh` no longer trusts the third-party source downloads. The `get.ASL` and `get.Mumps` scripts fetch their tarballs with a bare `curl -L -O`, with no integrity check and no retry, so a stalled transfer lands a short file that is then gunzipped anyway — which is how an aarch64 build failed with `gzip: MUMPS_5.8.2.tar.gz: not in gzip format` after receiving 54 kB of a 4.3 MB tarball. The script now points `CURL_HOME`/`WGETRC` at a config that abandons and retries a stalled transfer, and checks that `ThirdParty/ASL/solvers` and `ThirdParty/Mumps/MUMPS` actually exist, re-running the `get.*` script that came up empty (a second `coinbrew fetch` will not: it only runs `get.*` when the clone's revision changed).

### Removed

* Removed the separate `compas_cra_native` distribution and its `native/pyproject.toml`; the binding source stays in `native/`, built by the root `CMakeLists.txt`. An install from source (`pip install .`, editable included) now compiles the extension and needs `packaging/build_ipopt.sh` to have run first.
* Removed `requirements.txt`, `requirements-dev.txt` and `requirements-viz.txt`; `pip install -e ".[dev]"` replaces `pip install -r requirements-dev.txt`.
* Removed the `build_ghuser_components` task and its `ghuser` configuration from `tasks.py`; the paths pointed at `src/compas_notebook`, and compas_cra ships no Grasshopper components.

## [0.7.2] 2026-08-23

### Added

### Changed

* One ordered release pipeline in `release.yml`: native wheels build first (IPOPT compiled from source), the pure package is tested against the just-built native wheel, then publish native → publish cra → GitHub release → docs. `native.yml` is build-only push CI; `docs.yml` builds on main pushes only; the lint job runs without installing the package.

### Removed

## [0.7.1] 2026-08-23

### Added

### Changed

* `compas_cra_native` 0.1.1: the audit-hardened binding (sparsity index range checks, transient-callback-exception recovery, automatic quasi-Newton fallback, wall-time status) actually ships to PyPI — 0.1.0 predated those fixes and `skip-existing` had kept re-uploads out.
* Rhino scripts pin `compas_cra>=0.7.0` so cached ScriptEditor environments upgrade off the removed executable path.

### Removed

## [0.7.0] 2026-08-23

### Added

* Added `cra_penalty_problem` and `rbe_problem` NLP formulations with exact analytic derivatives, and the matching native solvers `cra_penalty_solve_native` and `rbe_solve_native`, validated against the pyomo implementations before those were removed (max force difference 1.9e-14 for RBE, 7.1e-9 for the penalty formulation).

### Changed

* `cra_solve`, `cra_penalty_solve` and `rbe_solve` are now the in-process native solvers. `compas_cra` is a pure Python package again — one universal wheel — with the compiled solver coming from the `compas_cra_native` dependency.
* The native wheel workflow builds the Linux wheels with cibuildwheel inside the manylinux containers, in the style of compas_cgal.

### Removed

* Removed the bundled IPOPT executable, the pyomo formulations and the pyomo dependency. There is no solver executable anywhere anymore: solving happens in-process through the nanobind extension, which also removes the Windows antivirus false-positive problem structurally.

## [0.6.6] 2026-08-23

### Added

### Changed

* `compas_cra_native` is now a dependency of `compas_cra` (on Python < 3.14), so `pip install compas_cra` — or a single `# r: compas_cra` line in Rhino — brings the in-process solver along automatically.

### Removed

## [0.6.5] 2026-08-23

### Added

### Changed

* Hardened the native binding after a three-stage independent audit: sparsity index range checks (out-of-range values from Python now raise instead of crashing), a solve that converges despite a transient callback exception keeps its solution, `hessian_approximation=limited-memory` is set automatically when no Hessian callback is given, and `Maximum_WallTime_Exceeded` is reported by name.
* The native backend no longer accepts iteration-capped points; only near-converged endings (`Restoration_Failed`, `Search_Direction_Becomes_Too_Small`) qualify for the feasible-point acceptance check, matching the executable path's failure behavior.
* The Rhino example scripts solve with `cra_solve_native` (in-process, no executable involved).
* macOS arm64 native wheels are built on macOS 14, so they install on macOS 14+.

### Removed

## [0.6.4] 2026-08-23

### Added

* Added `compas_cra.nlp`, a solver-agnostic sparse NLP layer, and `compas_cra.equilibrium.cra_problem`, the CRA optimisation problem formulated directly in numpy/scipy with exact analytic gradient, Jacobian and Lagrangian Hessian (no pyomo involved).
* Added `compas_cra_native` (in `native/`): IPOPT + MUMPS compiled into a Python extension module with nanobind, so CRA problems solve in-process — no bundled executable, no subprocess, no `.nl` files. Results match the pyomo + executable path on the test suite (bit-identical on the cube fixtures, < 1e-5 relative force difference on the arch).
* Added `cra_solve_native`, a drop-in alternative to `cra_solve` using the binding.

### Changed

### Removed

## [0.6.3] 2026-08-23

### Added

* Added the `COMPAS_CRA_IPOPT` environment variable to override the solver location, for machines where antivirus or application-control policies block the bundled binary.
* The Windows `ipopt.exe` now carries a version resource (product, publisher, version, license), and the release workflow supports Authenticode signing via Azure Trusted Signing when the signing secrets are configured.

### Changed

### Removed

## [0.6.2] 2026-08-23

### Added

### Changed

* Lowered `requires-python` to `>= 3.9` so the package installs into Rhino 8's bundled CPython 3.9.

### Removed

## [0.6.1] 2026-08-23

### Added

### Changed

### Removed

## [0.6.0] 2026-08-23

### Added

* Added a bundled, statically linked IPOPT executable to the platform wheels, so `pip install compas_cra` no longer needs conda, homebrew or a manually downloaded solver. Wheels are built for Windows, macOS (Apple Silicon and Intel) and manylinux (x86_64 and aarch64) by `.github/workflows/release.yml`.
* Added `packaging/` with the scripts that build IPOPT 3.14.19 from source with coinbrew (MUMPS linear solver, no HSL), pack the platform wheels and test them in a clean environment.
* Added `compas_cra._ipopt` to locate the solver, and a `compas-cra-ipopt` console script to check which solver will be used.
* Added the `viz` optional dependencies (`pip install compas_cra[viz]`).
* Added `packaging/check_release.py`, run by the publish job before uploading: a release is rejected unless it carries a wheel for every supported platform and each wheel actually contains an ipopt executable.

### Changed

* Declared the runtime dependencies (`compas`, `compas_assembly`, `numpy`, `pyomo`, `scipy`, `shapely`) in `requirements.txt`, which was empty.
* `cra_solve`, `cra_penalty_solve` and `rbe_solve` now build their solver with `compas_cra.equilibrium._solver.ipopt_solver`, which points pyomo at the bundled executable. Solver options and tolerances are unchanged.
* Rewrote the installation instructions around pip; the manual Windows workaround for missing ipopt binaries is no longer needed.

* Replaced kernel-layer `MatrixConstraint` with standard `pyo.Constraint` rules in `static_equilibrium_constraints` for compatibility with recent Pyomo versions.
* Pinned Python to `< 3.14` for viewer compatibility.

### Removed


## [0.5.0] 2024-10-21

### Added

### Changed

* Changed compas_view2 to compas_viewer.
* Changed sample files to COMPAS 2 format.
* Fixed bug in temp viewer arrow solution.

### Removed

* Removed `matplotlib` from env files.
* Removed `pip` requirements.
* Removed incompatible interface info.


## [0.4.0] 2024-03-02

### Added

* Add delete block and blocks methods in CRA_Assembly class. 
* A script to export mesh to json in Rhino. 

### Changed

### Removed


## [0.3.0] 2022-11-06

### Added

### Changed

### Removed


## [0.2.2] 2022-09-29

### Added

* Add example folder directory to tutorial docs for easy access. 

### Changed

* Fix some typos and wrong url links. 
* Change ipopt installation guide using conda.

### Removed


## [0.2.1] 2022-09-02

### Added

### Changed

### Removed


## [0.2.0] 2022-09-02

### Added

### Changed

### Removed

