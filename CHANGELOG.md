# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* `invoke setup` builds the solver and installs the package in one command, on Windows, macOS and Linux. It installs the toolchain where that is possible without root (MSYS2 packages via pacman, Homebrew formulae including the `gfortran` symlink the `gcc` formula does not create; on Linux it prints the `apt`/`dnf` line, which needs root, and stops), stages IPOPT unless `build/ipopt/stage` already holds it, and installs editable with the right environment - on Windows that means running `build_ipopt.sh` inside the MSYS2 UCRT64 shell and pointing `CMAKE_ARGS` and `EXTRA_LINK_DIRS` at its toolchain, which was four exports to get right by hand. It is idempotent, and defaults to `JOBS=1` because MUMPS races under a parallel build.

* Local development instructions in the installation docs: the compilers and libraries `packaging/build_ipopt.sh` needs on Linux, macOS and Windows, and the Docker/cibuildwheel route for building a wheel without installing any of them.
* `cra_solve`, `cra_penalty_solve` and `rbe_solve` appear in the API reference. They were assignments (`cra_solve = cra_solve_native`) rather than imports, so the documentation generator saw untyped attributes and rendered nothing for the three names that are the package's primary entry points. They are import aliases now, which changes nothing at runtime — each name is still the same function object as its `_native` counterpart.
* `invoke docs-serve` serves the documentation locally with live reload.

### Changed

* The documentation navigation is grouped into sections that Material renders as tabs in the header bar (`navigation.tabs`), the way `compas_cgal` does it: Getting Started, Examples, API Reference and About. The nav was a flat list before, so every page sat in the sidebar at the same level. The entries under the old `Other` section had no titles and were rendering as bare file names, and the reference entries are spelled with their full module names. No page moved or was renamed.
* The CRA solvers use IPOPT's adaptive barrier update (`mu_strategy="adaptive"`). CRA's complementarity constraints are degenerate and the monotone update crawls on them: the 20-block arch needed 1964 of IPOPT's 3000 permitted iterations, leaving almost no margin, and on macOS — where the wheel links Accelerate rather than OpenBLAS, and the arithmetic differs slightly — it ran out of iterations altogether and raised `solve failed: failed (Maximum_Iterations_Exceeded)`. The adaptive update reaches the same answer in 668 iterations (interface resultants agree to 1e-10), and 307 rather than 2757 on a 40-block arch. The tolerances are deliberately left alone: relaxing them makes it worse, turning a solve that stops at a feasible point into a solve that exhausts the cap.
* `test_cra_arch` no longer skips itself when the arch stalls, and a companion test asserts the solve keeps its iteration headroom. That skip is what let the macOS failure reach users while CI stayed green.
* IPOPT is built serially (`JOBS=1`). MUMPS' makefiles are missing dependencies, so a parallel build races and dies linking `libcoinmumps` — always just after `dtype3_root.F`, and with coinbrew reporting nothing but `Build failed, see error output above` and no error above it. It is intermittent, and it failed three of six pipeline runs in compas_sandbox before being pinned down.
* The PyPI upload and the GitHub release moved from `pipeline.yml` into `release.yml`. Trusted publishing matches the `job_workflow_ref` claim, which for a job defined inside a reusable workflow names the reusable file, so uploading from `pipeline.yml` presents `pipeline.yml` to a publisher configured for `release.yml` and fails with `invalid-publisher`. This is measured, not theoretical: compas_sandbox hit exactly that on the first release that reached the upload after the same refactor.
* **What the wheels actually bundle.** IPOPT 3.14.19 is compiled from source by `packaging/build_ipopt.sh` (via coinbrew) as *static* libraries, and `libipopt` and `libcoinmumps` are linked into the `compas_cra._native._core` extension. There is no `ipopt` executable, no `.nl` files and no subprocess: a solve is a function call. MUMPS is the linear solver, which is what CRA relies on; HSL is deliberately absent because it is not redistributable. BLAS/LAPACK is the one piece that differs by platform — a static OpenBLAS on Linux and Windows, Apple's Accelerate on macOS — and that difference is exactly why the arch stalled on Mac and not elsewhere. The remaining runtimes (the Fortran runtime, OpenBLAS where it is dynamic) stay shared and are grafted into the wheel by the platform repair tools: auditwheel on Linux, delocate on macOS, delvewheel on Windows. Licences travel with the binaries: IPOPT is EPL-2.0, MUMPS and the ASL are permissive, OpenBLAS is BSD-3-Clause, and the statically linked libgfortran/libgcc are covered by the GCC Runtime Library Exception.
* One package again: the IPOPT binding is built into `compas_cra` itself as `compas_cra._native`, and `compas_cra` ships as platform wheels for CPython 3.9-3.13 on Windows, macOS (Apple Silicon and Intel) and manylinux (x86_64 and aarch64). `pip install compas_cra` installs exactly one distribution, solver included.
* The build backend is scikit-build-core, which drives the CMake project at the repository root. It reads no `requirements.txt`, so the dependencies now live in `pyproject.toml`.
* The cibuildwheel settings live in `[tool.cibuildwheel]` in `pyproject.toml` rather than in workflow environment variables, so `cibuildwheel` reproduces the CI wheel locally with no environment set. Only the runner-dependent parts (workspace paths, the macOS deployment target, the MSYS2 location) stay in `pipeline.yml`.
* One CI pipeline for every push and every release: `build.yml` (pushes and PRs to main) and `release.yml` (`v*` tags) both call the reusable `pipeline.yml`, which differs only by its `publish` input. A release therefore runs the exact chain that was already green on main — IPOPT from source, all five platforms, the smoke tests, the test suite against a built wheel and the docs — and then publishes it. `wheels.yml` and `docs.yml` are gone, folded into that pipeline.
* `packaging/check_release.py` now checks for the `compas_cra._native._core` extension rather than a vendored `ipopt` executable, which the statically linked build no longer ships, and additionally requires every supported CPython on every platform. Verified against the artifacts of the PR build: it passes the complete set and rejects a partial one.
* `invoke release` is defined in `tasks.py` instead of coming from `compas_invocations2`, without its local `python -m build` step: building a wheel here needs a staged IPOPT tree, and the published artifacts are built by CI on all five platforms regardless.
* The documentation is built with mkdocs (mkdocs-material + mkdocstrings) instead of sphinx, following `compas_model`. Every page is markdown: the API reference is five one-line `::: compas_cra.<module>` pages instead of six hand-maintained `.rst` files plus 80 checked-in `autosummary` stubs, and the changelog and licence pages include the repository's own files rather than duplicating them. The docs dependencies are a `docs` extra in `pyproject.toml`; `invoke docs` and the CI docs job both run `mkdocs build --strict`, so a broken link or a missing page now fails the build.
* Two docstrings that documented parameters the functions do not take: `CRA_Assembly.add_to_interfaces` described a `type` argument it never had, and `Arch` documented `n` instead of `num_blocks` and omitted `extra_support`.
* `packaging/build_ipopt.sh` no longer trusts the third-party source downloads. The `get.ASL` and `get.Mumps` scripts fetch their tarballs with a bare `curl -L -O`, with no integrity check and no retry, so a stalled transfer lands a short file that is then gunzipped anyway — which is how an aarch64 build failed with `gzip: MUMPS_5.8.2.tar.gz: not in gzip format` after receiving 54 kB of a 4.3 MB tarball. The script now points `CURL_HOME`/`WGETRC` at a config that abandons and retries a stalled transfer, and checks that `ThirdParty/ASL/solvers` and `ThirdParty/Mumps/MUMPS` actually exist, re-running the `get.*` script that came up empty (a second `coinbrew fetch` will not: it only runs `get.*` when the clone's revision changed).
* `invoke release` no longer deletes the staged IPOPT tree. It called `compas_invocations2.build.clean`, whose default removes `build/` - which in this repository is where `packaging/build_ipopt.sh` stages IPOPT, and what `CMakeLists.txt` defaults `IPOPT_PREFIX` to. Every release therefore threw away about fifteen minutes of build and left the next `pip install -e .` failing with `No IPOPT build at IPOPT_PREFIX=...`. It is called with `builds=False` now, so `dist/` and the egg-info still go and the stage tree stays.
* `pip install -e ".[dev]"` installs the documentation stack too, because `[dev]` now includes `[docs]`. Neither extra could run `invoke docs` on its own - `[dev]` has no mkdocs and `[docs]` has no invoke - so the command the installation instructions gave produced an environment the documented task did not run in.
* `pymdown-extensions` is declared in the `docs` extra. `mkdocs.yml` configures `pymdownx.*` directly, so it is a direct dependency; it arrived only as a transitive dependency of `mkdocs-material`.
* CI builds the documentation with `invoke docs` rather than its own copy of the `mkdocs build` command line, so the local task and the pipeline cannot drift apart.
* The development instructions start from a `uv` environment and state that `pip install -e ".[dev]"` is the only install step, list the `invoke` tasks and which of them need the compiled solver (only `test` and `release` - the documentation is generated from `src/`, not from an imported package), and note that `conda install ipopt` cannot stand in for `build_ipopt.sh` because conda-forge builds IPOPT shared against a separate `mumps-seq` while `CMakeLists.txt` links it statically.
* The documentation is deployed on every push to `main`, not only from a release tag. The deploy step was gated on `inputs.publish`, which only `release.yml` sets, so merging a documentation fix built the site and threw it away and the published site tracked the last tag rather than `main`. It costs no extra CI: `build.yml` already runs the whole pipeline on every push to `main`, so the wheel the docs are built against exists either way.
* Every URL points at `BlockResearchGroup/compas_cra`. `site_url`, `repo_url`, the `[project.urls]` block, the README badges, `CITATION.cff` and the links in the tutorial and contribution pages all named the `petrasvestartas` fork, so the published site advertised a canonical URL that is not where it is served from - GitHub Pages serves it at `blockresearchgroup.github.io/compas_cra`.
* `requirements-docs.txt` installs the documentation and linting toolchain in one command, without building the package. An extra cannot do this: installing `compas_cra[docs]` builds `compas_cra`, which compiles the solver. The documentation never needs it.
* Two links in the README that were broken independently of the rename: the banner image pointed into `docs/_images`, which has not existed since the mkdocs migration moved it to `docs/assets/images`, and the examples link pointed at a `latest/examples.html` path on github.com rather than at the documentation site.

### Removed

* Removed the separate `compas_cra_native` distribution and its `native/pyproject.toml`; the binding source stays in `native/`, built by the root `CMakeLists.txt`. An install from source (`pip install .`, editable included) now compiles the extension and needs `packaging/build_ipopt.sh` to have run first.
* Removed `requirements.txt`, `requirements-dev.txt` and `requirements-viz.txt`; `pip install -e ".[dev]"` replaces `pip install -r requirements-dev.txt`.
* Removed the `build_ghuser_components` task and its `ghuser` configuration from `tasks.py`; the paths pointed at `src/compas_notebook`, and compas_cra ships no Grasshopper components.
* Removed `docs/conf.py`, every `.rst` page, the generated `docs/api/` tree and the `sphinx_compas2_theme` dependency.

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

## [0.6.0] 2026-08-26

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

