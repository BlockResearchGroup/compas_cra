# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.0] 2026-08-29

### Added

* Platform-logic tests for the invoke tasks (`tests/test_tasks.py`): the macOS and Linux environment construction of `invoke setup` is pure logic, pinned per-OS with mocks so a refactor cannot silently break a platform nobody is sitting at. Skipped where the dev tooling is not installed.

* Rhino 8 ports of the three curved-interface examples (`scripts/rhino_cra_curve_3_blocks.py`, `rhino_cra_cube_curve_short.py`, `rhino_cra_cube_curve_tall.py`) with a shared drawing module (`scripts/rhino_cra_view.py`) that mirrors the desktop visualization - the same mesh arrows (cylinder shaft, conical head), colors and element grouping, on toggleable Rhino layers, drawn through the compas 2 Scene API. Open in the ScriptEditor and run; the `# r:` header installs compas_cra on first use.

* A viewer smoke test (`tests/test_viewers.py`): the full `cra_view` scene is constructed headless against the installed `compas_viewer`, with only the blocking `show()` stubbed, so an API drift in the viewer fails in CI instead of on a user's screen. Skipped when the viz extra is absent.
* A regression test for the penalty solver on the arch, both variants of the `06_arch_penalty` example, pinning the penalty solution to `cra_solve` on the standard arch.

* `invoke setup` builds the solver and installs the package in one command, on Windows, macOS and Linux. It installs the toolchain where that is possible without root (MSYS2 packages via pacman, Homebrew formulae including the `gfortran` symlink the `gcc` formula does not create; on Linux it prints the `apt`/`dnf` line, which needs root, and stops), stages IPOPT unless `build/ipopt/stage` already holds it, and installs editable with the right environment - on Windows that means running `build_ipopt.sh` inside the MSYS2 UCRT64 shell and pointing `CMAKE_ARGS` and `EXTRA_LINK_DIRS` at its toolchain, which was four exports to get right by hand. It is idempotent, and defaults to `JOBS=1` because MUMPS races under a parallel build.

* Local development instructions in the installation docs: the compilers and libraries `packaging/build_ipopt.sh` needs on Linux, macOS and Windows, and the Docker/cibuildwheel route for building a wheel without installing any of them.
* `cra_solve`, `cra_penalty_solve` and `rbe_solve` appear in the API reference. They were assignments (`cra_solve = cra_solve_native`) rather than imports, so the documentation generator saw untyped attributes and rendered nothing for the three names that are the package's primary entry points. They are import aliases now, which changes nothing at runtime — each name is still the same function object as its `_native` counterpart.
* `invoke docs` serves the documentation at `localhost:8000` with live reload - serving is the default because mkdocs links pages by directory, which only a web server resolves, so a site opened from disk shows folder listings on every click. `invoke docs --no-serve` builds into `dist/docs` instead, which is what CI deploys.

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
* `requirements.txt` installs the documentation and linting toolchain in one command, without building the package. An extra cannot do this: installing `compas_cra[docs]` builds `compas_cra`, which compiles the solver. The documentation never needs it.
* Two links in the README that were broken independently of the rename: the banner image pointed into `docs/_images`, which has not existed since the mkdocs migration moved it to `docs/assets/images`, and the examples link pointed at a `latest/examples.html` path on github.com rather than at the documentation site.
* An editable install is self-contained on Windows, the way a repaired wheel is. The extension links MSYS2's shared runtimes (libgfortran, libopenblas, libstdc++), and Python does not search PATH for extension DLLs, so importing `compas_cra._native` died with `DLL load failed` anywhere delvewheel had not vendored them - which is every editable install. The CMake install now copies the runtime DLL closure next to the module when `BUNDLE_RUNTIME_DLLS=ON`, which `invoke setup` passes; the package code does not know MSYS2 exists, and the CI wheels keep delvewheel's repair (the option defaults to OFF, so they never carry a second copy).
* `pytest` is part of the `dev` extra. `invoke test` was unusable from the documented install - CI happened to install pytest separately.
* After `invoke setup`, a plain `uv pip install -e .` rebuilds the extension with no special environment, and nothing is configured globally. `CMakeLists.txt` finds the MinGW toolchain itself on Windows (and refuses MSVC with an explanation - the GCC-built static archives cannot be linked by it, and CMake's Visual Studio default produced a LNK1181 ten minutes in), the venv's activation scripts carry the toolchain PATH and the Ninja generator (written once by `invoke setup`, venv-scoped), and the runtime DLL bundling defaults to on for Windows builds - the CI wheels pass `-DBUNDLE_RUNTIME_DLLS=OFF` and keep delvewheel's repair. Python edits need no reinstall at all: the editable install serves `src/` directly.
* The viewers work with compas_viewer 2.x again. `cra_view` passed plain Python lists to `viewer.scene.add`, which raised `SceneObjectNotRegisteredError: <class 'list'>` on the first example anyone ran with a viewer; every list is wrapped in `Collection` now, and two stale `size=`/`color=` kwargs became `pointsize=`/`pointcolor=`.
* Force and weight arrows are meshes again - cylinder shaft, conical head, widths scaling with length - restoring the compas_view2 `Arrow` shape and its `head_portion=0.2, head_width=0.07, body_width=0.02` parameters, verified against the documentation screenshots. The port had replaced them with compas_viewer's `VectorObject`, which draws a fixed pixel-width line with a four-sided pyramid head that degenerates to a sliver for axis-parallel directions - the weight arrows, drawn at the old world-unit width as a pixel width, rendered as disembodied heads. Arrows are added at opacity 0.999 on purpose: the transparency pass draws after the semi-transparent block faces, so an arrow inside a block stays bold as it did in compas_view2.
* The uneven interface resultants on the curved examples are the solver's true optimum, not an artifact: the fully converged solution (`tol=1e-8`, status optimal, 356 iterations) matches the fast acceptable-level solution to 1.7e-6 across all 72 sub-interface resultants, and their sum equals the supported weight. The visual difference against the old screenshots was the absolute-width arrow heads exaggerating small forces.
* The CRA solver's acceptable-level tolerances relaxed from 1e-8 to 1e-6. The wedge and curved-blocks examples converge to ~1e-8 violation and then cannot reach the hard tolerances - at 1e-8 acceptance they died a hair above the line (`Restoration_Failed` at 4.6e-8, `Maximum_Iterations_Exceeded` hovering at 1.6e-8). The hard tolerances are unchanged; relaxing acceptance only adds a stopping opportunity.
* The penalty solver uses the adaptive barrier update too - same degenerate complementarity constraints as the CRA solver, same crawl under the monotone update.
* The `06_arch_penalty` example passes `d_bnd=1e-2`: with the default 1e-3 the solve exhausts even a 9000-iteration cap under both barrier strategies, at 1e-2 it converges in ~65 iterations, and on the standard arch the penalty resultants then agree with `cra_solve` to 1e-3 - same physics, feasible bound. All 18 example scripts now run to completion, verified headless.
* Resultant force arrows are red only for interfaces in net tension (`sum_n < 0`) - the rule the published screenshots were made with (commit `15e5edc`, 2022-08-17). The per-vertex `<= -1e-5` check added later turned whole interfaces red over complementarity-relaxation noise (vertex normals reach -2e-3 on shelf interfaces that are in overall compression), which is why the shelf example showed red arrows where its screenshot shows green. Genuine tension still shows: net-tension interfaces render red, and per-vertex tension renders as red nodal arrows, unchanged.
* Block and support edges draw at the same width (1.5), block edges in dark gray, support edges in the support salmon - matching the reference screenshots instead of the 4-wide support outlines.
* `09_bridge` runs `cra_solve` again - the example had been switched to `cra_penalty_solve`, producing a visibly different solution from its screenshot. Verified against the 2022 code in an era-exact environment: identical resultants to five decimals, no tension anywhere. The view flags also now match the screenshot (`forcesdirect=True`): the image was made with flags that were never committed, and the image is the contract.
* The solver restores the published force distributions. The CRA system is square, and until IPOPT 3.14.11 square problems were special-cased: the convergence check ignored dual feasibility and bound complementarity, so IPOPT stopped at the first primal-feasible point of the barrier path - an interior point where every contact face carries force. Every figure in the paper (Kao et al. 2022, doi:10.1016/j.cad.2022.103216) and the docs is such a point. IPOPT 3.14.12 removed the special case, and the same problems then converge to a degenerate vertex - on the curved-interface examples the whole load lands on 8 of 72 sub-faces and the rest go to zero, which is what this package has shipped since the bundled IPOPT became newer than that. Root-caused by building IPOPT 3.14.9, 3.14.11 and 3.14.14 from source against the same MUMPS 5.9: 3.14.11 still spreads (identically to 3.14.9), 3.14.14 concentrates - the flip is exactly at 3.14.12, the release whose changelog documents the square-problem change; not MUMPS, and not this repository - every code generation from the screenshot-day commit to today produces both solutions depending only on the binary. The fix sets `mu_target=1e-5` (with tolerances to match and a raised iteration cap): the barrier stops at the interior point deliberately instead of by accident of an old binary. Calibrated against the published figures themselves rather than a rebuilt era binary - the old square-problem shortcut stopped wherever its check happened to trip, so a from-source 3.14.11 stops at a lower barrier level (force staircase of 4.5x) than the author's 2022 conda/macOS build, whose figures show near-uniform arrows. 1e-5 reproduces the figures: cube-curve-tall max/median 1.11, cube-curve-short 1.20, all 72 faces loaded. Equilibrium stays exact at 1e-8, and the arch reaches the same resultants (1.9571/0.8437) in 270 iterations. Pinning the old IPOPT was evaluated and rejected: correctly built it links and runs fine, but its stopping level is platform luck - the very thing that made the published results irreproducible in the first place.
* `13_curve-3-blocks` runs its original published parameters again (density=0.1, scale=5). They only solve in the barrier-interior regime; under strict optimization they exhaust the iteration cap on every IPOPT generation, which is why the author had switched to density=1 as a workaround days after the screenshot.
* `invoke setup` carries `MACOSX_DEPLOYMENT_TARGET=11.0` through the IPOPT build and the editable install on macOS - the same value CI pins, for Apple Silicon and Intel alike - so the static libraries and the extension agree on a target instead of failing to link against the Python build's older default. It also checks for the Xcode Command Line Tools up front (build_ipopt.sh needs xcrun for the Accelerate stub) rather than failing fifteen minutes into the build.
* The venv activation patch re-prepends the venv's own Scripts directory after the MSYS2 toolchain path: ucrt64in carries a python.exe of its own, and it was shadowing the venv's - a bare `python` in an activated shell ran the wrong interpreter.

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

