# Notes for working on compas_cra

Things that cost real time to discover and are not visible from the code alone. The
changelog says what changed; this says what will bite you.

## The solver

**`mu_strategy="adaptive"` in `_CRA_OPTIONS` is load bearing.** CRA's complementarity
constraints are degenerate and IPOPT's monotone barrier update crawls on them. Without
it the 20-block arch takes 1964 of its 3000 permitted iterations and reaches its answer
only through the restoration phase — enough margin on some platforms and not on others,
which is how macOS users got `solve failed: failed (Maximum_Iterations_Exceeded)` while
CI stayed green. With it: 668 iterations, same answer to 1e-10.

**Do not "fix" the tolerances.** `tol=1e-10 / constr_viol_tol=1e-12` look absurd for
double precision and they are, but they are what drives the solve into restoration,
which is where it finds its answer. Measured: `tol=1e-8 / constr_viol_tol=1e-9` turns a
solve that stops at a feasible point into one that exhausts the 3000-iteration cap.
Looser is worse, not safer.

**Never let the arch test skip itself again.** `test_cra_arch` used to `pytest.skip` on
a stall. That is exactly why the macOS failure reached users while CI was green.
`test_cra_arch_has_iteration_headroom` guards the margin, not just the answer — a solve
creeping back toward the cap is the actual regression.

## The IPOPT build

**MUMPS races under parallel make.** Its makefiles are missing dependencies, so a `-j`
build intermittently dies linking `libcoinmumps` — always just after `dtype3_root.F`,
and coinbrew reports only `Build failed, see error output above` with no error above it.
Hence `JOBS=1`: on the `before-all` command in `pyproject.toml` for Linux (the container
does not inherit the workflow env) and in `pipeline.yml`'s `env:` for macOS and Windows.

**The caches hide a broken from-scratch build.** `hashFiles('packaging/build_ipopt.sh')`
is part of the ipopt cache key, so *any* edit to that file evicts the macOS and Windows
stage trees and makes them compile for the first time in months. That is when the race
above surfaces. Keep `JOBS=1` on those platforms even though it looks inert today.

**`hashFiles()` is not a plain sha256.** It hashes the file, then hashes the digest —
`sha256(sha256(content))`. Comparing a plain `sha256sum` against a cache key will tell
you the cache misses when it actually hits.

**macOS links Accelerate, everything else OpenBLAS.** That difference is not cosmetic:
it changes IPOPT's iterate path, which is why the arch stalled on Mac and nowhere else.
Moving macOS onto OpenBLAS for bit-reproducibility was tried in compas_sandbox and
reverted — libtool will not link a static `libopenblas.a` into the shared
`libcoinmumps` coinbrew builds regardless of `--disable-shared`, and the `.dylib` fails
silently. Linux only escapes it because `-lcralapack` is a linker script, so libtool
sees an ordinary `-l` flag rather than a library path.

## Releasing

**The PyPI upload must live in `release.yml`, not in the reusable `pipeline.yml`.**
Trusted publishing matches the `job_workflow_ref` claim, which names the file a job is
*defined* in. A job in the reusable workflow presents `pipeline.yml` to a publisher
configured for `release.yml` and fails with
`invalid-publisher: valid token, but no corresponding publisher`. Reusable workflows are
not supported by trusted publishing at all. This was measured on compas_sandbox, whose
first release to reach the upload after the same refactor failed exactly this way — an
earlier note here claimed the opposite, reasoning from provenance that predated that
refactor.

On PyPI, compas_cra needs a trusted publisher for this repository with `release.yml` as
the workflow and `pypi` as the environment.

**A flaky wheel job does not need a new tag.** `gh run rerun <id> --failed` re-runs just
the failed job on the existing tag; the other jobs keep their results and `publish`
cascades if it passes.

## Local development

A source checkout has no compiled `compas_cra._native._core`, so every solver test
skips (`pytest.importorskip`). To run them, `packaging/build_ipopt.sh` then
`pip install .`. Building a wheel always needs a staged IPOPT tree first, which is why
`invoke release` has no local `python -m build` step.
