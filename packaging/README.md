# Packaging

`compas_cra` solves its models with [IPOPT](https://coin-or.github.io/Ipopt/)
compiled into the package itself as the `compas_cra._native` extension module (the
binding lives in [`native/`](../native)). `build_ipopt.sh` builds IPOPT (with the MUMPS
linear solver) from source with coinbrew as **static libraries**, staged into
`build/ipopt/stage`; the extension links against that stage tree.
`.github/workflows/pipeline.yml` runs the build per platform and packs the wheels for
CPython 3.9–3.13 with cibuildwheel.

| script | what it does |
| --- | --- |
| `build_ipopt.sh` | builds IPOPT + MUMPS from source and stages libraries and headers into `build/ipopt/stage` |

## What goes into the stage tree

- **IPOPT 3.14.19**, Eclipse Public License 2.0, as `libipopt.a`
- **MUMPS** as the linear solver (`libcoinmumps.a`). This is IPOPT's default and the
  one `compas_cra` relies on, since none of the formulations set `linear_solver`.
- **BLAS/LAPACK**: a static OpenBLAS on Linux and Windows, Accelerate on macOS — named
  by its SDK stub path, because `-framework Accelerate` does not survive libtool.
- **No HSL.** The HSL linear solvers (`ma27`, `ma57`, ...) are not redistributable.

IPOPT and MUMPS are linked statically into the extension module. The remaining
runtimes (the Fortran runtime, OpenBLAS) are linked dynamically and grafted into the
wheels by the platform repair tools (`auditwheel` / `delocate` / `delvewheel`) —
the standard approach of the scientific Python wheels. The result is that `compas_cra`
ships as platform wheels rather than one universal wheel: one per CPython version and
platform, each carrying its own solver.

## Building locally

```bash
packaging/build_ipopt.sh    # ~15 minutes, stages into build/ipopt/stage
pip install .               # builds the extension against the stage tree
python native/smoke.py
```

On Linux, build inside the manylinux container so the result runs on older
distributions than your own (this is what CI does through cibuildwheel):

```bash
docker run --rm -v "$PWD":/io -w /io quay.io/pypa/manylinux_2_28_x86_64 bash -c '
  dnf install -y --enablerepo=powertools glibc-static libstdc++-static openblas-static openblas-devel
  CFLAGS="-O2 -fPIC" CXXFLAGS="-O2 -fPIC" FFLAGS="-O2 -fPIC" FCFLAGS="-O2 -fPIC" packaging/build_ipopt.sh'
```

On Windows the build runs in an MSYS2 UCRT64 shell; on macOS it needs `gfortran`
(`brew install gcc`).

## Releasing

`.github/workflows/pipeline.yml` is the one build pipeline: it compiles IPOPT, builds
and smoke-tests a `compas_cra` wheel for every platform and CPython, runs the test
suite against one of them, and builds the docs — leaving everything as run artifacts.
`build.yml` runs it on every push. `release.yml` runs it on a `v*` tag and then
publishes those artifacts to PyPI with trusted publishing.

`invoke release <patch|minor|major>` is the way to cut one: it runs the tests, bumps
the version, tags it and pushes. It uploads nothing itself — the tag is what triggers
the publish, and the wheels are built on the runners, never locally.
