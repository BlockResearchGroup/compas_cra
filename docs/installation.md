# Installation

## Stable

```bash
pip install compas_cra
```

That is all that is needed, on Windows, macOS (Apple Silicon and Intel) and Linux.
The [IPOPT](https://coin-or.github.io/Ipopt/) solver is compiled into the package
itself as an extension module, so solving happens in-process: no conda environment, no
homebrew and no solver executables are involved. The wheels cover CPython 3.9 to 3.13.

To also install the viewers:

```bash
pip install compas_cra[viz]
```

Verify the solver is available with:

```bash
python -c "from compas_cra import _native; print(_native.IPOPT_VERSION)"
```

## Rhino 8

Start a Python 3 script in the ScriptEditor with this header and run it — Rhino
installs everything on the first run:

```python
#! python3
# venv: compas-cra
# r: compas_cra
```

Ready-to-run examples are in the repository under `scripts/`.

## Development

Working on `compas_cra` itself means building the solver, because the solver is part
of the package: `pip install .` and `pip install -e .` both compile the
`compas_cra._native` extension, and that links against a static IPOPT tree produced
by `packaging/build_ipopt.sh`. Without such a tree the install stops with
`No IPOPT build at IPOPT_PREFIX=...`.

### Prerequisites

`invoke setup` installs these wherever that is possible without root — MSYS2 packages
on Windows, Homebrew formulae on macOS — so on those two platforms this section is
background rather than something to work through. On Linux the packages need root, so
`invoke setup` prints the line below and stops rather than running it.

The IPOPT build needs a Fortran compiler (its MUMPS linear solver is Fortran), a static
BLAS/LAPACK, a C/C++ toolchain, and `git`, `curl`, `make`, `patch` and
`pkg-config` for coinbrew. CMake and Ninja are *not* needed up front — the build
backend brings its own.

Debian / Ubuntu:

```bash
sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config
```

Fedora / RHEL / AlmaLinux:

```bash
sudo dnf install gcc gcc-c++ gcc-gfortran openblas-devel openblas-static \
    glibc-static libstdc++-static git curl make patch pkgconf
```

macOS:

```bash
brew install gcc bash
```

`gfortran` comes with the `gcc` formula, and BLAS/LAPACK is the system Accelerate
framework, so nothing else is needed. The `bash` formula is: coinbrew refuses to run
under bash 3, which is still what macOS ships.

Windows: the IPOPT build runs in an [MSYS2](https://www.msys2.org) UCRT64 shell.

```bash
pacman -S git patch make diffutils pkgconf \
    mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran \
    mingw-w64-ucrt-x86_64-openblas
```

### Building

[uv](https://docs.astral.sh/uv) manages the environment.

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/Scripts/activate   # .venv/bin/activate on macOS and Linux

uv pip install invoke compas_invocations2
invoke setup
invoke test
```

`invoke setup` is the whole build, on Windows, macOS and Linux alike: it installs the
toolchain where that can be done without root, stages IPOPT if it is not staged
already, and installs the package editable with the compilers and link directories the
extension needs — on Windows that means running `build_ipopt.sh` inside the MSYS2
UCRT64 shell and pointing CMake at its `gcc`, `g++` and `gfortran`, which is the part
that is tedious to get right by hand.

It takes about fifteen minutes the first time, almost all of it IPOPT. It is
idempotent: run it again and it reinstalls the package and leaves the staged IPOPT tree
alone. `invoke setup --jobs=N` raises the IPOPT build parallelism, but see the warning
in `--help` first — MUMPS races and dies under a parallel build.

It installs with `[dev]`, which pulls in `[docs]` as well, so `invoke docs` works
straight afterwards without a second command. `.venv` is git-ignored.

Doing it by hand instead: run `packaging/build_ipopt.sh` from the shell your platform
needs, then `pip install -e ".[dev]"` with `IPOPT_PREFIX` pointed at the resulting
`build/ipopt/stage`.

`build_ipopt.sh` takes about fifteen minutes and is a one-off: it stages IPOPT into
`build/ipopt/stage`, which every later install reuses. Rebuild it only when the script
or the IPOPT version changes. Set `IPOPT_PREFIX` to use a stage tree from elsewhere,
`IPOPT_EXTRA_LINK` and `EXTRA_LINK_DIRS` for extra link flags and directories.

There is no packaged IPOPT that can stand in for that step. conda-forge, for instance,
ships IPOPT as shared libraries resolved against a separate `mumps-seq` package, while
`CMakeLists.txt` links `libipopt` and `libcoinmumps` as *static* archives — that is what
makes the extension self-contained, and it is what the stage tree provides.

### Tasks

`tasks.py` is the entry point for everything else, via
[invoke](https://www.pyinvoke.org). The same tasks run in CI, so a green run locally
means a green run there.

| Command | What it does |
| --- | --- |
| `invoke setup` | Build the solver and install the package, on any platform |
| `invoke docs` | Build the site into `dist/docs`, failing on a broken reference |
| `invoke docs-serve` | Serve the site locally with live reload |
| `invoke lint` | `ruff check --fix src tests` |
| `invoke format` | `ruff format src tests` |
| `invoke test` | Run the test suite |
| `invoke release <major\|minor\|patch>` | Bump, tag and push, which triggers the release |

Only `invoke test` and `invoke release` need the compiled solver. The documentation is
generated from the sources in `src/` rather than from an imported package, so
`invoke docs` and `invoke lint` work in an environment where `pip install -e .` has
never run, which is the quick way in when the change is to the documentation or to
Python code only:

```bash
uv venv --python 3.12
source .venv/Scripts/activate
uv pip install -r requirements-docs.txt
invoke docs
```

`invoke release` pushes a `v*` tag, and that tag is what `.github/workflows/release.yml`
triggers on; it needs push access to the repository.

### Without a local toolchain

On Linux with Docker, `cibuildwheel` builds exactly what CI builds — IPOPT and all —
inside the manylinux container, so nothing has to be installed on the host:

```bash
pip install cibuildwheel
CIBW_BUILD="cp312-*" cibuildwheel --output-dir dist
```

Everything else — the manylinux image, the `dnf` line that installs gfortran and
OpenBLAS in the container, the IPOPT build itself and the smoke test — comes from
`[tool.cibuildwheel]` in `pyproject.toml`, so this is exactly what CI runs.
`CIBW_BUILD` is only there to limit the run to one CPython; drop it to build all five.

That takes about five minutes and leaves an installable wheel in `dist/`. It is also
the quickest way to get a working solver into a local environment when you only want to
change Python code.

See `native/README.md` and `packaging/README.md` for what is being built and why.
