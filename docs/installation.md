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

Run `packaging/build_ipopt.sh` from that shell, then `pip` from your normal Python
with `IPOPT_PREFIX` pointed at the resulting `build/ipopt/stage`.

### Building

```bash
git clone https://github.com/petrasvestartas/compas_cra.git
cd compas_cra
packaging/build_ipopt.sh
pip install -e ".[dev]"
invoke test
```

`build_ipopt.sh` takes about fifteen minutes and is a one-off: it stages IPOPT into
`build/ipopt/stage`, which every later install reuses. Rebuild it only when the script
or the IPOPT version changes. Set `IPOPT_PREFIX` to use a stage tree from elsewhere,
`IPOPT_EXTRA_LINK` and `EXTRA_LINK_DIRS` for extra link flags and directories.

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
