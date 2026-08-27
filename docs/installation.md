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

The solver is part of the package: `pip install -e .` compiles the `compas_cra._native`
extension against a static IPOPT tree. `invoke setup` does everything — toolchain,
IPOPT, editable install — on Windows, macOS and Linux.

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/Scripts/activate   # .venv/bin/activate on macOS and Linux

uv pip install invoke compas_invocations2
invoke setup                    # ~15 minutes the first time, seconds after that
invoke test
```

On Linux the toolchain needs root, so `invoke setup` prints the `apt`/`dnf` line and
stops — run it, then re-run `invoke setup`.

After that, on all three platforms:

- **Python edits need nothing** — the editable install serves `src/` directly.
- **C++ edits, or any reinstall**: `uv pip install -e .` in the activated venv, plain.

Nothing global is configured. The venv's own activation scripts carry the toolchain
`PATH` (written once by `invoke setup`, Windows only), `CMakeLists.txt` finds the
compilers itself, the staged IPOPT lives in the repository under `build/ipopt/stage`,
and the runtime libraries are installed next to the extension — importing `compas_cra`
needs no environment at all.

## Documentation only

No solver build needed — the docs are generated from `src/`, not from an imported
package.

```bash
uv venv --python 3.12
source .venv/Scripts/activate
uv pip install -r requirements.txt

invoke docs         # build the site into dist/docs — open dist/docs/index.html
invoke docs-serve   # or serve it live at http://localhost:8000
```

## Tasks

| Command | What it does |
| --- | --- |
| `invoke setup` | Build the solver and install the package, on any platform |
| `invoke docs` | Build the site into `dist/docs`, failing on a broken reference |
| `invoke docs-serve` | Serve the site locally with live reload |
| `invoke lint` | `ruff check --fix src tests` |
| `invoke format` | `ruff format src tests` |
| `invoke test` | Run the test suite |
| `invoke release <major\|minor\|patch>` | Bump, tag and push, which triggers the release |

Only `invoke test` and `invoke release` need the compiled solver. `invoke release`
pushes a `v*` tag — that is what triggers `.github/workflows/release.yml` — and needs
push access to the repository.

## Manual build

What `invoke setup` runs, step by step.

Prerequisites — a Fortran compiler, a static BLAS/LAPACK, a C/C++ toolchain, and
`git`, `curl`, `make`, `patch`, `pkg-config` for coinbrew:

```bash
# Debian / Ubuntu
sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config

# Fedora / RHEL / AlmaLinux
sudo dnf install gcc gcc-c++ gcc-gfortran openblas-devel openblas-static     glibc-static libstdc++-static git curl make patch pkgconf

# macOS (gfortran comes with gcc; bash because coinbrew refuses to run under bash 3)
brew install gcc bash

# Windows, in an MSYS2 UCRT64 shell (https://www.msys2.org)
pacman -S git patch make diffutils pkgconf     mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran     mingw-w64-ucrt-x86_64-openblas
```

Build IPOPT, then install:

```bash
JOBS=1 packaging/build_ipopt.sh   # stages into build/ipopt/stage; ~15 minutes, once
pip install -e ".[dev]"
```

Keep `JOBS=1`: MUMPS races under a parallel build. On Windows run `build_ipopt.sh`
from the MSYS2 UCRT64 shell; the install itself needs no compiler flags —
`CMakeLists.txt` finds the MinGW toolchain on its own. Set `IPOPT_PREFIX` to
use a stage tree from elsewhere, `IPOPT_EXTRA_LINK` and `EXTRA_LINK_DIRS` for extra
link flags and directories. Packaged IPOPTs (conda-forge etc.) do not work here: they
are shared builds, and `CMakeLists.txt` links `libipopt` and `libcoinmumps` statically —
that is what makes the extension self-contained.

## Without a local toolchain

On Linux with Docker, `cibuildwheel` builds exactly what CI builds — IPOPT and all —
inside the manylinux container:

```bash
pip install cibuildwheel
CIBW_BUILD="cp312-*" cibuildwheel --output-dir dist
```

About five minutes, and everything — the image, the container `dnf` line, the IPOPT
build, the smoke test — comes from `[tool.cibuildwheel]` in `pyproject.toml`, so this
is exactly what CI runs. Drop `CIBW_BUILD` to build all five CPythons.

See `native/README.md` and `packaging/README.md` for what is being built and why.
