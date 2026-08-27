# Installation

## Install

```bash
pip install compas_cra
```

Windows, macOS and Linux, CPython 3.9-3.13. The [IPOPT](https://coin-or.github.io/Ipopt/)
solver is compiled into the package — no conda, no homebrew, no solver executables.

```bash
pip install compas_cra[viz]                                              # with viewers
python -c "from compas_cra import _native; print(_native.IPOPT_VERSION)" # verify
```

In Rhino 8, run a script with this header — Rhino installs everything on the first run
(examples in `scripts/`):

```python
#! python3
# venv: compas-cra
# r: compas_cra
```

## Development

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/Scripts/activate   # .venv/bin/activate on macOS and Linux

uv pip install invoke compas_invocations2
invoke setup                    # toolchain + IPOPT (~15 minutes, once) + editable install
invoke test
```

On Linux the toolchain needs root: `invoke setup` prints the `apt`/`dnf` line — run it,
re-run `invoke setup`. After that, on all platforms:

- **Python edits need nothing** — the editable install serves `src/` directly.
- **C++ edits, or any reinstall**: `uv pip install -e .`, plain.

Nothing global is configured: the venv's activation scripts carry the toolchain `PATH`
(Windows), `CMakeLists.txt` finds the compilers, IPOPT lives in `build/ipopt/stage`,
and the runtime libraries are installed next to the extension.

### Documentation only

No solver build needed — the docs are generated from `src/`:

```bash
uv venv --python 3.12
source .venv/Scripts/activate
uv pip install -r requirements.txt
invoke docs
```

### Tasks

| Command | What it does |
| --- | --- |
| `invoke setup` | Build the solver and install the package, on any platform |
| `invoke docs` | Serve the site at `localhost:8000` with live reload; `--no-serve` builds `dist/docs` |
| `invoke lint` | `ruff check --fix src tests` |
| `invoke format` | `ruff format src tests` |
| `invoke test` | Run the test suite |
| `invoke release <major\|minor\|patch>` | Bump, tag and push, which triggers the release |

Only `test` and `release` need the compiled solver; `release` needs push access to the
repository.

### Manual build

What `invoke setup` runs. Prerequisites — Fortran + C/C++ compilers, static
BLAS/LAPACK, and `git curl make patch pkg-config` for coinbrew:

```bash
# Debian / Ubuntu
sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config

# Fedora / RHEL / AlmaLinux
sudo dnf install gcc gcc-c++ gcc-gfortran openblas-devel openblas-static     glibc-static libstdc++-static git curl make patch pkgconf

# macOS
brew install gcc bash

# Windows, in an MSYS2 UCRT64 shell (https://www.msys2.org)
pacman -S git patch make diffutils pkgconf     mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran     mingw-w64-ucrt-x86_64-openblas
```

```bash
JOBS=1 packaging/build_ipopt.sh   # stages build/ipopt/stage; keep JOBS=1 - MUMPS races
pip install -e ".[dev]"
```

`IPOPT_PREFIX` points at a stage tree elsewhere; `EXTRA_LINK_DIRS` and
`IPOPT_EXTRA_LINK` add link directories and flags. Packaged IPOPTs (conda-forge etc.)
do not work: they are shared builds, and the extension links `libipopt` and
`libcoinmumps` statically.

### Without a local toolchain

On Linux with Docker, `cibuildwheel` builds exactly what CI builds, IPOPT included:

```bash
pip install cibuildwheel
CIBW_BUILD="cp312-*" cibuildwheel --output-dir dist
```

Everything comes from `[tool.cibuildwheel]` in `pyproject.toml`; drop `CIBW_BUILD` to
build all five CPythons. See `native/README.md` and `packaging/README.md` for what is
built and why.
