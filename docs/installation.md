# Installation

## Install

```bash
pip install compas_cra
```

Windows, macOS and Linux, CPython 3.9-3.13. The [IPOPT](https://coin-or.github.io/Ipopt/)
solver is compiled into the package.

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

### Windows

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/Scripts/activate

uv pip install invoke compas_invocations2
invoke setup
invoke test
```

### macOS

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/bin/activate

uv pip install invoke compas_invocations2
invoke setup
invoke test
```

### Linux

```bash
sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config

git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra

uv venv --python 3.12
source .venv/bin/activate

uv pip install invoke compas_invocations2
invoke setup
invoke test
```

### After the first build

```bash
uv pip install -e .   # rebuild after C++ changes; Python changes need no reinstall
```

### Documentation only

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

### Without a local toolchain

```bash
pip install cibuildwheel
CIBW_BUILD="cp312-*" cibuildwheel --output-dir dist   # Linux/Docker; exactly what CI runs
```

For building IPOPT by hand, see `packaging/README.md`.
