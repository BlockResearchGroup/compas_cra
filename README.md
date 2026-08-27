# COMPAS CRA

![build](https://github.com/BlockResearchGroup/compas_cra/workflows/build/badge.svg)
[![GitHub - License](https://img.shields.io/github/license/BlockResearchGroup/compas_cra.svg)](./LICENSE)
[![pip downloads](https://img.shields.io/pypi/dm/compas_cra)](https://pypi.python.org/project/compas_cra)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/compas_cra.svg)](https://pypi.python.org/project/compas_cra)
[![PyPI - Latest Release](https://img.shields.io/pypi/v/compas_cra.svg)](https://pypi.python.org/project/compas_cra)
[![DOI](https://zenodo.org/badge/374677757.svg)](https://zenodo.org/badge/latestdoi/374677757)

[Coupled Rigid-Block Analysis (CRA)](https://doi.org/10.1016/j.cad.2022.103216) implementation using [COMPAS](https://compas.dev/) framework.

> developed with <span style="color: #e25555;">&#9829;</span> by [Gene Ting-Chun Kao](https://geneatcg.com)

## Installation

```bash
pip install compas_cra
```

Nothing else is needed on Windows, macOS (Apple Silicon and Intel) or Linux. The
[IPOPT](https://coin-or.github.io/Ipopt/) 3.14.19 solver (Eclipse Public License 2.0,
MUMPS linear solver, no HSL) is compiled into the package itself as an extension
module built with [nanobind](https://github.com/wjakob/nanobind), so solving happens
in-process — one `pip install`, no separate solver package, no solver executables, no
subprocess, no conda. See [`native/README.md`](./native/README.md) for how the
extension is built and [`packaging/README.md`](./packaging/README.md) for the IPOPT
build itself.

### Local development

The solver is compiled into the package, so an install from source builds it, and that
needs IPOPT staged first:

```bash
git clone https://github.com/BlockResearchGroup/compas_cra.git
cd compas_cra
packaging/build_ipopt.sh          # ~15 minutes, once
pip install -e ".[dev]"
invoke test
```

`build_ipopt.sh` needs a Fortran compiler and a static BLAS/LAPACK — on Debian/Ubuntu
`sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config`,
on macOS `brew install gcc bash`, on Windows an [MSYS2](https://www.msys2.org) UCRT64
shell. Without a local toolchain, `CIBW_BUILD="cp312-*" cibuildwheel` builds the wheel
inside the manylinux container instead — the build is configured in `pyproject.toml`,
so that one command reproduces CI. Full instructions, per platform, are in
[the installation docs](./docs/installation.rst).

### Rhino 8

The easiest way is to let the ScriptEditor install everything for you. Start a new
Python 3 script with this header and run it — on the first run Rhino installs
`compas_cra` and all of its dependencies into its own environment (this takes a
couple of minutes, watch the progress in the ScriptEditor console):

```python
#! python3
# venv: compas-cra
# r: compas_cra
```

No extra `# r:` lines are needed for `compas`, `numpy` or `shapely` — they are
dependencies of `compas_cra` and pip installs them automatically, and the solver is
inside the `compas_cra` wheel itself.

Ready-to-run examples for the ScriptEditor are in [`scripts/`](./scripts):

- [`scripts/rhino_cra_cubes.py`](./scripts/rhino_cra_cubes.py) — three stacked cubes,
  baked as wireframe blocks with interface outlines and resultant contact forces.
- [`scripts/rhino_cra_arch.py`](./scripts/rhino_cra_arch.py) — a parametric masonry
  arch (span, rise, thickness, number of voussoirs, friction) solved and baked the
  same way.

Alternatively, install manually into Rhino's Python with the shell that ships with
Rhino (`Tools > Options > Plugins > Rhino Code > Open Shell`, or
`~/.rhinocode/py39-rh8/shell/open-shell` on macOS):

```bash
pip install compas_cra
```

Note that Rhino 8 ships CPython 3.9, which is fully supported: solver wheels exist for
CPython 3.9–3.13 on all platforms. Since the solver is a regular Python extension
module rather than an executable, it is also not affected by the antivirus policies
that quarantine unknown `.exe` files on managed Windows machines. To verify the solver
in any environment, run [`scripts/rhino_ipopt_check.py`](./scripts/rhino_ipopt_check.py).

To find out more about CRA, please refer to our paper in the CAD Computer-Aided Design journal:
[https://doi.org/10.1016/j.cad.2022.103216](https://doi.org/10.1016/j.cad.2022.103216 )

## Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies

![image](https://github.com/BlockResearchGroup/compas_cra/blob/main/docs/assets/images/cra_bridge.png?raw=true)

### Abstract

The rigid-block equilibrium (RBE) method uses a penalty formulation to
measure structural infeasibility or to guide the design of stable
discrete-element assemblies from unstable geometry.
However, RBE is a purely force-based formulation,
and it incorrectly describes stability when
complex interface geometries are involved.
To overcome this issue, this paper introduces
the coupled rigid-block analysis (CRA) method,
a more robust approach building upon RBE’s strengths.
The CRA method combines equilibrium and kinematics in a penalty formulation
in a nonlinear programming problem.
An extensive benchmark campaign is used to show how CRA enables
accurate modelling of complex three-dimensional discrete-element assemblies
formed by rigid blocks.
In addition, an interactive stability-aware design process to
guide user design towards structurally-sound assemblies is proposed.
Finally, the potential of our method for real-world problems are demonstrated
by designing complex and scaffolding-free physical models.

### Please cite our work if you use CRA in your research

#### Paper

```latex
@article{kao2022coupled,
    title     = {Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies},
    author    = {Kao, Gene Ting-Chun and Iannuzzo, Antonino and Thomaszewski, Bernhard and Coros, Stelian and Van Mele, Tom and Block, Philippe},
    journal   = {Computer-Aided Design},
    volume    = {146},
    pages     = {103216},
    year      = {2022},
    publisher = {Elsevier},
    doi       = {10.1016/j.cad.2022.103216},
    url       = {https://doi.org/10.1016/j.cad.2022.103216}
}
```

#### Software implementation

```latex
@misc{compas-cra,
    title  = {{COMPAS CRA}: Coupled Rigid-Block Analysis ({CRA}) for the {COMPAS} framework},
    author = {Kao, Gene Ting-Chun},
    note   = {https://github.com/BlockResearchGroup/compas\_cra},
    year   = {2020-2022},
    doi    = {10.5281/zenodo.7043135},
    url    = {https://doi.org/10.5281/zenodo.7043135},
}
```

### Read the docs

[https://github.com/BlockResearchGroup/compas_cra](https://github.com/BlockResearchGroup/compas_cra)

### Examples to reproduce our paper results

See examples in [docs](https://blockresearchgroup.github.io/compas_cra/examples/) or try them in [docs/examples](https://github.com/BlockResearchGroup/compas_cra/blob/main/docs/examples).
