# The native solver

`native/src/ipopt_nb.cpp` is IPOPT (with the MUMPS linear solver) bound into a Python
extension module with [nanobind](https://github.com/wjakob/nanobind). It is built by
the CMake project at the repository root and installed into the package as
`compas_cra._native`, so `compas_cra` solves CRA problems in-process: no bundled
executable, no subprocess, no `.nl` files, and no separate solver distribution to
install.

| file | what it is |
| --- | --- |
| `src/ipopt_nb.cpp` | the binding: a sparse NLP in numpy arrays and Python callbacks, solved by IPOPT |
| `smoke.py` | minimal numpy-only check that a built wheel loads and solves; run by cibuildwheel on every wheel |

## Building

The extension links the static IPOPT tree produced by `packaging/build_ipopt.sh`, so
that has to exist first:

```bash
packaging/build_ipopt.sh    # ~15 minutes, stages into build/ipopt/stage
pip install .               # picks the stage tree up automatically
python native/smoke.py
```

Point `IPOPT_PREFIX` at a different stage tree to override. Extra link flags (e.g. a
static Fortran runtime) go in `IPOPT_EXTRA_LINK`, extra link directories in
`EXTRA_LINK_DIRS`.

A plain `pip install -e .` needs the same stage tree — an editable install still
compiles the extension. Without it CMake stops with
`No IPOPT build at IPOPT_PREFIX=...`.

## Usage

The `nlp` backend discovers the extension automatically, so nothing here is called
directly:

```python
from compas_cra.equilibrium import cra_solve
cra_solve(assembly)
```
