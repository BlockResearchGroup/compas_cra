# Packaging

`compas_cra` solves its models with [pyomo](https://pyomo.readthedocs.io), which drives
[IPOPT](https://coin-or.github.io/Ipopt/) through the AMPL `.nl` file interface: it
shells out to the `ipopt` **command line executable** rather than linking against the
library. That executable has never been installable with pip - `cyipopt` publishes an
sdist only and needs a system IPOPT to compile against, and coin-or's own releases carry
Windows binaries only - which is why installing `compas_cra` used to require conda, plus
a manual download and copy of solver binaries on Windows.

These scripts remove that requirement by building `ipopt` from source in CI and shipping
it inside the platform wheels.

Because no compiled Python extension is involved, **one `py3-none-<platform>` wheel per
platform serves every supported Python version**:

| wheel | built on |
| --- | --- |
| `py3-none-manylinux_2_28_x86_64` | `quay.io/pypa/manylinux_2_28_x86_64` |
| `py3-none-manylinux_2_28_aarch64` | `quay.io/pypa/manylinux_2_28_aarch64` on an arm runner |
| `py3-none-macosx_11_0_arm64` | `macos-14` |
| `py3-none-macosx_11_0_x86_64` | `macos-13` |
| `py3-none-win_amd64` | `windows-latest` + MSYS2 UCRT64 |

## Scripts

| script | what it does |
| --- | --- |
| `build_ipopt.sh` | builds `ipopt` with coinbrew and stages it into `src/compas_cra/_ipopt/bin/` |
| `check_binary.sh` | runs the binary and asserts it has no non-system dynamic dependencies |
| `build_wheel.sh <tag>` | builds the wheel and forces the `py3-none-<tag>` tag |
| `test_wheel.sh [python]` | installs the wheel into a clean venv and runs the test suite there |

## What goes into the binary

- **IPOPT 3.14.19**, Eclipse Public License 2.0
- **MUMPS** as the linear solver. This is IPOPT's default and the one `compas_cra`
  relies on, since none of the formulations set `linear_solver`.
- **AMPL ASL**, without which the `ipopt` executable is not built at all.
- **BLAS/LAPACK**: serial OpenBLAS on Linux and Windows, Accelerate on macOS.
- **No HSL.** The HSL linear solvers (`ma27`, `ma57`, ...) are not redistributable.

Everything that can be is linked statically, so the shipped binary has no dynamic
dependencies beyond the platform's own system libraries. That is what
`check_binary.sh` enforces, and it is why no `auditwheel`/`delocate` step is needed.

## Building locally

```bash
packaging/build_ipopt.sh                          # ~15 minutes
packaging/build_wheel.sh manylinux_2_28_x86_64    # or your platform tag
packaging/test_wheel.sh
```

On Linux, build inside the manylinux container so the result runs on older
distributions than your own:

```bash
docker run --rm -v "$PWD":/io -w /io quay.io/pypa/manylinux_2_28_x86_64 bash -c '
  dnf install -y --enablerepo=powertools glibc-static libstdc++-static openblas-static
  packaging/build_ipopt.sh'
```

On Windows the build runs in an MSYS2 UCRT64 shell; on macOS it needs `gfortran`
(`brew install gcc`).

## Releasing

`.github/workflows/wheels.yml` builds and tests all five wheels plus an sdist on every
push, and publishes them to PyPI with trusted publishing when a `v*` tag is pushed.
