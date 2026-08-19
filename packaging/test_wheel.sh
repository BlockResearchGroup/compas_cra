#!/usr/bin/env bash
#
# Install the built wheel into a clean virtual environment and check that a real solve
# works there - i.e. that pyomo finds the bundled ipopt and MUMPS converges. Testing
# the installed wheel rather than the source tree is the point: it is the only way to
# catch a binary that was not packaged, or was packaged without its executable bit.
#
# Usage:  packaging/test_wheel.sh [python]
#
set -euo pipefail

PYTHON="${1:-python3}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
cd "$REPO_ROOT"

WHEEL="$(ls dist/*.whl | head -1)"
[ -n "$WHEEL" ] || { echo "ERROR: no wheel in dist/" >&2; exit 1; }

VENV="$REPO_ROOT/build/venv-test"
rm -rf "$VENV"
"$PYTHON" -m venv "$VENV"

if [ -d "$VENV/Scripts" ]; then
    VBIN="$VENV/Scripts"
else
    VBIN="$VENV/bin"
fi

"$VBIN/python" -m pip install -q --upgrade pip
"$VBIN/python" -m pip install -q "$WHEEL" pytest

echo "--- the ipopt console script is on PATH of the environment"
"$VBIN/ipopt" --version

echo "--- the bundled executable is found from Python"
"$VBIN/python" -c "
from compas_cra._ipopt import bundled, executable, ipopt_version
assert bundled() is not None, 'the wheel does not contain an ipopt executable'
print(executable())
print(ipopt_version())
"

echo "--- solving with the installed wheel, from a directory without the sources"
cd "$REPO_ROOT/build"
"$VBIN/python" -m pytest "$REPO_ROOT/tests" -v

echo "--- ok"
