#!/usr/bin/env bash
#
# Build a platform wheel for the current runner.
#
# The wheel contains no compiled Python extension - only a plain executable as package
# data - so a single wheel works for every supported Python version. Setuptools cannot
# know that, and tags the wheel `py3-none-any`, so the platform tag is applied here.
#
# Usage:  packaging/build_wheel.sh <platform-tag>
#     e.g. packaging/build_wheel.sh manylinux_2_28_x86_64
#
set -euo pipefail

PLATFORM_TAG="${1:?usage: build_wheel.sh <platform-tag>}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
cd "$REPO_ROOT"

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) EXE_NAME=ipopt.exe ;;
    *)                    EXE_NAME=ipopt ;;
esac

BINARY="src/compas_cra/_ipopt/bin/$EXE_NAME"
if [ ! -f "$BINARY" ]; then
    echo "ERROR: $BINARY is missing - run packaging/build_ipopt.sh first" >&2
    exit 1
fi

rm -rf dist build/lib build/bdist.*
python -m build --wheel
python -m wheel tags --remove \
    --python-tag py3 --abi-tag none --platform-tag "$PLATFORM_TAG" dist/*.whl

echo "--- wheel contents"
python -m zipfile --list dist/*.whl | grep -E "_ipopt|RECORD|WHEEL" || true
ls -lh dist/
