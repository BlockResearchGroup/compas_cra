#!/usr/bin/env bash
#
# Audit a staged `ipopt` binary: it must run, and it must not depend on anything
# beyond the platform's own system libraries. A dependency on libgfortran, MUMPS,
# OpenBLAS or a mingw runtime DLL means the wheel would break on a user's machine.
#
set -euo pipefail

BINARY="${1:?usage: check_binary.sh <path to ipopt>}"

echo "--- $BINARY"
ls -lh "$BINARY"

echo "--- version"
"$BINARY" --version

echo "--- dynamic dependencies"
case "$(uname -s)" in
    Linux*)
        if ldd "$BINARY" 2>&1 | grep -q "not a dynamic executable"; then
            echo "statically linked - no dynamic dependencies"
        else
            ldd "$BINARY"
            # Only the C library and friends may stay dynamic. Anything else - a Fortran
            # runtime, a BLAS, libstdc++ - would have to be present on the user's machine.
            if ldd "$BINARY" | awk '{print $1}' | grep -vE \
                "^(linux-vdso\.so\.1|libc\.so\.6|libm\.so\.6|libdl\.so\.2|librt\.so\.1|libpthread\.so\.0|libz\.so\.1|/lib64/ld-linux.*|/lib/ld-linux.*)$" \
                | grep -q .; then
                echo "ERROR: binary depends on a library outside the system set" >&2
                exit 1
            fi
        fi
        ;;
    Darwin*)
        otool -L "$BINARY"
        if otool -L "$BINARY" | tail -n +2 | grep -vqE "/usr/lib/|/System/Library/"; then
            echo "ERROR: binary depends on a library outside /usr/lib and /System/Library" >&2
            exit 1
        fi
        ;;
    MINGW*|MSYS*|CYGWIN*)
        # objdump is always present in MSYS2; dumpbin is not.
        objdump -p "$BINARY" | grep "DLL Name:" | sort -u
        if objdump -p "$BINARY" | grep -i "DLL Name:" | grep -Eqi "libgcc|libstdc\+\+|libgfortran|libquadmath|libwinpthread|libgomp"; then
            echo "ERROR: binary depends on a mingw runtime DLL" >&2
            exit 1
        fi
        ;;
esac

echo "--- ok"
