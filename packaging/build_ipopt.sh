#!/usr/bin/env bash
#
# Build a self-contained `ipopt` executable for the current platform and stage it
# into src/compas_cra/_ipopt/bin/.
#
# COMPAS CRA reaches IPOPT exclusively through `pyomo.SolverFactory("ipopt")`, which
# shells out to the IPOPT command line executable (the AMPL/.nl interface). So all we
# need to ship is that one executable - no Python bindings, no CPython ABI involved.
#
# The build uses coinbrew with:
#   - ThirdParty-ASL    (mandatory: the `ipopt` executable only exists with ASL)
#   - ThirdParty-Mumps  (the default linear solver, which is what CRA relies on)
#   - the platform's own static BLAS/LAPACK (see LAPACK_FLAGS below)
#   - no HSL            (not redistributable)
#
# Everything is linked statically so the resulting wheel is a single self-contained
# file: no rpath surgery, no bundled runtime DLLs, no auditwheel/delocate pass.
#
# Usage:  packaging/build_ipopt.sh
# Env:    IPOPT_VERSION, WORK_DIR, DEST_DIR, JOBS
#
set -euo pipefail

IPOPT_VERSION="${IPOPT_VERSION:-3.14.19}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
WORK_DIR="${WORK_DIR:-$REPO_ROOT/build/ipopt}"
STAGE_DIR="$WORK_DIR/stage"
DEST_DIR="${DEST_DIR:-$REPO_ROOT/src/compas_cra/_ipopt/bin}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}"

case "$(uname -s)" in
    Linux*)               PLATFORM=linux;   EXE_NAME=ipopt ;;
    Darwin*)              PLATFORM=macos;   EXE_NAME=ipopt ;;
    MINGW*|MSYS*|CYGWIN*) PLATFORM=windows; EXE_NAME=ipopt.exe ;;
    *) echo "unsupported platform: $(uname -s)" >&2; exit 1 ;;
esac

# BLAS/LAPACK comes from the platform, statically linked: the serial OpenBLAS archive
# on Linux and Windows (libopenblas.a provides LAPACK as well; the serial rather than the
# pthread/openmp build, since MUMPS is called sequentially here), and the system
# Accelerate framework on macOS.
#
# A static libopenblas.a also needs libgfortran, libquadmath, libm and libpthread, and
# they have to come *after* it on the link line. Neither LIBS nor a multi-word
# --with-lapack survives coinbrew (it clears the first and splits the second on
# whitespace), so on Linux the whole chain is wrapped in a GNU ld linker script that
# looks like an archive and is therefore a single token. See make_lapack_shim below.

LAPACK_SHIM_DIR="$WORK_DIR/lapack-shim"

make_lapack_shim() {
    # A text file named lib<name>.a is read by GNU ld as a script, so -lcralapack pulls
    # in every archive listed here, in this order.
    local openblas
    openblas="$(${CC:-gcc} -print-file-name=libopenblas.a)"
    if [ ! -f "$openblas" ]; then
        echo "ERROR: no static openblas found (install openblas-static)" >&2
        exit 1
    fi
    local threads="-lpthread"
    if [ "$PLATFORM" = "windows" ]; then
        threads=""
    fi
    mkdir -p "$LAPACK_SHIM_DIR"
    echo "INPUT($openblas$RUNTIME_LINK_LIBS -lm $threads)" \
        > "$LAPACK_SHIM_DIR/libcralapack.a"
    echo "  lapack shim: $(cat "$LAPACK_SHIM_DIR/libcralapack.a")"
}

case "$PLATFORM" in
    linux|windows)
        LAPACK_FLAGS="${LAPACK_FLAGS:--lcralapack}"
        LAPACK_LABEL="serial OpenBLAS, statically linked        BSD 3-clause"
        ;;
    macos)
        LAPACK_FLAGS="${LAPACK_FLAGS:--Wl,-framework,Accelerate}"
        LAPACK_LABEL="Apple Accelerate framework"
        ;;
esac

echo "=============================================================================="
echo " building ipopt $IPOPT_VERSION for $PLATFORM ($(uname -m)) with $JOBS jobs"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# static linking flags
# ------------------------------------------------------------------------------
# The goal is a binary that runs on a machine which has never seen a Fortran compiler.
# Plain `-static` is not an option: the final link goes through libtool, which swallows
# it as one of its own flags instead of passing it to the linker. So the compiler
# runtimes are pinned statically one by one instead:
#
#   -static-libgcc / -static-libstdc++   understood by gcc, passed through by libtool
#   -L<dir of .a only>                   ld takes the first match in -L order, so a
#                                        directory holding only libgfortran.a and
#                                        libquadmath.a forces those to link statically
#
# What remains dynamic is the platform's own C library and friends (libc, libm, libdl,
# libpthread, libz), which exist everywhere. check_binary.sh enforces exactly that.

# -static-libgcc/-static-libstdc++ are gcc options. On macOS the C compiler is clang,
# which rejects them outright - and does not need them: libc++ is a system library
# there, so only the Fortran runtime has to be pinned, which the -L trick below does.
if [ "$PLATFORM" = "macos" ]; then
    STATIC_LDFLAGS=""
else
    STATIC_LDFLAGS="-static-libgcc -static-libstdc++"
fi

# Which of these exist is platform dependent - libquadmath is x86 only, libwinpthread
# is mingw only - so the list is discovered rather than assumed, and reused when the
# lapack shim is written. libgomp is in there because MSYS2 ships the OpenMP build of
# OpenBLAS; on the other platforms nothing references it and ld pulls in nothing.
RUNTIME_LINK_LIBS=""

force_static_runtime_libs() {
    local libdir="$WORK_DIR/static-runtime"
    rm -rf "$libdir" && mkdir -p "$libdir"
    local lib path
    for lib in libgfortran libquadmath libgcc libemutls_w libgomp libwinpthread; do
        path="$(${FC:-gfortran} -print-file-name=${lib}.a)"
        if [ -f "$path" ]; then
            ln -sf "$path" "$libdir/${lib}.a"
            RUNTIME_LINK_LIBS="$RUNTIME_LINK_LIBS -l${lib#lib}"
        fi
    done
    echo "  static runtime libs:$RUNTIME_LINK_LIBS"
    STATIC_LDFLAGS="$STATIC_LDFLAGS -L$libdir"
}

force_static_runtime_libs

if [ "$PLATFORM" = "macos" ]; then
    : "${MACOSX_DEPLOYMENT_TARGET:?MACOSX_DEPLOYMENT_TARGET must be set}"
    export MACOSX_DEPLOYMENT_TARGET
fi

# ------------------------------------------------------------------------------
# fetch
# ------------------------------------------------------------------------------
# coinbrew refuses to run under bash 3, which is still what /bin/bash is on macOS.
find_modern_bash() {
    local candidate version
    for candidate in bash /opt/homebrew/bin/bash /usr/local/bin/bash; do
        version="$("$candidate" -c 'echo ${BASH_VERSINFO[0]}' 2>/dev/null)" || continue
        if [ "${version:-0}" -ge 4 ]; then
            echo "$candidate"
            return 0
        fi
    done
    echo "ERROR: coinbrew needs bash >= 4; install one (brew install bash)" >&2
    exit 1
}

BASH_BIN="${BASH_BIN:-$(find_modern_bash)}"
echo "  coinbrew runs under: $BASH_BIN ($("$BASH_BIN" --version | head -1))"

# coinbrew drives git, which refuses to touch a tree owned by another user - which is
# what happens as soon as the build runs in a container over a mounted checkout. Set the
# exception through the environment so no global git config is written.
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0=safe.directory
export GIT_CONFIG_VALUE_0="*"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

if [ ! -f ./coinbrew ]; then
    curl -fsSL -o coinbrew \
        https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
fi

# Pulls Ipopt itself plus ThirdParty/ASL and ThirdParty/Mumps, whose get.* scripts
# download the ASL and MUMPS sources. ThirdParty/HSL is cloned too but coinbrew skips
# building it because we never provide the (non-redistributable) HSL sources.
if [ ! -d Ipopt ]; then
    "$BASH_BIN" ./coinbrew fetch "Ipopt@releases/$IPOPT_VERSION" --no-prompt --skip-update
fi

# ------------------------------------------------------------------------------
# build
# ------------------------------------------------------------------------------
COMMON_CONFIG=(
    --no-prompt
    --skip-update
    --prefix="$STAGE_DIR"
    --parallel-jobs="$JOBS"
    --tests=none
    --disable-shared
    --enable-static
)

if [ "$LAPACK_FLAGS" = "-lcralapack" ]; then
    make_lapack_shim
    STATIC_LDFLAGS="$STATIC_LDFLAGS -L$LAPACK_SHIM_DIR"
fi

export LDFLAGS="${LDFLAGS:-} $STATIC_LDFLAGS"

"$BASH_BIN" ./coinbrew build Ipopt "${COMMON_CONFIG[@]}" --with-lapack="$LAPACK_FLAGS"

# ------------------------------------------------------------------------------
# stage
# ------------------------------------------------------------------------------
collect_licenses() {
    # Everything bundled must carry its license: Ipopt is EPL-2.0, MUMPS and the ASL
    # are permissive, OpenBLAS is BSD-3, and the statically linked libgfortran/libgcc
    # are covered by the GCC Runtime Library Exception.
    local dest="$1"
    mkdir -p "$dest"
    local sources=(
        "Ipopt/LICENSE:Ipopt-EPL-2.0.txt"
        "ThirdParty/Mumps/MUMPS/LICENSE:MUMPS-LICENSE.txt"
        "ThirdParty/ASL/LICENSE:ASL-LICENSE.txt"
    )
    # Accelerate is a system framework, so only the OpenBLAS builds bundle its license.
    if [ "$LAPACK_FLAGS" = "-lcralapack" ]; then
        sources+=("$HERE/licenses/OpenBLAS-LICENSE.txt:OpenBLAS-LICENSE.txt")
    fi
    local src
    for src in "${sources[@]}"
    do
        local from="${src%%:*}" to="${src##*:}"
        if [ -f "$WORK_DIR/$from" ]; then
            cp "$WORK_DIR/$from" "$dest/$to"
        elif [ -f "$from" ]; then
            cp "$from" "$dest/$to"
        fi
    done
    cat > "$dest/README.txt" <<EOF
The ipopt executable shipped in ../bin was built from source by
packaging/build_ipopt.sh in the compas_cra repository.

    IPOPT     $IPOPT_VERSION   Eclipse Public License 2.0    https://github.com/coin-or/Ipopt
    MUMPS               linear solver, permissive     https://mumps-solver.org
    AMPL ASL            solver library, permissive    https://github.com/ampl/asl
    BLAS/LAPACK         $LAPACK_LABEL

The GNU Fortran and C++ runtimes are linked statically under the GCC Runtime Library
Exception. HSL linear solvers are deliberately not included: they are not
redistributable. Sources for the EPL-2.0 licensed parts are available from the URLs
above at the versions stated.
EOF
    ls "$dest"
}

BUILT="$STAGE_DIR/bin/$EXE_NAME"
[ -f "$BUILT" ] || { echo "ERROR: $BUILT was not produced" >&2; exit 1; }

mkdir -p "$DEST_DIR"
cp "$BUILT" "$DEST_DIR/$EXE_NAME"
chmod +x "$DEST_DIR/$EXE_NAME"

# An unstripped ipopt is ~50 MB, most of it debug symbols nobody can act on from a
# wheel. Stripping takes it to a few MB, which matters when it ships in every wheel.
echo "  size before strip: $(du -h "$DEST_DIR/$EXE_NAME" | cut -f1)"
if [ "$PLATFORM" = "macos" ]; then
    strip -x "$DEST_DIR/$EXE_NAME"
else
    strip "$DEST_DIR/$EXE_NAME"
fi
echo "  size after strip:  $(du -h "$DEST_DIR/$EXE_NAME" | cut -f1)"

# The linker ad-hoc signs arm64 Mach-O binaries automatically, but a `cp` across
# filesystems has been known to strip it; re-signing is cheap and idempotent.
if [ "$PLATFORM" = "macos" ]; then
    codesign --force --sign - "$DEST_DIR/$EXE_NAME" || true
fi

collect_licenses "$(dirname "$DEST_DIR")/licenses"

"$HERE/check_binary.sh" "$DEST_DIR/$EXE_NAME"

echo
echo "staged $DEST_DIR/$EXE_NAME"
