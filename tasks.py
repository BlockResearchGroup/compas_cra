from __future__ import print_function

import os
import platform
import shutil

import invoke
from compas_invocations2 import build
from compas_invocations2 import style
from compas_invocations2 import tests
from compas_invocations2.console import chdir
from compas_invocations2.console import confirm
from invoke import Collection


@invoke.task(
    help={"release_type": "Type of release follows semver rules. Must be one of: major, minor, patch, pre_l, pre_n."}
)
def release(ctx, release_type):
    """Releases the project in one swift command!

    This is `compas_invocations2.build.release` without its local `python -m build`
    step. compas_cra is not pure Python any more: building a wheel here needs a staged
    IPOPT tree (packaging/build_ipopt.sh), and CMakeLists.txt hard-fails without one.
    The artifacts that actually get published are built by .github/workflows/pipeline.yml
    on all five platforms anyway, so a local wheel would only be thrown away.
    """
    if release_type not in ("patch", "minor", "major", "pre_l", "pre_n"):
        raise invoke.Exit("The release type parameter is invalid.\nMust be one of: major, minor, patch.")

    ctx.run("invoke format")
    ctx.run("invoke test")

    # bump-my-version also commits and tags, as v{new_version} - which is what
    # .github/workflows/release.yml triggers on
    ctx.run("bump-my-version bump %s --verbose" % release_type)

    build.prepare_changelog(ctx)
    # builds=False: the default `clean` removes build/, which here holds the staged IPOPT
    # tree packaging/build_ipopt.sh spends about fifteen minutes producing. CMakeLists.txt
    # defaults IPOPT_PREFIX to build/ipopt/stage, so removing it makes the next
    # `pip install -e .` fail with `No IPOPT build at IPOPT_PREFIX=...`.
    build.clean(ctx, builds=False)

    if confirm(
        "Everything is ready. You are about to push to git, which will build every wheel "
        "and release to pypi.org. Are you sure?",
        assume_yes=False,
    ):
        ctx.run("git push --tags && git push")
    else:
        raise invoke.Exit("You need to manually revert the tag/commits created.")


# ============================================================================
# building the solver
# ============================================================================
# compas_cra is not pure Python: `pip install -e .` compiles compas_cra._native against
# a static IPOPT tree that packaging/build_ipopt.sh stages into build/ipopt/stage. The
# `setup` task below is that whole sequence - toolchain, IPOPT, editable install - as
# one command, on all three platforms.

IPOPT_STAGE = os.path.join("build", "ipopt", "stage")
# what CMakeLists.txt itself probes for, so "staged" here means exactly what it means there
IPOPT_MARKER = os.path.join(IPOPT_STAGE, "include", "coin-or", "IpTNLP.hpp")

MSYS2_PACKAGES = (
    "git patch make diffutils pkgconf mingw-w64-ucrt-x86_64-gcc "
    "mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-openblas"
)

LINUX_HINT = """No Fortran compiler found. The IPOPT build needs one (MUMPS is Fortran),
plus a static BLAS/LAPACK and a C/C++ toolchain:

  Debian / Ubuntu:
    sudo apt install build-essential gfortran libopenblas-dev git curl make patch pkg-config

  Fedora / RHEL / AlmaLinux:
    sudo dnf install gcc gcc-c++ gcc-gfortran openblas-devel openblas-static \
        glibc-static libstdc++-static git curl make patch pkgconf

These need root, so they are not installed for you. Re-run `invoke setup` afterwards."""


def _msys2_root():
    """Locate the MSYS2 installation the Windows IPOPT build runs inside."""
    candidates = [os.environ.get("MSYS2_ROOT"), r"C:\msys64", r"D:\msys64"]
    for root in candidates:
        if root and os.path.isfile(os.path.join(root, "usr", "bin", "bash.exe")):
            return root
    raise invoke.Exit(
        "No MSYS2 installation found. The Windows IPOPT build runs in an MSYS2 UCRT64 "
        "shell; install it from https://www.msys2.org, or set MSYS2_ROOT if it lives "
        r"somewhere other than C:\msys64 or D:\msys64."
    )


def _native_build_env():
    """Environment `pip install -e .` needs to find the compilers and the runtimes.

    Only Windows needs anything: elsewhere the toolchain is on PATH and CMakeLists.txt
    defaults IPOPT_PREFIX to the stage tree by itself.
    """
    if platform.system() != "Windows":
        return {}
    root = _msys2_root()
    ucrt = root.replace("\\", "/") + "/ucrt64"
    return {
        "EXTRA_LINK_DIRS": os.path.join(root, "ucrt64", "lib"),
        "CMAKE_GENERATOR": "Ninja",
        "CMAKE_ARGS": (
            "-DCMAKE_C_COMPILER={0}/bin/gcc.exe "
            "-DCMAKE_CXX_COMPILER={0}/bin/g++.exe "
            "-DFORTRAN_COMPILER={0}/bin/gfortran.exe".format(ucrt)
        ),
    }


def _install_toolchain(ctx):
    """Install what build_ipopt.sh needs, where that can be done without root."""
    system = platform.system()

    if system == "Windows":
        pacman = os.path.join(_msys2_root(), "usr", "bin", "pacman.exe")
        # --needed makes this a no-op once installed, --noconfirm keeps it non-interactive
        ctx.run('"%s" -S --needed --noconfirm %s' % (pacman, MSYS2_PACKAGES))

    elif system == "Darwin":
        if not shutil.which("brew"):
            raise invoke.Exit("Homebrew is required on macOS: https://brew.sh")
        # coinbrew refuses to run under bash 3, which is what macOS still ships
        ctx.run("brew list bash >/dev/null 2>&1 || brew install bash")
        ctx.run("brew list gcc >/dev/null 2>&1 || brew install gcc")
        # gfortran arrives from the gcc formula as gfortran-14 or similar; build_ipopt.sh
        # looks for the unsuffixed name
        ctx.run(
            'test -x "$(brew --prefix)/bin/gfortran" || '
            'ln -sf "$(ls "$(brew --prefix)"/bin/gfortran-* | tail -1)" "$(brew --prefix)/bin/gfortran"'
        )

    elif not shutil.which("gfortran"):
        raise invoke.Exit(LINUX_HINT)


def _build_ipopt(ctx, jobs):
    """Stage IPOPT into build/ipopt/stage, unless it is already there."""
    if os.path.isfile(os.path.join(ctx.base_folder, IPOPT_MARKER)):
        print("IPOPT already staged in %s, skipping." % IPOPT_STAGE)
        return

    print("Building IPOPT. This takes about fifteen minutes, and only has to happen once.")

    if platform.system() == "Windows":
        bash = os.path.join(_msys2_root(), "usr", "bin", "bash.exe")
        # cygpath rather than string surgery on the drive letter: MSYS2 can be mounted
        # with a prefix other than /c, and only it knows
        unix_path = ctx.run(
            '"{0}" -c \'cygpath -u "{1}"\''.format(bash, ctx.base_folder), hide=True
        ).stdout.strip()
        # MSYSTEM=UCRT64 with a login shell is what puts the ucrt64 toolchain on PATH
        ctx.run(
            '"{0}" -lc \'cd "{1}" && JOBS={2} packaging/build_ipopt.sh\''.format(bash, unix_path, jobs),
            env={"MSYSTEM": "UCRT64", "CHERE_INVOKING": "1"},
        )
    else:
        ctx.run("packaging/build_ipopt.sh", env={"JOBS": str(jobs)})


@invoke.task(
    help={
        "jobs": (
            "Parallel jobs for the IPOPT build. Defaults to 1, and should stay there: "
            "MUMPS' makefiles are missing dependencies, so a parallel build races and "
            "dies linking libcoinmumps, intermittently and with no error output."
        ),
        "toolchain": "False to skip the toolchain check and install.",
    }
)
def setup(ctx, jobs=1, toolchain=True):
    """Build the solver and install compas_cra into the active environment.

    Everything docs/installation.md spells out by hand, in one command, on Windows,
    macOS and Linux: install the toolchain where that is possible without root, stage
    IPOPT if it is not staged already, then install the package editable with the
    compilers and link directories the extension needs.

    Idempotent - a second run reinstalls the package and leaves the IPOPT tree alone.
    """
    with chdir(ctx.base_folder):
        if toolchain:
            _install_toolchain(ctx)

        _build_ipopt(ctx, jobs)

        # uv where the contributor is using it, pip otherwise; both install into whatever
        # environment is active
        pip = "uv pip" if shutil.which("uv") else "python -m pip"
        ctx.run('%s install -e ".[dev]"' % pip, env=_native_build_env())

    print("\nDone. `invoke test` should now pass.")


@invoke.task(help={"strict": "Fail the build on a broken reference or a missing page."})
def docs(ctx, strict=True):
    """Build the documentation with mkdocs into dist/docs."""
    ctx.run("mkdocs build --site-dir dist/docs" + (" --strict" if strict else ""))


@invoke.task(help={"port": "Port to serve on. Defaults to 8000."})
def docs_serve(ctx, port=8000):
    """Serve the documentation locally, rebuilding on every change."""
    ctx.run("mkdocs serve --dev-addr localhost:%s" % port)


ns = Collection(
    style.check,
    style.lint,
    style.format,
    docs,
    docs_serve,
    tests.test,
    tests.testdocs,
    build.prepare_changelog,
    build.clean,
    setup,
    release,
)
ns.configure({"base_folder": os.path.dirname(__file__)})
