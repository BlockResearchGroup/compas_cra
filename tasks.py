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

    Windows needs the most: the MSYS2 toolchain paths. macOS needs one thing - the
    same deployment target CI builds IPOPT with (pipeline.yml pins 11.0), on arm64 and
    intel alike; without it the static libs can refuse to link into an extension
    targeting the Python build's older default. Linux needs nothing: the toolchain is
    on PATH and CMakeLists.txt defaults IPOPT_PREFIX to the stage tree by itself.
    """
    system = platform.system()
    if system == "Darwin":
        return {"MACOSX_DEPLOYMENT_TARGET": os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0")}
    if system != "Windows":
        return {}
    root = _msys2_root()
    return {
        # ucrt64\bin on PATH is not optional: g++ spawns cc1plus.exe from ucrt64\lib,
        # and cc1plus finds its DLLs (gmp, mpfr, isl, zstd) in ucrt64\bin. Without it
        # the compiler dies instantly with exit 1 and no output at all, and CMake
        # reports the C++ compiler as "broken".
        "PATH": os.path.join(root, "ucrt64", "bin") + os.pathsep + os.environ.get("PATH", ""),
        "CMAKE_GENERATOR": "Ninja",
    }


VENV_PATCH_MARK = "# compas_cra: MinGW toolchain for rebuilding the extension"


def _patch_venv_activation():
    """Write the build environment into the venv's activation scripts, once.

    After this, a plain `uv pip install -e .` works in the activated venv with no
    special environment: PATH provides the DLLs the compiler's own cc1plus needs, and
    CMAKE_GENERATOR keeps CMake off the Visual Studio generator (MSVC cannot link the
    GCC-built static archives; scikit-build-core also reads it to fetch ninja). The
    compilers themselves are found by CMakeLists.txt. Venv-scoped on purpose - nothing
    global is touched.
    """
    if platform.system() != "Windows":
        return
    scripts = os.path.join(".venv", "Scripts")
    if not os.path.isdir(scripts):
        print("no .venv to patch, skipping activation setup")
        return
    root = _msys2_root()
    ucrt_win = os.path.join(root, "ucrt64", "bin")
    # bash cannot hold C:\...-style entries in PATH (the colon is its separator)
    ucrt_posix = "/" + root.replace(":", "").replace("\\", "/").lower() + "/ucrt64/bin"
    # the venv's own Scripts dir is re-prepended after the toolchain: ucrt64in
    # carries a python.exe of its own, and it must never shadow the venv's
    patches = {
        "activate": [
            VENV_PATCH_MARK,
            'export PATH="%s:$PATH"' % ucrt_posix,
            'export PATH="$VIRTUAL_ENV/Scripts:$PATH"',
            "export CMAKE_GENERATOR=Ninja",
        ],
        "activate.bat": [
            "rem " + VENV_PATCH_MARK[2:],
            'set "PATH=%s;%%PATH%%"' % ucrt_win,
            'set "PATH=%%VIRTUAL_ENV%%\Scripts;%%PATH%%"',
            'set "CMAKE_GENERATOR=Ninja"',
        ],
        "activate.ps1": [
            VENV_PATCH_MARK,
            '$env:PATH = "%s;" + $env:PATH' % ucrt_win,
            '$env:PATH = "$env:VIRTUAL_ENV\Scripts;" + $env:PATH',
            '$env:CMAKE_GENERATOR = "Ninja"',
        ],
    }
    for name, lines in patches.items():
        path = os.path.join(scripts, name)
        if not os.path.isfile(path):
            continue
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
        if "compas_cra: MinGW toolchain" in content:
            continue
        with open(path, "a", encoding="utf-8", newline="") as f:
            f.write("\n" + "\n".join(lines) + "\n")
        print("patched %s" % path)


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
        # build_ipopt.sh resolves Accelerate through xcrun; fail here with a clear
        # message rather than fifteen minutes into the build
        if not shutil.which("xcrun"):
            raise invoke.Exit("Xcode Command Line Tools are required: xcode-select --install")
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


def _verbosity():
    """coinbrew verbosity for interactive builds: stream the compile lines.

    An interactive fifteen-minute build with no output is indistinguishable from a hung
    one, so `invoke setup` defaults to 3. Not 2: coinbrew demotes ThirdParty projects
    one level (`invoke_make $((verbosity-1))`), and ASL and MUMPS - almost the whole
    build - are ThirdParty, so at 2 their make output still goes to /dev/null. At 3
    they stream; only their short autoconf configures stay quiet (those need 4). CI
    calls build_ipopt.sh directly and keeps the script's own quiet default. VERBOSITY
    in the environment still wins, in either direction.
    """
    return os.environ.get("VERBOSITY", "3")


def _build_ipopt(ctx, jobs):
    """Stage IPOPT into build/ipopt/stage, unless it is already there."""
    if os.path.isfile(os.path.join(ctx.base_folder, IPOPT_MARKER)):
        print("IPOPT already staged in %s, skipping." % IPOPT_STAGE)
        return

    print("Building IPOPT. This takes about fifteen minutes, and only has to happen once.")

    if platform.system() == "Windows":
        bash = os.path.join(_msys2_root(), "usr", "bin", "bash.exe")
        # No `cd`, no `&&`, no inline JOBS=: invoke runs this through cmd.exe, which
        # splits on && before bash ever sees the line (single quotes do not protect in
        # cmd). Everything those provided comes through the environment instead:
        # MSYSTEM=UCRT64 with a login shell puts the ucrt64 toolchain on PATH,
        # CHERE_INVOKING keeps the login shell in the invoking directory (we are in
        # base_folder courtesy of the chdir in `setup`), and build_ipopt.sh reads JOBS.
        ctx.run(
            '"{0}" -lc "packaging/build_ipopt.sh"'.format(bash),
            env={"MSYSTEM": "UCRT64", "CHERE_INVOKING": "1", "JOBS": str(jobs), "VERBOSITY": _verbosity()},
        )
    else:
        env = {"JOBS": str(jobs), "VERBOSITY": _verbosity()}
        if platform.system() == "Darwin":
            # same deployment target the CI IPOPT step pins, for arm64 and intel alike
            env["MACOSX_DEPLOYMENT_TARGET"] = os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0")
        ctx.run("packaging/build_ipopt.sh", env=env)


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

        _patch_venv_activation()

        # uv where the contributor is using it, pip otherwise; both install into whatever
        # environment is active
        pip = "uv pip" if shutil.which("uv") else "python -m pip"
        ctx.run('%s install -e ".[dev]"' % pip, env=_native_build_env())

    print("\nDone. `invoke test` should now pass.")


@invoke.task(
    help={
        "serve": "Serve at localhost with live reload (the default). --no-serve builds into dist/docs instead - the artifact CI deploys.",
        "port": "Port to serve on. Defaults to 8000.",
        "strict": "Fail on a broken reference or a missing page.",
    }
)
def docs(ctx, serve=True, port=8000, strict=True):
    """Serve the documentation at http://localhost:8000, rebuilding on every edit.

    Serving is the default because it is the only way the site is browsable locally:
    mkdocs links pages by directory (`examples/`), which only a web server resolves to
    the page inside - opened from file:// every click shows a folder listing. CI runs
    `invoke docs --no-serve` for the deployable build.
    """
    flag = " --strict" if strict else ""
    if serve:
        ctx.run("mkdocs serve --dev-addr localhost:%s%s" % (port, flag))
    else:
        ctx.run("mkdocs build --site-dir dist/docs" + flag)


ns = Collection(
    style.check,
    style.lint,
    style.format,
    docs,
    tests.test,
    tests.testdocs,
    build.prepare_changelog,
    build.clean,
    setup,
    release,
)
ns.configure({"base_folder": os.path.dirname(__file__)})
