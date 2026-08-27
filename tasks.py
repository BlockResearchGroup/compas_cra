from __future__ import print_function

import glob
import os
import platform
import re
import shutil
import sysconfig

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


def _ipopt_work_dir(base_folder):
    """Where build_ipopt.sh builds and stages IPOPT: build/ipopt, space permitting.

    IPOPT is autotools, and autotools flatly refuses a source or prefix path containing
    whitespace - configure stops at `unsafe srcdir value`, and coinbrew's own test
    expressions fall apart the same way. A checkout under, say, "Desktop/untitled
    folder" can never build IPOPT in-tree. The script takes WORK_DIR from the
    environment and CMakeLists.txt takes IPOPT_PREFIX, so in that case the build goes
    to ~/.compas_cra/ipopt instead - shared between such checkouts on purpose, the
    staged tree is a function of the pinned IPOPT version, not of the clone.
    """
    default = os.path.join(base_folder, "build", "ipopt")
    if not re.search(r"\s", default):
        return default
    fallback = os.path.join(os.path.expanduser("~"), ".compas_cra", "ipopt")
    if re.search(r"\s", fallback):
        raise invoke.Exit(
            "Both the repository path and %s contain whitespace, and the IPOPT build is "
            "autotools, which cannot build under either. Move the repository to a path "
            "without spaces, or set WORK_DIR to one." % fallback
        )
    print("NOTE: the repository path contains a space, which autotools cannot build under.")
    print("      IPOPT will be built and staged in %s instead." % fallback)
    return fallback


def _ipopt_stage_dir(base_folder):
    return os.path.join(_ipopt_work_dir(base_folder), "stage")


def _ipopt_marker(base_folder):
    """What CMakeLists.txt itself probes for, so "staged" means exactly what it means there."""
    return os.path.join(_ipopt_stage_dir(base_folder), "include", "coin-or", "IpTNLP.hpp")


def _ipopt_prefix_override(base_folder):
    """The IPOPT_PREFIX the extension build needs, or None when the default is right.

    None whenever the stage tree is at build/ipopt/stage, where CMakeLists.txt looks by
    itself; a path only when _ipopt_work_dir redirected the build away from a
    whitespace repository path.
    """
    stage = _ipopt_stage_dir(base_folder)
    if stage == os.path.join(base_folder, "build", "ipopt", "stage"):
        return None
    return stage


MSYS2_PACKAGES = (
    "git patch make diffutils pkgconf mingw-w64-ucrt-x86_64-gcc "
    "mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-openblas"
)

# The official installer. Run with NONINTERACTIVE=1 it asks nothing, but it still needs
# sudo once to create the prefix, so macOS prompts for the login password a single time.
# Set COMPAS_CRA_NO_BREW_INSTALL to refuse the bootstrap and just fail with a hint.
HOMEBREW_INSTALL_URL = "https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh"
# where brew lands: Apple silicon first, then Intel. Checked directly because a brew that
# was installed moments ago - or by someone who never edited their shell profile - is on
# disk but not yet on PATH.
HOMEBREW_PREFIXES = ("/opt/homebrew", "/usr/local")

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


def _brew():
    """Locate brew, on PATH or at either standard prefix. None if it is not installed."""
    found = shutil.which("brew")
    if found:
        return found
    for prefix in HOMEBREW_PREFIXES:
        # these are macOS paths by definition; os.path.join would backslash them on a
        # Windows machine running the test suite
        candidate = prefix + "/bin/brew"
        if os.access(candidate, os.X_OK):
            return candidate
    return None


def _use_brew(brew):
    """Put brew's bin directory on PATH for the rest of this process.

    A brew installed by _install_homebrew during this very run is not on the PATH we
    inherited, and ctx.run passes os.environ down, so without this every later `brew`
    call would fail until the contributor opened a new terminal.
    """
    bin_dir = os.path.dirname(brew)
    path = os.environ.get("PATH", "")
    if bin_dir not in path.split(os.pathsep):
        os.environ["PATH"] = bin_dir + os.pathsep + path
    return brew


def _refuse_root():
    """Stop if the task is being run under sudo.

    `sudo invoke setup` is the obvious thing to try after the Homebrew installer asks for
    administrator rights, and it is worse than the problem: Homebrew flatly refuses to
    install as root, and were it not for that, build/, the staged IPOPT tree and the venv
    would all come out owned by root and unwritable afterwards. Correct is to run as
    yourself; the Homebrew install asks for the password itself.
    """
    # no geteuid on Windows, where there is nothing equivalent to guard against
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        raise invoke.Exit(
            "Do not run `invoke setup` with sudo. Homebrew refuses to install as root, and "
            "the IPOPT build, the venv and build/ would end up owned by root. Run it as "
            "yourself - it asks for your password at the one point it needs it."
        )


def _install_homebrew(ctx):
    """Install Homebrew, and return the path to the brew it produced.

    macOS ships no package manager and no Fortran compiler, so there is nothing to build
    IPOPT with until this has run. It is the one step that needs root - the installer
    creates /opt/homebrew and chowns it - which is a single password prompt, not a
    separate thing to go and do by hand. pty=True so that prompt is actually visible.
    """
    if os.environ.get("COMPAS_CRA_NO_BREW_INSTALL"):
        raise invoke.Exit(
            "Homebrew is required on macOS and COMPAS_CRA_NO_BREW_INSTALL is set. "
            "Install it from https://brew.sh and re-run `invoke setup`."
        )
    print("Homebrew not found. Installing it - macOS asks for your login password once.")
    installer = '/bin/bash -c "$(curl -fsSL %s)"' % HOMEBREW_INSTALL_URL
    # `sudo -v` has to share a pty with the installer, hence one ctx.run and not two.
    # macOS sudo keys its credential cache to the controlling tty, and pty=True allocates
    # a fresh pty per call - so a password entered in an earlier ctx.run is invisible to
    # the installer's own `sudo -n` probe, which then reports the account is not an
    # administrator when it is. One run, one tty, one password.
    result = ctx.run("sudo -v && " + installer, env={"NONINTERACTIVE": "1"}, pty=True, warn=True)
    if not result.ok and not _brew():
        # NONINTERACTIVE makes every prompt fatal, a timestamp that expired during a slow
        # download included. Let the installer do its own asking; it wants a RETURN too.
        print("\nNon-interactive install did not finish. Retrying - press RETURN when asked.")
        ctx.run(installer, pty=True, warn=True)
    brew = _brew()
    if not brew:
        raise invoke.Exit(
            "The Homebrew installer finished but no brew turned up at %s. "
            "Install it manually from https://brew.sh and re-run `invoke setup`." % " or ".join(HOMEBREW_PREFIXES)
        )
    print("Homebrew installed at %s" % brew)
    return brew


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
VENV_PATCH_MARK_MACOS = "# compas_cra: Homebrew toolchain for rebuilding the extension"
# both marks share this prefix, so one containment test keeps either patch from being
# appended twice - including to venvs patched by the Windows-only version of this task
VENV_PATCH_SENTINEL = "compas_cra:"


def _append_once(path, lines):
    """Append lines to a shell script, unless this task has already written to it."""
    if not os.path.isfile(path):
        return
    with open(path, "r", encoding="utf-8") as f:
        if VENV_PATCH_SENTINEL in f.read():
            return
    with open(path, "a", encoding="utf-8", newline="") as f:
        f.write("\n" + "\n".join(lines) + "\n")
    print("patched %s" % path)


def _patch_venv_activation_macos():
    """Put Homebrew and the deployment target into the venv's activate script.

    The macOS mirror of the Windows patch below, and for the same reason: after it, a
    plain `uv pip install -e .` in the activated venv rebuilds the extension with no
    special environment. It matters more here than it looks, because `invoke setup` may
    have installed Homebrew moments ago into a shell profile the current terminal never
    read. Venv-scoped on purpose - no dotfile outside the repository is touched.
    """
    activate = os.path.join(".venv", "bin", "activate")
    if not os.path.isfile(activate):
        print("no .venv to patch, skipping activation setup")
        return
    brew = _brew()
    if not brew:
        return
    lines = [
        VENV_PATCH_MARK_MACOS,
        'eval "$(%s shellenv)"' % brew,
        # brew shellenv prepends its own bin, which carries a python of its own; the
        # venv must win, exactly as Scripts must on Windows
        'export PATH="$VIRTUAL_ENV/bin:$PATH"',
        "export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}",
    ]
    stage = _ipopt_prefix_override(os.getcwd())
    if stage:
        lines.append('export IPOPT_PREFIX="%s"' % stage)
    _append_once(activate, lines)


def _patch_venv_activation():
    """Write the build environment into the venv's activation scripts, once.

    Windows below, macOS in _patch_venv_activation_macos, nothing on Linux - there the
    toolchain is already on PATH and needs no help.

    After this, a plain `uv pip install -e .` works in the activated venv with no
    special environment: PATH provides the DLLs the compiler's own cc1plus needs, and
    CMAKE_GENERATOR keeps CMake off the Visual Studio generator (MSVC cannot link the
    GCC-built static archives; scikit-build-core also reads it to fetch ninja). The
    compilers themselves are found by CMakeLists.txt. Venv-scoped on purpose - nothing
    global is touched.
    """
    system = platform.system()
    if system == "Darwin":
        _patch_venv_activation_macos()
        return
    if system != "Windows":
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
            'set "PATH=%%VIRTUAL_ENV%%\\Scripts;%%PATH%%"',
            'set "CMAKE_GENERATOR=Ninja"',
        ],
        "activate.ps1": [
            VENV_PATCH_MARK,
            '$env:PATH = "%s;" + $env:PATH' % ucrt_win,
            '$env:PATH = "$env:VIRTUAL_ENV\\Scripts;" + $env:PATH',
            '$env:CMAKE_GENERATOR = "Ninja"',
        ],
    }
    stage = _ipopt_prefix_override(os.getcwd())
    if stage:
        patches["activate"].append('export IPOPT_PREFIX="%s"' % _msys_path(stage))
        patches["activate.bat"].append('set "IPOPT_PREFIX=%s"' % stage)
        patches["activate.ps1"].append('$env:IPOPT_PREFIX = "%s"' % stage)
    for name, lines in patches.items():
        _append_once(os.path.join(scripts, name), lines)


def _install_toolchain(ctx):
    """Install what build_ipopt.sh needs, per platform.

    Windows and macOS are both bootstrapped end to end - pacman packages into an existing
    MSYS2 on the one, Homebrew and its formulae on the other, installing Homebrew itself
    if it is missing. Linux is the exception: its package managers need root for every
    install, not just the first, so it prints LINUX_HINT and stops.
    """
    system = platform.system()

    if system == "Windows":
        pacman = os.path.join(_msys2_root(), "usr", "bin", "pacman.exe")
        # --needed makes this a no-op once installed, --noconfirm keeps it non-interactive
        ctx.run('"%s" -S --needed --noconfirm %s' % (pacman, MSYS2_PACKAGES))

    elif system == "Darwin":
        # installed on demand: a stock macOS has no package manager, and requiring the
        # contributor to go and get one first is the one manual step this task exists to
        # remove. Addressed by absolute path throughout - a brew installed just now is
        # not on the PATH this process started with.
        brew = _use_brew(_brew() or _install_homebrew(ctx))
        # build_ipopt.sh resolves Accelerate through xcrun. The Homebrew installer pulls
        # the Command Line Tools in itself, so by here this should hold; it is checked
        # anyway rather than failing fifteen minutes into the build.
        if not shutil.which("xcrun"):
            raise invoke.Exit("Xcode Command Line Tools are required: xcode-select --install")
        # coinbrew refuses to run under bash 3, which is what macOS still ships
        ctx.run('"{0}" list bash >/dev/null 2>&1 || "{0}" install bash'.format(brew))
        ctx.run('"{0}" list gcc >/dev/null 2>&1 || "{0}" install gcc'.format(brew))
        # not strictly required, but coinbrew warns loudly without it, and configure
        # resolves dependencies more reliably with it; MSYS2 and the Linux hint both
        # already include their pkgconf
        ctx.run('"{0}" list pkgconf >/dev/null 2>&1 || "{0}" install pkgconf'.format(brew))
        # gfortran arrives from the gcc formula as gfortran-14 or similar; build_ipopt.sh
        # looks for the unsuffixed name
        ctx.run(
            'test -x "$("{0}" --prefix)/bin/gfortran" || '
            'ln -sf "$(ls "$("{0}" --prefix)"/bin/gfortran-* | tail -1)" "$("{0}" --prefix)/bin/gfortran"'.format(brew)
        )

    elif not shutil.which("gfortran"):
        raise invoke.Exit(LINUX_HINT)


def _codesign_extension(ctx):
    """Ad-hoc re-sign the freshly installed extension modules. Apple silicon only.

    CMake's install step edits the load commands of the .so after the linker ad-hoc
    signed it, which invalidates the signature - and arm64 macOS answers an invalid
    signature at import with SIGKILL: the process dies with exit 137 and no message at
    all. The wheels CI publishes are re-signed by delocate; the local editable install
    gets the same treatment here. Idempotent, codesign --force replaces cheerfully.
    """
    if platform.system() != "Darwin":
        return
    pattern = os.path.join(sysconfig.get_path("platlib"), "compas_cra", "**", "*.so")
    for so in sorted(glob.glob(pattern, recursive=True)):
        ctx.run('codesign --force --sign - "%s"' % so)


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


def _msys_path(path):
    """C:\\a\\b as MSYS2 bash needs it: /c/a/b. Colon and backslashes both break bash.

    By hand rather than os.path.splitdrive, which only understands drive letters when
    running on Windows itself - this way the helper means the same thing everywhere.
    """
    drive, colon, rest = path.partition(":")
    if colon and len(drive) == 1:
        return "/" + drive.lower() + rest.replace("\\", "/")
    return path.replace("\\", "/")


def _build_ipopt(ctx, jobs):
    """Stage IPOPT, in-tree or in the whitespace fallback, unless already staged."""
    if os.path.isfile(_ipopt_marker(ctx.base_folder)):
        print("IPOPT already staged in %s, skipping." % _ipopt_stage_dir(ctx.base_folder))
        return

    work_dir = _ipopt_work_dir(ctx.base_folder)
    in_tree = work_dir == os.path.join(ctx.base_folder, "build", "ipopt")

    print("Building IPOPT. This takes about fifteen minutes, and only has to happen once.")

    if platform.system() == "Windows":
        bash = os.path.join(_msys2_root(), "usr", "bin", "bash.exe")
        # No `cd`, no `&&`, no inline JOBS=: invoke runs this through cmd.exe, which
        # splits on && before bash ever sees the line (single quotes do not protect in
        # cmd). Everything those provided comes through the environment instead:
        # MSYSTEM=UCRT64 with a login shell puts the ucrt64 toolchain on PATH,
        # CHERE_INVOKING keeps the login shell in the invoking directory (we are in
        # base_folder courtesy of the chdir in `setup`), and build_ipopt.sh reads JOBS.
        env = {"MSYSTEM": "UCRT64", "CHERE_INVOKING": "1", "JOBS": str(jobs), "VERBOSITY": _verbosity()}
        if not in_tree:
            # WORK_DIR is consumed by bash, so it has to be a path bash can hold
            env["WORK_DIR"] = _msys_path(work_dir)
        ctx.run('"{0}" -lc "packaging/build_ipopt.sh"'.format(bash), env=env)
    else:
        env = {"JOBS": str(jobs), "VERBOSITY": _verbosity()}
        if not in_tree:
            env["WORK_DIR"] = work_dir
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
    macOS and Linux: install the toolchain - Homebrew included, on a Mac that has
    none - stage
    IPOPT if it is not staged already, then install the package editable with the
    compilers and link directories the extension needs.

    Idempotent - a second run reinstalls the package and leaves the IPOPT tree alone.
    """
    _refuse_root()

    with chdir(ctx.base_folder):
        if toolchain:
            _install_toolchain(ctx)

        _build_ipopt(ctx, jobs)

        _patch_venv_activation()

        # uv where the contributor is using it, pip otherwise; both install into whatever
        # environment is active
        pip = "uv pip" if shutil.which("uv") else "python -m pip"
        env = _native_build_env()
        stage = _ipopt_prefix_override(ctx.base_folder)
        if stage:
            # the stage tree is not at build/ipopt/stage, so CMakeLists.txt has to be told
            env["IPOPT_PREFIX"] = stage
        ctx.run('%s install -e ".[dev]"' % pip, env=env)

        _codesign_extension(ctx)

    print("\nDone. `invoke test` should now pass.")


@invoke.task(
    help={
        "serve": (
            "Serve at localhost with live reload (the default). --no-serve builds "
            "into dist/docs instead - the artifact CI deploys."
        ),
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
