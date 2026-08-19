"""Refuse to publish a release that is missing a platform, or a solver.

The wheels of this project are only useful because they carry an ipopt executable. A
pure `py3-none-any` wheel installs cleanly and then fails at the first solve, and a
release missing one platform silently leaves those users on the previous version. Both
are easy to cause by accident - a renamed job, a dropped artifact, a local `python -m
build` uploaded by hand - and neither is visible in a green CI run.

So the release is checked before it is uploaded, not after.

Usage:  python packaging/check_release.py dist
"""

import sys
import zipfile
from pathlib import Path

# Every platform we promise a working solver on.
EXPECTED_PLATFORMS = {
    "manylinux_2_28_x86_64",
    "manylinux_2_28_aarch64",
    "macosx_11_0_arm64",
    "macosx_11_0_x86_64",
    "win_amd64",
}

# A stripped ipopt is several MB; anything tiny is a stub or a broken copy.
MIN_BINARY_BYTES = 1_000_000


def check(dist):
    """Return a list of problems with the distributions in `dist`."""
    problems = []
    wheels = sorted(dist.glob("*.whl"))
    sdists = sorted(dist.glob("*.tar.gz"))

    if not wheels:
        return ["no wheels at all in %s" % dist]

    # 1. a pure wheel would be installed in preference to the sdist on any platform we
    #    do not build for, and would carry no solver
    pure = [w.name for w in wheels if w.name.endswith("-py3-none-any.whl")]
    if pure:
        problems.append("pure py3-none-any wheel(s) present, these carry no solver: %s" % ", ".join(pure))

    # 2. exactly one version, so a stale build cannot ride along
    versions = {w.name.split("-")[1] for w in wheels} | {s.name[: -len(".tar.gz")].split("-")[-1] for s in sdists}
    if len(versions) > 1:
        problems.append("more than one version in dist: %s" % ", ".join(sorted(versions)))

    # 3. every platform accounted for
    found = {w.name.rsplit("-", 1)[-1][: -len(".whl")] for w in wheels}
    missing = EXPECTED_PLATFORMS - found
    if missing:
        problems.append("no wheel for: %s" % ", ".join(sorted(missing)))
    unexpected = found - EXPECTED_PLATFORMS - {"any"}
    if unexpected:
        problems.append("unexpected platform tag(s): %s" % ", ".join(sorted(unexpected)))

    # 4. an sdist, which is what unsupported platforms fall back to
    if not sdists:
        problems.append("no sdist")

    # 5. and the point of the whole exercise: each wheel actually contains a solver
    for wheel in wheels:
        with zipfile.ZipFile(wheel) as zf:
            binaries = [
                i
                for i in zf.infolist()
                if i.filename in ("compas_cra/_ipopt/bin/ipopt", "compas_cra/_ipopt/bin/ipopt.exe")
            ]
            if not binaries:
                problems.append("%s contains no ipopt executable" % wheel.name)
            elif binaries[0].file_size < MIN_BINARY_BYTES:
                problems.append("%s has a suspiciously small ipopt (%d bytes)" % (wheel.name, binaries[0].file_size))
            else:
                print("  ok  %-52s ipopt %5.1f MB" % (wheel.name, binaries[0].file_size / 1e6))

    for sdist in sdists:
        print("  ok  %-52s (no binary, by design)" % sdist.name)

    return problems


def main():
    dist = Path(sys.argv[1] if len(sys.argv) > 1 else "dist")
    print("checking %s" % dist.resolve())
    problems = check(dist)
    if problems:
        print("\nthis release must not be published:")
        for p in problems:
            print("  - %s" % p)
        return 1
    print("\nall %d platforms present, every wheel carries a solver" % len(EXPECTED_PLATFORMS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
