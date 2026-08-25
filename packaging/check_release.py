"""Refuse to publish a release that is missing a platform, or a solver.

The wheels of this project are only useful because IPOPT is linked into the
`compas_cra._native._core` extension. A pure `py3-none-any` wheel installs cleanly and
then fails at the first solve, and a release missing one platform or one CPython
silently leaves those users on the previous version. Both are easy to cause by
accident - a renamed job, a dropped artifact, a local `python -m build` uploaded by
hand - and neither is visible in a green CI run.

So the release is checked before it is uploaded, not after.

Usage:  python packaging/check_release.py dist
"""

import sys
import zipfile
from pathlib import Path

# Every platform we promise a working solver on, matched by substring so that bumping a
# macOS deployment target or picking up an extra manylinux alias does not fail the
# release. The compound linux tags really do look like
# `manylinux_2_27_x86_64.manylinux_2_28_x86_64`.
EXPECTED_PLATFORMS = {
    "linux x86_64": "manylinux_2_28_x86_64",
    "linux aarch64": "manylinux_2_28_aarch64",
    "macos arm64": "_arm64",
    "macos x86_64": "macosx_",  # refined below, see `matches`
    "windows x86_64": "win_amd64",
}

# Kept in step with `build` in [tool.cibuildwheel]; Rhino 8 is the reason cp39 is here.
EXPECTED_PYTHONS = ("cp39", "cp310", "cp311", "cp312", "cp313")

# The extension carries a statically linked IPOPT and MUMPS; it weighs 4-6 MB. Anything
# tiny is a stub, or a build that failed to link the solver in.
MIN_EXTENSION_BYTES = 1_000_000

EXTENSION_PREFIX = "compas_cra/_native/_core"


def matches(platform, tag):
    """Return True if the platform tag `tag` belongs to `platform`."""
    if platform == "macos arm64":
        return tag.startswith("macosx_") and tag.endswith("_arm64")
    if platform == "macos x86_64":
        return tag.startswith("macosx_") and tag.endswith("_x86_64")
    return EXPECTED_PLATFORMS[platform] in tag


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

    # 3. every platform accounted for, on every CPython we claim to support
    tags = {w.name.rsplit("-", 1)[-1][: -len(".whl")] for w in wheels}
    for platform in sorted(EXPECTED_PLATFORMS):
        found = {t for t in tags if matches(platform, t)}
        if not found:
            problems.append("no wheel for %s" % platform)
            continue
        pythons = {w.name.split("-")[2] for w in wheels if w.name.rsplit("-", 1)[-1][: -len(".whl")] in found}
        missing = [p for p in EXPECTED_PYTHONS if p not in pythons]
        if missing:
            problems.append("%s is missing %s" % (platform, ", ".join(missing)))

    unclaimed = sorted(t for t in tags if not any(matches(p, t) for p in EXPECTED_PLATFORMS) and t != "any")
    if unclaimed:
        problems.append("unexpected platform tag(s): %s" % ", ".join(unclaimed))

    # 4. an sdist, which is what unsupported platforms fall back to
    if not sdists:
        problems.append("no sdist")

    # 5. and the point of the whole exercise: each wheel actually carries the extension
    for wheel in wheels:
        with zipfile.ZipFile(wheel) as zf:
            extensions = [i for i in zf.infolist() if i.filename.startswith(EXTENSION_PREFIX)]
            if not extensions:
                problems.append("%s contains no %s* extension" % (wheel.name, EXTENSION_PREFIX))
            elif extensions[0].file_size < MIN_EXTENSION_BYTES:
                problems.append(
                    "%s has a suspiciously small extension (%d bytes)" % (wheel.name, extensions[0].file_size)
                )
            else:
                print("  ok  %-62s _core %5.1f MB" % (wheel.name, extensions[0].file_size / 1e6))

    for sdist in sdists:
        print("  ok  %-62s (no extension, by design)" % sdist.name)

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
    print(
        "\nall %d platforms present on %s, every wheel carries a solver"
        % (len(EXPECTED_PLATFORMS), ", ".join(EXPECTED_PYTHONS))
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
