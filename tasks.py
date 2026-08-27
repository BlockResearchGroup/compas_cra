from __future__ import print_function

import os

import invoke
from compas_invocations2 import build
from compas_invocations2 import style
from compas_invocations2 import tests
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
    build.clean(ctx)

    if confirm(
        "Everything is ready. You are about to push to git, which will build every wheel "
        "and release to pypi.org. Are you sure?",
        assume_yes=False,
    ):
        ctx.run("git push --tags && git push")
    else:
        raise invoke.Exit("You need to manually revert the tag/commits created.")


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
    release,
)
ns.configure({"base_folder": os.path.dirname(__file__)})
