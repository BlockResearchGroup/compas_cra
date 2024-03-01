from __future__ import print_function

import os

from compas_invocations2 import build
from compas_invocations2 import docs
from compas_invocations2 import style
from compas_invocations2 import tests
from invoke import Collection

ns = Collection(
    docs.help,
    style.check,
    style.lint,
    style.format,
    docs.docs,
    docs.linkcheck,
    tests.test,
    tests.testdocs,
    build.build_ghuser_components,
    build.prepare_changelog,
    build.clean,
    build.release,
)
ns.configure(
    {
        "base_folder": os.path.dirname(__file__),
        "ghuser": {
            "source_dir": "src/compas_notebook/ghpython/components",
            "target_dir": "src/compas_notebook/ghpython/components/ghuser",
        },
    }
)
