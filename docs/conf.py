# -*- coding: utf-8 -*-

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = "1.0"

import sys
import os
import inspect
import importlib

import sphinx_compas_theme
from sphinx.ext.napoleon.docstring import NumpyDocstring

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

# -- General configuration ------------------------------------------------

project = "COMPAS CRA"
copyright = "Gene Ting-Chun Kao, Block Research Group"
author = "Gene Ting-Chun Kao"
release = "0.1.0"
version = ".".join(release.split(".")[0:2])

master_doc = "index"
source_suffix = [
    ".rst",
]
templates_path = sphinx_compas_theme.get_autosummary_templates_path()
exclude_patterns = []

pygments_style = "sphinx"
show_authors = True
add_module_names = True
language = None


# -- Extension configuration ------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.linkcode",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "matplotlib.sphinxext.plot_directive",
]

# autodoc options

autodoc_type_aliases = {}

# this does not work properly yet
autodoc_typehints = "none"
autodoc_typehints_format = "short"
autodoc_typehints_description_target = "documented"

autodoc_default_options = {
    "undoc-members": True,
    "show-inheritance": True,
}

autodoc_member_order = "alphabetical"

autoclass_content = "class"


def skip(app, what, name, obj, would_skip, options):
    if name.startswith("_"):
        return True
    return would_skip


def setup(app):
    app.connect("autodoc-skip-member", skip)


# autosummary options

autosummary_generate = True

# napoleon options

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = False
napoleon_use_rtype = False

# plot options

plot_html_show_source_link = False
plot_html_show_formats = False

# docstring sections


def parse_attributes_section(self, section):
    return self._format_fields("Attributes", self._consume_fields())


NumpyDocstring._parse_attributes_section = parse_attributes_section


def patched_parse(self):
    self._sections["attributes"] = self._parse_attributes_section
    self._unpatched_parse()


NumpyDocstring._unpatched_parse = NumpyDocstring._parse
NumpyDocstring._parse = patched_parse

# intersphinx options

intersphinx_mapping = {
    "python": ("https://docs.python.org/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "compas": ("https://compas.dev/compas/latest/", None),
    "compas_assembly": (
        "https://blockresearchgroup.github.io/compas_assembly/0.4.1/",
        None,
    ),
}

# linkcode


def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None
    if not info["fullname"]:
        return None

    package = info["module"].split(".")[0]
    if not package.startswith("compas_cra"):
        return None

    module = importlib.import_module(info["module"])
    parts = info["fullname"].split(".")

    if len(parts) == 1:
        obj = getattr(module, info["fullname"])
        filename = inspect.getmodule(obj).__name__.replace(".", "/")
        lineno = inspect.getsourcelines(obj)[1]
    elif len(parts) == 2:
        obj_name, attr_name = parts
        obj = getattr(module, obj_name)
        attr = getattr(obj, attr_name)
        if inspect.isfunction(attr):
            filename = inspect.getmodule(obj).__name__.replace(".", "/")
            lineno = inspect.getsourcelines(attr)[1]
        else:
            return None
    else:
        return None

    return f"https://github.com/BlockResearchGroup/compas_cra/blob/main/src/{filename}.py#L{lineno}"


# extlinks

extlinks = {}

# -- Options for HTML output ----------------------------------------------

html_theme = "compaspkg"
html_theme_path = sphinx_compas_theme.get_html_theme_path()

html_theme_options = {
    "package_name": "compas_cra",
    "package_title": project,
    "package_version": release,
    "package_author": "Gene Ting-Chun Kao",
    "package_docs": "https://BlockResearchGroup.github.io/compas_cra/",
    "package_repo": "https://github.com/BlockResearchGroup/compas_cra",
    "package_old_versions_txt": "https://BlockResearchGroup.github.io/compas_cra/doc_versions.txt",
}

html_context = {}
html_static_path = sphinx_compas_theme.get_html_static_path()
html_extra_path = []
html_last_updated_fmt = ""
html_copy_source = False
html_show_sourcelink = False
html_permalinks = False
html_permalinks_icon = ""
html_experimental_html5_writer = False
html_compact_lists = True