#! python3
# venv: compas-cra
# r: compas_cra>=0.7.2
"""Check that the in-process IPOPT solver is available in this Python environment."""

import os

import compas_cra

print("compas_cra:", compas_cra.__version__)
print("package at:", os.path.dirname(compas_cra.__file__))

try:
    from compas_cra import _native

    print("solver: OK, IPOPT", _native.IPOPT_VERSION)
except ImportError as e:
    print("solver: NOT AVAILABLE ({})".format(e))
    print("fix: pip install --force-reinstall compas_cra  (the solver ships inside the platform wheel)")
