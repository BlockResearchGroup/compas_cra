#! python3
# venv: compas-cra
# r: compas_cra>=0.7.2
"""Tall cube stack with curved interfaces: CRA equilibrium, drawn into the Rhino document.

Rhino port of docs/examples/15_cube-curve-tall.py. Keep rhino_cra_view.py next to this
script (and save the script) so the shared drawing module can be imported.
"""

import os
import sys

import compas

import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:  # unsaved ScriptEditor document: assume the cwd holds the helper
    HERE = os.getcwd()
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from rhino_cra_view import cra_view_rhino  # noqa: E402

density = 0.1

FILE_I = os.path.join(compas_cra.SAMPLE, "cube-curve-tall.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

cra_solve(assembly, verbose=True, timer=True, density=density)
cra_view_rhino(
    assembly,
    resultant=True,
    nodal=False,
    scale=40,
    density=density,
    layer_prefix="CRA-CubeCurveTall",
)
