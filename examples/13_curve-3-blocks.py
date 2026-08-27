"""Example to calculate three block with curved interfaces"""

import os

import compas

import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

# The published screenshot's original parameters. They only solve in the barrier-
# interior regime the solver now runs in (mu_target - see cra_native.py); under strict
# optimization they exhaust the iteration cap on every IPOPT generation, which is why
# the author had switched to density=1 as a workaround days after the screenshot.
density = 0.1

FILE_I = os.path.join(compas_cra.SAMPLE, "curve-3-blocks.json")

assembly = compas.json_load(FILE_I)
assembly: CRA_Assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

cra_solve(assembly, verbose=True, timer=True, density=density)
cra_view(
    assembly,
    resultant=False,
    nodal=True,
    grid=True,
    displacements=True,
    dispscale=0,
    scale=5,
    density=density,
)
