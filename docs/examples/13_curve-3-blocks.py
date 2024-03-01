"""Example to calculate three block with curved interfaces"""

import os

import compas
import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

density = 1

FILE_I = os.path.join(compas_cra.SAMPLE, "curve-3-blocks.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

cra_solve(assembly, verbose=True, timer=True, density=density)
cra_view(
    assembly,
    resultant=False,
    nodal=True,
    grid=True,
    displacements=True,
    dispscale=0,
    scale=1,
    density=density,
)
