"""Example to calculate cube tall with curved interfaces"""

import os

import compas

import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

density = 0.1

FILE_I = os.path.join(compas_cra.SAMPLE, "cube-curve-tall.json")

assembly = compas.json_load(FILE_I)
assembly: CRA_Assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

cra_solve(assembly, verbose=True, timer=True, density=density)
cra_view(
    assembly,
    resultant=True,
    nodal=False,
    grid=True,
    displacements=True,
    dispscale=0,
    scale=40,
    density=density,
)
