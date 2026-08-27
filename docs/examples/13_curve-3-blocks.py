"""Example to calculate three block with curved interfaces"""

import os

import compas

import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

# The published screenshot was made 2022-08-17 with density=0.1, scale=5 - parameters
# that hit IPOPT's iteration cap on MUMPS builds in the 2022 code and today alike
# (verified by running both), which is presumably why the author moved to density=1
# days later. density=1 solves; scale=0.5 draws the arrows at the same lengths the
# screenshot shows (10x the force at a twentieth of the scale).
density = 1

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
    scale=0.5,
    density=density,
)
