"""Example to calculate interlocking joint"""

import os

import compas

import compas_cra
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

mu = 0.84
density = 1
deg = 50  # rotation angle in degree
rotate_axis = [0, 1, 0]  # around y-axis

FILE_I = os.path.join(compas_cra.SAMPLE, "concave-long.json")

assembly = compas.json_load(FILE_I)
assembly: CRA_Assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)

cra_solve(assembly, verbose=True, timer=True, density=density, mu=mu, eps=1e-3, d_bnd=1e-2)
cra_view(
    assembly,
    resultant=True,
    nodal=True,
    grid=True,
    displacements=True,
    dispscale=1,
    scale=5,
    density=density,
)
