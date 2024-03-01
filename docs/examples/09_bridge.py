"""Simple example to calculate three stacked cubes"""

import os

import compas
import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_penalty_solve
from compas_cra.equilibrium import cra_solve
from compas_cra.equilibrium import density_setup
from compas_cra.viewers import cra_view

mu = 0.9
dispbnd = 1e-1
overlap = 0
d = 1

FILE_I = os.path.join(compas_cra.SAMPLE, "bridge.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0, 1])

density = {node: 3.51 if node in [11, 12, 13, 14, 15] else 1 for node in assembly.graph.nodes()}
density_setup(assembly, density)

assembly.delete_blocks([3, 4, 5, 6, 7, 8, 9, 10])
assembly.delete_blocks([11, 12, 13, 14, 15])

assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
# cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
cra_view(
    assembly,
    resultant=True,
    nodal=False,
    grid=True,
    weights=True,
    forcesdirect=False,
    forcesline=True,
    displacements=True,
    dispscale=1,
    scale=0.5 / d,
)
