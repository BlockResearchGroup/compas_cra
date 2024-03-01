"""Simple example to calculate three stacked cubes"""

import os

import compas
import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_penalty_solve
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

mu = 0.7
dispbnd = 1e-1
overlap = 1e-4 * 0
d = 1

FILE_I = os.path.join(compas_cra.SAMPLE, "snake.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)

assembly.set_boundary_conditions([0])
# assembly.delete_node(1)
# assembly.delete_node(2)
# assembly.delete_node(3)

assembly_interfaces_numpy(assembly, amin=1e-4, tmax=1e-2)

cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
# cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
cra_view(
    assembly,
    resultant=True,
    nodal=True,
    grid=True,
    weights=False,
    displacements=True,
    dispscale=1,
    scale=0.1 / d,
)
