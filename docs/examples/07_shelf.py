"""Simple example to calculate three stacked cubes"""

import os
import compas
import compas_cra

from compas_cra.datastructures import CRA_Assembly
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.equilibrium import cra_solve, cra_penalty_solve
from compas_cra.viewers import cra_view

mu = 0.9
dispbnd = 1e-1
overlap = 1e-3
d = 1

FILE_I = os.path.join(compas_cra.SAMPLE, "shelf.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])
# assembly.set_boundary_conditions([i+1 for i in range(9)])
# assembly.graph.delete_node(11)
# assembly.graph.delete_node(12)

assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

# cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
cra_view(
    assembly,
    resultant=True,
    nodal=False,
    grid=False,
    weights=False,
    displacements=False,
    dispscale=1,
    scale=10 / d,
)
