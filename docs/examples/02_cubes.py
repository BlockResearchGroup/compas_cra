"""Simple example to calculate three stacked cubes"""

import os

import compas
import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

FILE_I = os.path.join(compas_cra.SAMPLE, "cubes.json")

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0])

assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)

cra_solve(assembly, verbose=True, timer=True)
cra_view(
    assembly,
    resultant=False,
    nodal=True,
    grid=True,
    displacements=True,
    dispscale=0,
    scale=0.5,
)
