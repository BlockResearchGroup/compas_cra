"""Model H and A."""

import os

import compas
import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import rbe_solve
from compas_cra.viewers import cra_view

IS_A = True  # True: model A. False: model H

mu = 0.84  # friction coefficient
deg = 0  # rotation in degree

assembly = compas.json_load(os.path.join(compas_cra.SAMPLE, "A.json" if IS_A else "H.json"))
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0, 1])

assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)

rbe_solve(assembly, verbose=True, timer=True)
cra_view(
    assembly,
    resultant=False,
    nodal=True,
    grid=True,
    displacements=False,
    dispscale=0,
    scale=1,
)
