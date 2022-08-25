"""Example to simulate wedge example type-a, type-b, type-c, type-d"""

import os
import math as mt
import compas
import compas_cra

from compas_cra.datastructures import CRA_Assembly
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view


mu = 0.84  # friction coefficient
deg = 90  # rotation in degree
axis = "y-axis"  # y-axis, x-axis, xy30-axis

FILE_I = os.path.join(compas_cra.SAMPLE, "type-b.json")

rotate_axis = [0, 1, 0]
if axis == "y-axis":
    rotate_axis = [0, 1, 0]  # around y-axis
if axis == "x-axis":
    rotate_axis = [1, 0, 0]  # around x-axis
if axis == "xy30-axis":
    rotate_axis = [mt.sqrt(3), 1, 0]  # rotate around xy30-axis

assembly = compas.json_load(FILE_I)
assembly = assembly.copy(cls=CRA_Assembly)
assembly.set_boundary_conditions([0, 1])

assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)

assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)

cra_solve(assembly, verbose=True, timer=True, d_bnd=1e-2)
cra_view(
    assembly,
    resultant=False,
    nodal=True,
    grid=True,
    displacements=False,
    dispscale=0,
    scale=0.5,
)
