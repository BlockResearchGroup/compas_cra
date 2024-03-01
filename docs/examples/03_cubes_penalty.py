"""Simple example to calculate cubes"""

from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_penalty_solve
from compas_cra.viewers import cra_view

deg = 40  # rotation in degree
rotate_axis = [0, 1, 0]  # around y-axis

support = Box(1, 1, 1)  # supporting block
free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0.75, 0, 1])))  # block to analyse

assembly = CRA_Assembly()
assembly.add_block(Block.from_shape(support), node=0)
assembly.add_block(Block.from_shape(free1), node=1)
assembly.set_boundary_conditions([0])

assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)

assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

cra_penalty_solve(assembly, verbose=True, timer=True, density=1)
cra_view(assembly, resultant=False, nodal=True, grid=True, scale=1)
