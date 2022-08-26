"""This is the complete code for the tutorial"""

from compas.geometry import Box, Frame, Translation, Rotation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view


support = Box(Frame.worldXY(), 4, 2, 1)  # supporting block
free1 = Box(
    Frame.worldXY().transformed(
        Translation.from_vector([0, 0, 1])
        * Rotation.from_axis_and_angle([0, 0, 1], 0.2)
    ),
    1,
    3,
    1,
)  # block to analyse

assembly = CRA_Assembly()  # create empty assembly

assembly.add_block(Block.from_shape(support), node=0)  # add support to assembly
assembly.add_block(Block.from_shape(free1), node=1)  # add block to assembly

assembly.set_boundary_conditions([0])  # set support as boundary condition

assembly_interfaces_numpy(assembly)  # identify interface between two blocks

cra_solve(assembly, verbose=True, timer=True)  # solve equilibrium using cra_solve

cra_view(assembly, resultant=False, nodal=True, grid=True)  # visiualisation
