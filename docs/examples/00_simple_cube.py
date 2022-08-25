"""Simple example to calculate cubes"""

from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

support = Box(Frame.worldXY(), 1, 1, 1)  # supporting block
free1 = Box(
    Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])), 1, 1, 1
)  # block to analyse

assembly = CRA_Assembly()
assembly.add_block(Block.from_shape(support), node=0)
assembly.add_block(Block.from_shape(free1), node=1)
assembly.set_boundary_conditions([0])

interface1 = Mesh()
# interface corners
corners = [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5]]
for i, c in enumerate(corners):
    interface1.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
interface1.add_face([0, 1, 2, 3])

assembly.add_interfaces_from_meshes([interface1], 0, 1)

cra_solve(assembly, verbose=True, timer=True, density=1)
cra_view(assembly, resultant=False, nodal=True, grid=True, forcesline=False, scale=1)
