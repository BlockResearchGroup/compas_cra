"""Simple example to calculate cubes"""

from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import rbe_solve
from compas_cra.viewers import cra_view

deg = 48.00  # rotation in degree
rotate_axis = [0, 1, 0]  # around y-axis

support = Box(1, 1, 1)  # supporting block
free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])))  # block to analyse

assembly = CRA_Assembly()
assembly.add_block(Block.from_shape(support))
assembly.add_block(Block.from_shape(free1))
assembly.set_boundary_conditions([0])

interface1 = Mesh()
# interface corners
corners = [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5]]
for i, c in enumerate(corners):
    interface1.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
interface1.add_face([0, 1, 2, 3])

assembly.add_interfaces_from_meshes([interface1], 0, 1)

assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)  # around y-axis

rbe_solve(assembly, mu=0.1, verbose=True, timer=True)
cra_view(assembly, resultant=False, nodal=True, grid=True, displacements=True, dispscale=10)
