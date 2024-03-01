"""Simple example to calculate three stacked cubes"""

from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve
from compas_cra.viewers import cra_view

deg = 20  # rotation angle in degree
rotate_axis = [0, 1, 0]  # around y-axis

support = Box(1, 1, 1)
free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])))
free2 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 2])))

assembly = CRA_Assembly()
assembly.add_block(Block.from_shape(support), node=0)
assembly.add_block(Block.from_shape(free1), node=1)
assembly.add_block(Block.from_shape(free2), node=2)
assembly.set_boundary_conditions([0])

interface1 = Mesh()
corners = [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5]]
for i, c in enumerate(corners):
    interface1.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
interface1.add_face([0, 1, 2, 3])

interface2 = Mesh()
corners = [[0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, -0.5, 1.5], [0.5, -0.5, 1.5]]
for i, c in enumerate(corners):
    interface2.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
interface2.add_face([0, 1, 2, 3])

assembly.add_interfaces_from_meshes([interface1], 0, 1)
assembly.add_interfaces_from_meshes([interface2], 1, 2)

assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)

cra_solve(assembly, verbose=True, timer=True)
cra_view(assembly, resultant=False, nodal=True, grid=True, displacements=True, dispscale=10)
