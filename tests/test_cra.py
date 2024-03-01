from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_solve


def test_cra():
    support = Box(1, 1, 1)  # supporting block
    free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])))  # block to analyse

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

    cra_solve(assembly, density=1)

    IS_FORCE_CORRECT = True
    for edge in assembly.graph.edges():
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            corners = interface.points
            forces = interface.forces
            for i, corner in enumerate(corners):
                force = forces[i]["c_np"] - forces[i]["c_nn"]
                if round(force, 2) == 0.25:
                    IS_FORCE_CORRECT = True
                else:
                    IS_FORCE_CORRECT = False

    assert IS_FORCE_CORRECT
