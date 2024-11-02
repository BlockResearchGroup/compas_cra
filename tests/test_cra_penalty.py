from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.equilibrium import cra_penalty_solve


def test_cra_penalty():
    support = Box(1, 1, 1)  # supporting block
    free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0.75, 0, 1])))  # block to analyse

    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), node=0)
    assembly.add_block(Block.from_shape(free1), node=1)
    assembly.set_boundary_conditions([0])

    assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

    cra_penalty_solve(assembly, density=1, gravity=1)

    block = assembly.graph.node_attribute(1, "block")
    weight = block.volume()
    resultant = 0
    for edge in assembly.graph.edges():
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            corners = interface.points
            forces = interface.forces
            for i, corner in enumerate(corners):
                force = forces[i]["c_np"] - forces[i]["c_nn"]
                resultant += force

    assert round(weight, 2) == round(resultant, 2)
