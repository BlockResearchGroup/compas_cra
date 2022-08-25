"""This is the script to add assembly and interfaces from rhino meshes"""

import os
import compas
import rhinoscriptsyntax as rs

from compas_rhino.geometry import RhinoMesh
from compas_cra.datastructures import CRA_Assembly

HERE = os.path.abspath(os.path.dirname(__file__))
PATH = os.path.abspath(
    os.path.join(HERE, "..", "data")
)  # or change to your own directory

guid = rs.GetObjects(
    "select blocks",
    preselect=False,
    select=False,
    group=False,
    filter=rs.filter.mesh,
)

assembly = CRA_Assembly()
assembly.add_blocks_from_rhinomeshes(guids=guid)

node_labels = []
for node in assembly.nodes():
    block = assembly.graph.node_attribute(node, "block")
    c = block.centroid()
    node_labels.append(rs.AddTextDot(node, c))

IS_FINISHED = False
while not IS_FINISHED:

    interface_guids = rs.GetObjects(
        "select interfaces",
        preselect=False,
        select=False,
        group=False,
        filter=rs.filter.mesh,
    )

    interfaces = []
    for interface_guid in interface_guids:
        mesh = RhinoMesh.from_guid(interface_guid)
        interfaces.append(mesh.to_compas())

    edge_a = rs.GetInteger("assign interface from")
    edge_b = rs.GetInteger("assign interface to")

    assembly.add_interfaces_from_meshes(interfaces, edge_a, edge_b)

    IS_FINISHED = rs.GetBoolean(
        "Continue select interface?", ("Continue", "Continue", "Stop"), (False)
    )[0]

rs.DeleteObjects(node_labels)

filename = rs.GetString("file name (xxx.json):")
file_o = os.path.join(PATH, filename)
compas.json_dump(assembly, file_o)
print("file save to: ", file_o)
