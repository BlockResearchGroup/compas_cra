"""This is the script to add assembly without interfaces from rhino meshes"""

import os
import compas
import rhinoscriptsyntax as rs

from compas_rhino import select_meshes
from compas_cra.datastructures import CRA_Assembly

HERE = os.path.abspath(os.path.dirname(__file__))
PATH = os.path.abspath(
    os.path.join(HERE, "..", "data")
)  # or change to your own directory

guid = select_meshes()

assembly = CRA_Assembly()
assembly.add_blocks_from_rhinomeshes(guid)

filename = rs.GetString("file name (xxx.json):")
file_o = os.path.join(PATH, filename)
compas.json_dump(assembly, file_o)
print("file save to: ", file_o)
