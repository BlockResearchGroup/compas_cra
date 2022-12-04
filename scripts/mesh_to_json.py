"""This is the script to add assembly without interfaces from rhino meshes"""

import os
import compas
import compas_rhino.geometry

import rhinoscriptsyntax as rs

from compas_rhino import select_mesh

HERE = os.path.abspath(os.path.dirname(__file__))
PATH = os.path.abspath(
    os.path.join(HERE, "..", "data")
)  # or change to your own directory

guid = select_mesh()

mesh = compas_rhino.geometry.RhinoMesh.from_guid(guid).to_compas()

filename = rs.GetString("file name (xxx.json):")
file_o = os.path.join(PATH, filename)
compas.json_dump(mesh, file_o)
print("file save to: ", file_o)
