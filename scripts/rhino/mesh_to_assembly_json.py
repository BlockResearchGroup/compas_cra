#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is the script to add assembly without interfaces from rhino meshes
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':
    import os
    import rhinoscriptsyntax as rs

    from compas_rhino import select_meshes
    from compas_assembly.datastructures import Assembly

    HERE = os.path.abspath(os.path.dirname(__file__))
    DATA = os.path.abspath(os.path.join(HERE, "../../data/"))

    guid = select_meshes()

    assembly = Assembly()
    assembly.add_blocks_from_rhinomeshes(guid)

    filename = rs.GetString("file name (xxx.json):")
    print(os.path.join(DATA + '/' + filename))
    assembly.to_json(os.path.join(DATA + '/' + filename))
