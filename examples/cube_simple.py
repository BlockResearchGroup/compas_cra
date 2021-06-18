#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple example to calculate three stacked cubes
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

if __name__ == '__main__':

    import math as mt

    from compas.datastructures import Mesh
    from compas.geometry import Box
    from compas.geometry import Frame
    from compas.geometry import Translation
    from compas_assembly.datastructures import Block
    from compas_cra.datastructures import CRA_Assembly
    from compas_cra.equilibrium import cra_solve
    from compas_cra.viewers import cra_view

    support = Box(Frame.worldXY(), 1, 1, 1)
    free1 = Box(Frame.worldXY().transformed(
        Translation.from_vector([0, 0, 1])), 1, 1, 1)
    free2 = Box(Frame.worldXY().transformed(
        Translation.from_vector([0, 0, 2])), 1, 1, 1)

    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), key=0)
    assembly.add_block(Block.from_shape(free1), key=1)
    assembly.add_block(Block.from_shape(free2), key=2)
    assembly.set_boundary_conditions([0])

    interface1 = Mesh()
    corners = [[.5, .5, .5], [-.5, .5, .5], [-.5, -.5, .5], [.5, -.5, .5]]
    for i, c in enumerate(corners):
        interface1.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
    interface1.add_face([0, 1, 2, 3])

    interface2 = Mesh()
    corners = [[.5, .5, 1.5], [-.5, .5, 1.5], [-.5, -.5, 1.5], [.5, -.5, 1.5]]
    for i, c in enumerate(corners):
        interface2.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
    interface2.add_face([0, 1, 2, 3])

    assembly.add_interfaces_from_meshes([interface1], 0, 1)
    assembly.add_interfaces_from_meshes([interface2], 1, 2)

    deg = 20  # rotation in degree
    rad = deg * mt.pi / 180
    assembly.rotate_assembly([0, 0, 0], [0, 1, 0], rad)  # around y-axis

    cra_solve(assembly, verbose=True, timer=True)
    cra_view(assembly, resultant=False, nodal=True, grid=True,
             displacements=True, dispscale=10)
