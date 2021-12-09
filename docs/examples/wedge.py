#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example to simulate wedge example type-a, type-b, type-c, type-d
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':

    import compas
    import compas_cra
    import os
    import math as mt

    from compas_cra.datastructures import CRA_Assembly
    from compas_cra.datastructures import assembly_interfaces_numpy
    from compas_cra.equilibrium import cra_solve
    from compas_cra.viewers import cra_view

    assembly = compas.json_load(
            os.path.join(compas_cra.DATA, 'type-b.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0, 1])

    assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)

    print("blocks: ", assembly.number_of_nodes())
    print("interfaces: ", assembly.number_of_edges())

    mu = 0.84  # friction coefficient
    deg = 90  # rotation in degree
    rad = deg * mt.pi / 180
    # assembly.rotate_assembly([0, 0, 0], [1, 0, 0], rad)  # around x-axis
    assembly.rotate_assembly([0, 0, 0], [0, 1, 0], rad)  # around y-axis
    # assembly.rotate_assembly([0, 0, 0], [mt.sqrt(3), 1, 0], rad)
    # rotate around xy30-axis

    cra_solve(assembly, verbose=True, timer=True, d_bnd=1e-2)
    cra_view(assembly, resultant=False, nodal=True, grid=True,
             displacements=False, dispscale=0, scale=0.5)
