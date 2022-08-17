#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example to calculate interlocking joint
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':

    import compas
    import compas_cra
    import os

    from compas_cra.datastructures import CRA_Assembly
    from compas_cra.equilibrium import cra_solve
    from compas_cra.viewers import cra_view

    mu = 0.84
    density = 1
    deg = 50  # rotation angle in degree
    rotate_axis = [0, 1, 0]  # around y-axis

    assembly = compas.json_load(os.path.join(compas_cra.DATA, './concave-long.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])

    assembly.rotate_assembly([0, 0, 0], rotate_axis, deg)

    cra_solve(assembly, verbose=True, timer=True, density=density, mu=mu, eps=1e-3, d_bnd=1e-2)
    cra_view(assembly, resultant=True, nodal=True, grid=True,
             displacements=True, dispscale=1, scale=5, density=density)
