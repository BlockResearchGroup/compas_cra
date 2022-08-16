#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Parametric Arch
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':

    from compas_cra.viewers import cra_view
    from compas_cra.equilibrium import cra_penalty_solve, cra_solve
    from compas_cra.datastructures import assembly_interfaces_numpy
    from compas_cra.geometry import Arch

    height = 5
    span = 10
    thickness = 0.5
    depth = 0.5
    num_blocks = 20

    assembly = Arch(height=height, span=span, thickness=thickness, depth=depth, num_blocks=num_blocks, extra_support=False).assembly()

    assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)

    # cra_penalty_solve(assembly, mu=0.7, verbose=True, timer=True)
    cra_solve(assembly, mu=0.7, verbose=True, timer=True)
    cra_view(assembly, resultant=True, nodal=False, grid=True, forcesdirect=False, forcesline=True, displacements=True, dispscale=10, scale=0.5)
