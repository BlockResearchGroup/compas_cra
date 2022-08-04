#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple example to calculate three stacked cubes
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':

    import compas
    import compas_cra
    import os

    from compas_cra.datastructures import CRA_Assembly
    from compas_cra.datastructures import assembly_interfaces_numpy
    from compas_cra.equilibrium import cra_solve, cra_penalty_solve
    from compas_cra.viewers import cra_view

    assembly = compas.json_load(
            os.path.join(compas_cra.DATA, './snake.json'))
    assembly = assembly.copy(cls=CRA_Assembly)

    assembly.set_boundary_conditions([0])
    # assembly.delete_node(1)
    # assembly.delete_node(2)
    # assembly.delete_node(3)

    assembly_interfaces_numpy(assembly, amin=1e-4, tmax=1e-2)

    print(assembly)

    mu = 0.7
    dispbnd = 1e-1
    overlap = 1e-4 * 0
    d = 1
    cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
    # cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
    cra_view(assembly, resultant=True, nodal=True, grid=True, weights=False,
             displacements=True, dispscale=1, scale=.1/d)