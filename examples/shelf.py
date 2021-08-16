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
            os.path.join(compas_cra.DATA, './shelf-stable.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])

    assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

    print("blocks: ", assembly.number_of_nodes())
    print("interfaces: ", assembly.number_of_edges())

    dispbnd = 1e-2
    overlap = 1e-4
    d = 1
    # cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap)
    cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap)
    cra_view(assembly, resultant=True, nodal=True, grid=True,
             displacements=True, dispscale=1, scale=5)
