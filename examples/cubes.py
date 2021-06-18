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
    from compas_cra.equilibrium import cra_solve
    from compas_cra.viewers import cra_view

    # assembly = CRA_Assembly.from_json(
    #     os.path.join(compas_cra.DATA, './cubes.json'))
    assembly = compas.json_load(
            os.path.join(compas_cra.DATA, './cube.json'))
    # assembly = compas.json_load(
    #     os.path.join(compas_cra.DATA, './armadillo.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])

    assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)
    # assembly_interfaces_numpy(assembly, nmax=10, tmax=0.05, amin=0.0001)

    print("blocks: ", assembly.number_of_nodes())
    print("interfaces: ", assembly.number_of_edges())

    # compas.json_dump(assembly,
    #                  os.path.join(compas_cra.DATA, './cubes-test.json'))

    cra_solve(assembly, verbose=True, timer=True)
    cra_view(assembly, resultant=False, nodal=True, grid=True,
             displacements=True, dispscale=0, scale=0.5)
