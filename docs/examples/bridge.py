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
    from compas_cra.equilibrium import density_setup
    from compas_cra.viewers import cra_view

    mu = 0.9
    dispbnd = 1e-1
    overlap = 0
    d = 1

    assembly = compas.json_load(
            os.path.join(compas_cra.DATA, './bridge.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0, 1])

    density = {node: 3.51 if node in [11, 12, 13, 14, 15] else 1 for node in assembly.graph.nodes()}
    # for node in assembly.graph.nodes():
    #     block = assembly.graph.node_attribute(node, 'block')
    #     block.attributes["density"] = 1
    #     if node in [11, 12, 13, 14, 15]:
    #         block.attributes["density"] = 3.51
    density_setup(assembly, density)
    # assembly.graph.delete_node(2)
    # assembly.graph.delete_node(3)
    # assembly.graph.delete_node(4)
    # assembly.graph.delete_node(5)
    # assembly.graph.delete_node(6)
    # assembly.graph.delete_node(7)
    # assembly.graph.delete_node(8)
    # assembly.graph.delete_node(9)
    # assembly.graph.delete_node(10)

    # assembly.graph.delete_node(11)
    # assembly.graph.delete_node(12)
    # assembly.graph.delete_node(13)
    # assembly.graph.delete_node(14)
    # assembly.graph.delete_node(15)

    assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)

    cra_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
    # cra_penalty_solve(assembly, verbose=True, density=d, d_bnd=dispbnd, eps=overlap, mu=mu)
    cra_view(assembly, resultant=True, nodal=False, grid=True, weights=True, forcesdirect=False, forcesline=True,
             displacements=True, dispscale=1, scale=.5/d)
