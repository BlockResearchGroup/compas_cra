#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example to calculate three block with curved interfaces
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

    density = 0.1

    assembly = compas.json_load(os.path.join(compas_cra.DATA, './cube-curve-tall.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])

    cra_solve(assembly, verbose=True, timer=True, density=density)
    cra_view(assembly, resultant=True, nodal=False, grid=True,
             displacements=True, dispscale=0, scale=40, density=density)