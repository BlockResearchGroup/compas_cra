#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Load and visualise armadillo vault CRA result.
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


def load_result():
    assembly = compas.json_load(
        os.path.join(compas_cra.DATA, './armadillo_cra.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    return assembly


def identify_interfaces():
    assembly = compas.json_load(
        os.path.join(compas_cra.DATA, './armadillo.json'))
    assembly = assembly.copy(cls=CRA_Assembly)
    assembly_interfaces_numpy(assembly, tmax=0.05, amin=0.0001)
    return assembly


if __name__ == '__main__':

    import compas
    import compas_cra
    import os

    from compas_cra.datastructures import CRA_Assembly
    from compas_cra.datastructures import assembly_interfaces_numpy
    from compas_cra.viewers import cra_view

    # if True show precalculated result, False to identify interfaces
    load_from_result = True

    if load_from_result:
        assembly = load_result()
    else:
        assembly = identify_interfaces()

    print(assembly)

    cra_view(assembly, resultant=True, nodal=False, weights=False, grid=False,
             interfaces=not load_from_result, forcesline=True, forcesdirect=False,
             displacements=False, dispscale=0, scale=300)
