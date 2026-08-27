"""Load and visualise armadillo vault CRA result."""

import os

import compas

import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.viewers import cra_view

# if True show precalculated result, False to identify interfaces
load_from_result = True


def load_result():
    """Load armadillo vault with cra result."""
    armadillo_assembly = compas.json_load(os.path.join(compas_cra.SAMPLE, "armadillo_cra.json"))
    return armadillo_assembly


def identify_interfaces():
    """Load armadillo vault blocks."""
    armadillo_assembly = compas.json_load(os.path.join(compas_cra.SAMPLE, "armadillo.json"))
    assembly_interfaces_numpy(armadillo_assembly, tmax=0.05, amin=0.0001)
    return armadillo_assembly


if load_from_result:
    assembly = load_result()
else:
    assembly = identify_interfaces()

print(assembly)

cra_view(
    assembly,
    resultant=True,
    nodal=False,
    weights=False,
    grid=False,
    interfaces=not load_from_result,
    forcesline=True,
    forcesdirect=False,
    displacements=False,
    dispscale=0,
    scale=300,
)
