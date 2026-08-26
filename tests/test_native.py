"""Tests for the NLP formulations and the in-process IPOPT binding.

The reference values were produced by the historical pyomo + ipopt-executable path,
which the native solvers were validated against (machine-precision agreement on the
cube fixtures, < 1e-5 relative on the arch) before that path was removed.
"""

import os

import compas
import numpy as np
import pytest
from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import Translation
from compas_assembly.datastructures import Block

import compas_cra
from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.datastructures import CRA_Assembly
from compas_cra.equilibrium import cra_penalty_problem
from compas_cra.equilibrium import cra_problem
from compas_cra.equilibrium import cra_solve
from compas_cra.equilibrium import rbe_problem

pytest.importorskip("compas_cra._native")


def box_assembly():
    support = Box(1, 1, 1)
    free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])))
    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), node=0)
    assembly.add_block(Block.from_shape(free1), node=1)
    assembly.set_boundary_conditions([0])
    interface = Mesh()
    for i, c in enumerate([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5]]):
        interface.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
    interface.add_face([0, 1, 2, 3])
    assembly.add_interfaces_from_meshes([interface], 0, 1)
    return assembly


def overlap_assembly():
    support = Box(1, 1, 1)
    free1 = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0.75, 0, 1])))
    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), node=0)
    assembly.add_block(Block.from_shape(free1), node=1)
    assembly.set_boundary_conditions([0])
    assembly_interfaces_numpy(assembly, amin=1e-6, tmax=1e-4)
    return assembly


def interface_resultants(assembly):
    out = []
    for edge in assembly.graph.edges():
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            out.append(sum(f["c_np"] - f["c_nn"] for f in interface.forces))
    return out


@pytest.mark.parametrize(
    "builder",
    [
        lambda a: cra_problem(a, density=1),
        lambda a: cra_penalty_problem(a, density=1),
        lambda a: rbe_problem(a, density=1),
    ],
    ids=["cra", "cra_penalty", "rbe"],
)
def test_derivatives_match_finite_differences(builder):
    problem, _ = builder(overlap_assembly())
    rng = np.random.default_rng(7)
    for _ in range(3):
        x = rng.normal(scale=0.6, size=problem.n)
        lam = rng.normal(size=problem.m)
        # gradient/Hessian tolerances account for FD cancellation on the 1e6 weights
        assert problem.check_gradient(x) < 1e-2
        assert problem.check_jacobian(x) < 1e-6
        assert problem.check_hessian(x, sigma=0.9, lam=lam) < 1e-3


def test_cra_cubes_regression():
    """The cubes sample, against values validated bit-for-bit with the executable."""
    assembly = compas.json_load(os.path.join(compas_cra.SAMPLE, "cubes.json")).copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])
    assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)
    cra_solve(assembly, density=1)
    resultants = sorted(interface_resultants(assembly))
    assert resultants == pytest.approx([2.9897, 6.5960], abs=1e-3)


def test_cra_snake():
    """A 58-vertex real assembly that converges robustly on every platform."""
    assembly = compas.json_load(os.path.join(compas_cra.SAMPLE, "snake.json")).copy(cls=CRA_Assembly)
    assembly.set_boundary_conditions([0])
    assembly_interfaces_numpy(assembly, amin=1e-4, tmax=1e-2)
    cra_solve(assembly, mu=0.7, density=1, d_bnd=1e-1, eps=0)
    resultants = interface_resultants(assembly)
    assert len(resultants) > 0
    assert all(np.isfinite(r) for r in resultants)
    assert max(resultants) > 0


def arch_assembly(num_blocks=20):
    from compas_cra.geometry import Arch

    assembly = Arch(height=5.0, span=10.0, thickness=0.5, depth=0.5, num_blocks=num_blocks).assembly()
    assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)
    return assembly


def test_cra_arch():
    """The 20-block arch: a degenerate benchmark that used to be allowed to stall.

    Its convergence is steered by BLAS rounding, and with the monotone barrier update it
    took 1964 of the 3000 permitted iterations while running out of them altogether on
    the macOS wheels, which link Accelerate rather than OpenBLAS -- reported as
    ``solve failed: failed (Maximum_Iterations_Exceeded)``. The adaptive barrier update
    needs a third of the iterations, which is margin enough that the rounding no longer
    decides the outcome, so a stall is now a regression rather than a skip."""
    assembly = arch_assembly()
    cra_solve(assembly, mu=0.7)
    resultants = interface_resultants(assembly)
    assert len(resultants) == 19
    assert max(resultants) == pytest.approx(1.96, abs=0.05)
    assert min(resultants) > 0  # a standing arch is all compression


def test_cra_arch_has_iteration_headroom():
    """Guard the margin, not just the answer.

    The arch failed on macOS by exhausting IPOPT's 3000-iteration cap, and a solve that
    creeps back towards it is one rounding difference away from doing so again on some
    other wheel. Assert it finishes with room to spare: it takes ~670 iterations."""
    from compas_cra.equilibrium.cra_native import _CRA_OPTIONS
    from compas_cra.nlp import solve_nlp

    problem, _ = cra_problem(arch_assembly(), mu=0.7, density=1.0, d_bnd=1e-3, eps=1e-4)
    result = solve_nlp(problem, backend="native", options=_CRA_OPTIONS)
    assert result.success, result.status_message
    assert result.iterations < 1500, "lost iteration headroom: {} iterations".format(result.iterations)
