"""The viewer module against the installed compas_viewer, without a window.

The scene objects are fully constructed - only the blocking show() is stubbed - so an
API drift in compas_viewer (scene.add rejecting lists, renamed kwargs) fails here
instead of on a user's screen.
"""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402

pytest.importorskip("compas_viewer", reason="the viz extra is not installed")

from compas.datastructures import Mesh  # noqa: E402
from compas.geometry import Box  # noqa: E402
from compas.geometry import Frame  # noqa: E402
from compas.geometry import Translation  # noqa: E402
from compas_assembly.datastructures import Block  # noqa: E402
from compas_viewer import Viewer  # noqa: E402

from compas_cra.datastructures import CRA_Assembly  # noqa: E402
from compas_cra.equilibrium import cra_solve  # noqa: E402
from compas_cra.viewers import cra_view  # noqa: E402


def solved_cubes():
    support = Box(1, 1, 1)
    free = Box(1, 1, 1, frame=Frame.worldXY().transformed(Translation.from_vector([0, 0, 1])))
    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), node=0)
    assembly.add_block(Block.from_shape(free), node=1)
    assembly.set_boundary_conditions([0])
    interface = Mesh()
    for i, c in enumerate([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5]]):
        interface.add_vertex(key=i, x=c[0], y=c[1], z=c[2])
    interface.add_face([0, 1, 2, 3])
    assembly.add_interfaces_from_meshes([interface], 0, 1)
    cra_solve(assembly, density=1)
    return assembly


def test_cra_view_builds_the_scene(monkeypatch):
    seen = {}
    monkeypatch.setattr(Viewer, "show", lambda self: seen.setdefault("objects", len(self.scene.objects)))
    cra_view(solved_cubes(), forcesline=True, nodal=True)
    assert seen["objects"] > 0
