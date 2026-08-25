#! python3
# venv: compas-cra
# r: compas_cra>=0.7.2
"""Parametric masonry arch: CRA equilibrium, baked into the Rhino document."""

import rhinoscriptsyntax as rs

from compas_cra.algorithms import assembly_interfaces_numpy
from compas_cra.equilibrium import cra_solve
from compas_cra.geometry import Arch

HEIGHT = 5.0  # rise, at most SPAN / 2
SPAN = 10.0
THICKNESS = 0.5
DEPTH = 0.5
NUM_BLOCKS = 20
MU = 0.7  # friction coefficient
FORCE_SCALE = 0.5
FORCE_RADIUS = 0.05  # pipe radius of the drawn force lines

assembly = Arch(
    height=HEIGHT,
    span=SPAN,
    thickness=THICKNESS,
    depth=DEPTH,
    num_blocks=NUM_BLOCKS,
).assembly()

assembly_interfaces_numpy(assembly, nmax=10, amin=1e-2, tmax=1e-2)
cra_solve(assembly, mu=MU)


def ensure_layer(name, color):
    if not rs.IsLayer(name):
        rs.AddLayer(name, color)
    return name


def draw_block(block, layer):
    guids = [rs.AddLine(*block.edge_coordinates(edge)) for edge in block.edges()]
    rs.ObjectLayer(guids, layer)
    rs.AddObjectsToGroup(guids, rs.AddGroup())


blocks_layer = ensure_layer("CRA-Arch::Blocks", (0, 0, 0))
supports_layer = ensure_layer("CRA-Arch::Supports", (247, 157, 132))
interfaces_layer = ensure_layer("CRA-Arch::Interfaces", (0, 70, 139))
compression_layer = ensure_layer("CRA-Arch::Compression", (0, 120, 0))
tension_layer = ensure_layer("CRA-Arch::Tension", (200, 0, 0))

rs.EnableRedraw(False)

for node in assembly.graph.nodes():
    block = assembly.graph.node_attribute(node, "block")
    is_support = assembly.graph.node_attribute(node, "is_support")
    draw_block(block, supports_layer if is_support else blocks_layer)

for edge in assembly.graph.edges():
    for interface in assembly.graph.edge_attribute(edge, "interfaces") or []:
        corners = list(interface.points)
        rs.ObjectLayer(rs.AddPolyline(corners + [corners[0]]), interfaces_layer)

        forces = interface.forces
        if forces is None:
            continue

        w, u, v = interface.frame.zaxis, interface.frame.xaxis, interface.frame.yaxis
        normals = [f["c_np"] - f["c_nn"] for f in forces]
        sum_n = sum(normals)
        if sum_n == 0:
            continue
        sum_u = sum(f["c_u"] for f in forces)
        sum_v = sum(f["c_v"] for f in forces)
        pos = [sum(c[i] * n for c, n in zip(corners, normals)) / sum_n for i in range(3)]
        f = (w * sum_n + u * sum_u + v * sum_v) * 0.5 * FORCE_SCALE
        layer = compression_layer if sum_n >= 0 else tension_layer
        line = rs.AddLine([pos[i] + f[i] for i in range(3)], [pos[i] - f[i] for i in range(3)])
        pipe = rs.AddPipe(line, 0, FORCE_RADIUS, cap=1)
        rs.DeleteObject(line)
        rs.ObjectLayer(pipe, layer)
        rs.ObjectLayer(rs.AddTextDot("{:.2f}".format(abs(sum_n)), pos), layer)

rs.EnableRedraw(True)
rs.ZoomExtents(all=True)
