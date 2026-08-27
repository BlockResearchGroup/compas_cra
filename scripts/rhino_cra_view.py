#! python3
# venv: compas-cra
# r: compas_cra>=0.7.2
"""Rhino port of the compas_cra desktop viewer (compas_cra.viewers.cra_view).

Draws a solved CRA assembly into the active Rhino document through the COMPAS scene API,
mirroring the desktop viewer: block and support meshes with their edge lines, interface
outlines, resultant/nodal force arrows (cylinder shaft + conical head, widths scaled by
arrow length), and self-weight arrows - all in the same colors. Objects are sorted onto
layers ("<prefix>::Blocks", "<prefix>::Forces::Resultant", ...) so each group can be
toggled in the Rhino layer panel.

This module imports cleanly outside Rhino (only ``compas`` is needed), so it can be
syntax- and logic-checked headless; the actual drawing requires Rhino, where the COMPAS
scene resolves every item to its Rhino scene object.
"""

from math import sqrt

from compas.colors import Color
from compas.datastructures import Mesh
from compas.geometry import Cone
from compas.geometry import Cylinder
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Plane
from compas.geometry import Point
from compas.geometry import Polygon
from compas.geometry import Polyline
from compas.geometry import Vector
from compas.geometry import is_coplanar
from compas.scene import Scene

# desktop viewer colors (compas_cra.viewers.cra_view)
COLOR_BLOCK = Color(0.9, 0.9, 0.9)
COLOR_BLOCK_EDGE = Color(0.0, 0.0, 0.0)
COLOR_SUPPORT = Color.from_hex("#f79d84")
COLOR_INTERFACE_FACE = Color(0.8, 0.8, 0.8)
COLOR_INTERFACE_OUTLINE = Color.from_hex("#fac05e")
COLOR_INTERFACE_POINT = Color(0.0, 0.0, 0.0)
COLOR_COMPRESSION = Color.from_hex("#386641")
COLOR_TENSION = Color(0.8, 0.0, 0.0)
COLOR_NODAL_COMPRESSION = Color.from_hex("#00468b")
COLOR_NODAL_TENSION = Color(1.0, 0.0, 0.0)
COLOR_FRICTION = Color(1.0, 0.5, 0.0)
COLOR_WEIGHT = Color.from_hex("#59cd90")
COLOR_SUPPORT_NODE = Color.from_hex("#ee6352")
COLOR_BLOCK_NODE = Color.from_hex("#3284a0")


class Arrow:
    """A force arrow as a mesh: cylinder shaft plus conical head.

    Same parameters as the desktop viewer's ``Arrow``: the head takes ``head_portion``
    of the total length, and both widths are scaled by the arrow length, so a short
    force renders as a thin pin and a long one as a bold arrow.
    """

    def __init__(self, position=(0, 0, 0), direction=(0, 0, 1), head_portion=0.2, head_width=0.07, body_width=0.02):
        self.position = Vector(*position)
        self.direction = Vector(*direction)
        self.head_portion = head_portion
        self.head_width = head_width
        self.body_width = body_width

    def mesh(self):
        """Build the arrow mesh, or return None for a zero-length arrow."""
        length = self.direction.length
        if length == 0:
            return None
        zaxis = self.direction.unitized()
        body_length = length * (1 - self.head_portion)
        shaft = Cylinder(
            radius=self.body_width * length,
            height=body_length,
            frame=Frame.from_plane(Plane(self.position + zaxis * (body_length / 2), zaxis)),
        )
        head = Cone(
            radius=self.head_width * length,
            height=length * self.head_portion,
            frame=Frame.from_plane(Plane(self.position + zaxis * body_length, zaxis)),
        )
        mesh = Mesh.from_shape(shaft, u=32)
        mesh.join(Mesh.from_shape(head, u=32))
        return mesh


def _add_arrow(scene, arrow, color, layer):
    mesh = arrow.mesh()
    if mesh is None:
        return
    scene.add(mesh, color=color, layer=layer)


def _weighted_average(points, weights):
    total = sum(weights)
    if total == 0:  # degenerate: fall back to the plain centroid
        weights = [1.0] * len(points)
        total = float(len(points))
    return Point(
        sum(p[0] * w for p, w in zip(points, weights)) / total,
        sum(p[1] * w for p, w in zip(points, weights)) / total,
        sum(p[2] * w for p, w in zip(points, weights)) / total,
    )


def draw_blocks(assembly, scene, edge=True, tol=0.0, layer_prefix="CRA"):
    """Add block and support meshes (plus their edge lines) to the scene."""
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, "block")
        is_support = assembly.graph.node_attribute(node, "is_support")
        if is_support:
            scene.add(block, color=COLOR_SUPPORT, layer=layer_prefix + "::Supports")
        else:
            scene.add(block, color=COLOR_BLOCK, layer=layer_prefix + "::Blocks")
        if not edge:
            continue
        for block_edge in block.edges():
            if tol != 0.0:
                fkeys = block.edge_faces(block_edge)
                ps = [
                    block.face_center(fkeys[0]),
                    block.face_center(fkeys[1]),
                    *block.edge_coordinates(block_edge),
                ]
                if is_coplanar(ps, tol=tol):
                    continue
            line = Line(*block.edge_coordinates(block_edge))
            if is_support:
                scene.add(line, color=COLOR_SUPPORT, layer=layer_prefix + "::Supports")
            else:
                scene.add(line, color=COLOR_BLOCK_EDGE, layer=layer_prefix + "::Blocks")


def _draw_interface(scene, interface, layer, faces=True):
    corners = [Point(*point) for point in interface.points]
    if faces:
        scene.add(Mesh.from_polygons([Polygon(interface.points)]), color=COLOR_INTERFACE_FACE, layer=layer)
    scene.add(Polyline(corners + corners[:1]), color=COLOR_INTERFACE_OUTLINE, layer=layer)
    for corner in corners:
        scene.add(corner, color=COLOR_INTERFACE_POINT, layer=layer)


def draw_interfaces(assembly, scene, layer_prefix="CRA"):
    """Add interface faces, outlines, and corner points to the scene."""
    layer = layer_prefix + "::Interfaces"
    for edge in assembly.graph.edges():
        interface = assembly.graph.edge_attribute(edge, "interface")
        if interface is not None:
            support = assembly.graph.node_attribute(edge[0], "is_support") or assembly.graph.node_attribute(
                edge[1], "is_support"
            )
            _draw_interface(scene, interface, layer, faces=not support)
        for subinterface in assembly.graph.edge_attribute(edge, "interfaces") or []:
            _draw_interface(scene, subinterface, layer)


def draw_forcesdirect(assembly, scene, scale=1.0, resultant=True, nodal=False, layer_prefix="CRA"):
    """Add resultant and/or nodal force arrows to the scene, as in the desktop viewer."""
    layer_resultant = layer_prefix + "::Forces::Resultant"
    layer_nodal = layer_prefix + "::Forces::Nodal"
    layer_friction = layer_prefix + "::Forces::Friction"
    thres = 1e-6
    locs = []
    res_np = []
    res_nn = []
    fnp = []
    fnn = []
    ft = []
    for edge in assembly.graph.edges():
        support_0 = assembly.graph.node_attribute(edge[0], "is_support")
        support_1 = assembly.graph.node_attribute(edge[1], "is_support")
        flip = not (support_0 and not support_1)
        interfaces = assembly.graph.edge_attribute(edge, "interfaces")
        if interfaces is None:
            continue
        for interface in interfaces:
            forces = interface.forces
            if forces is None:
                continue
            corners = [Point(*point) for point in interface.points]
            frame = interface.frame
            w, u, v = frame.zaxis, frame.xaxis, frame.yaxis
            if nodal:
                for i, pt in enumerate(corners):
                    force = forces[i]["c_np"] - forces[i]["c_nn"]
                    fn = w * force * scale
                    if fn.length == 0:
                        continue
                    arrow = Arrow(pt, fn * -1 if flip else fn)
                    if force >= 0:
                        fnp.append(arrow)
                    else:
                        fnn.append(arrow)
                    ft_uv = (u * forces[i]["c_u"] + v * forces[i]["c_v"]) * scale
                    if ft_uv.length == 0:
                        continue
                    ft.append(Arrow(pt, ft_uv * -1 if flip else ft_uv))
            if resultant:
                sum_n = sum(force["c_np"] - force["c_nn"] for force in forces)
                # net tension decides the color - the rule the published screenshots used;
                # per-vertex epsilon noise from the relaxed complementarity must not flip it
                is_tension = sum_n < 0
                sum_u = sum(force["c_u"] for force in forces)
                sum_v = sum(force["c_v"] for force in forces)
                if abs(sum_n) <= thres:  # pure friction interface
                    weights = [sqrt(force["c_u"] ** 2 + force["c_v"] ** 2) for force in forces]
                    friction = True
                else:
                    weights = [force["c_np"] - force["c_nn"] for force in forces]
                    friction = False
                pos = _weighted_average(corners, weights)
                resultant_f = (w * sum_n + u * sum_u + v * sum_v) * scale
                if resultant_f.length >= thres:
                    locs.append(pos)
                arrow = Arrow(pos, resultant_f * -1 if flip else resultant_f)
                if friction:
                    _add_arrow(scene, arrow, COLOR_FRICTION, layer_friction)
                if not is_tension:
                    res_np.append(arrow)
                else:
                    res_nn.append(arrow)
    for loc in locs:
        scene.add(loc, color=COLOR_COMPRESSION, layer=layer_resultant)
    for arrow in res_np:
        _add_arrow(scene, arrow, COLOR_COMPRESSION, layer_resultant)
    for arrow in res_nn:
        _add_arrow(scene, arrow, COLOR_TENSION, layer_resultant)
    for arrow in fnp:
        _add_arrow(scene, arrow, COLOR_NODAL_COMPRESSION, layer_nodal)
    for arrow in fnn:
        _add_arrow(scene, arrow, COLOR_NODAL_TENSION, layer_nodal)
    for arrow in ft:
        _add_arrow(scene, arrow, COLOR_FRICTION, layer_friction)


def draw_weights(assembly, scene, scale=1.0, density=1.0, layer_prefix="CRA"):
    """Add self-weight arrows and block/support center points to the scene."""
    layer = layer_prefix + "::Weights"
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, "block")
        if assembly.graph.node_attribute(node, "is_support"):
            scene.add(Point(*block.center()), color=COLOR_SUPPORT_NODE, layer=layer)
            continue
        d = block.attributes["density"] if "density" in block.attributes else density
        _add_arrow(scene, Arrow(block.center(), [0, 0, -block.volume() * d * scale]), COLOR_WEIGHT, layer)
        scene.add(Point(*block.center()), color=COLOR_BLOCK_NODE, layer=layer)


def _rhino_redraw(enable):
    try:
        import rhinoscriptsyntax as rs
    except ImportError:  # not running inside Rhino
        return
    rs.EnableRedraw(enable)
    if enable:
        rs.ZoomExtents(all=True)


def cra_view_rhino(
    assembly,
    scale=1.0,
    density=1.0,
    tol=1e-5,
    resultant=True,
    nodal=False,
    edge=True,
    blocks=True,
    interfaces=True,
    forcesdirect=True,
    weights=True,
    layer_prefix="CRA",
):
    """Draw a solved CRA assembly into the Rhino document, like the desktop ``cra_view``.

    Parameters
    ----------
    assembly : :class:`~compas_cra.datastructures.CRA_Assembly`
        The solved rigid block assembly.
    scale : float, optional
        Force scale.
    density : float, optional
        Density of the block material (for the self-weight arrows).
    tol : float, optional
        Tolerance below which coplanar block edges are not drawn.
    resultant : bool, optional
        Plot resultant force arrows.
    nodal : bool, optional
        Plot nodal force arrows.
    edge : bool, optional
        Plot block edges.
    blocks : bool, optional
        Plot blocks and supports.
    interfaces : bool, optional
        Plot interfaces.
    forcesdirect : bool, optional
        Plot forces as arrows.
    weights : bool, optional
        Plot block self-weight as arrows.
    layer_prefix : str, optional
        Rhino layers are created as "<layer_prefix>::<group>".

    Returns
    -------
    :class:`compas.scene.Scene`
        The scene that was drawn.
    """
    scene = Scene()
    if blocks:
        draw_blocks(assembly, scene, edge=edge, tol=tol, layer_prefix=layer_prefix)
    if interfaces:
        draw_interfaces(assembly, scene, layer_prefix=layer_prefix)
    if forcesdirect:
        draw_forcesdirect(assembly, scene, scale=scale, resultant=resultant, nodal=nodal, layer_prefix=layer_prefix)
    if weights:
        draw_weights(assembly, scene, scale=scale, density=density, layer_prefix=layer_prefix)
    _rhino_redraw(False)
    try:
        scene.draw()
    finally:
        _rhino_redraw(True)
    return scene
