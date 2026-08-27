"""CRA view style using compas_viewer"""

from math import sqrt

import numpy as np
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
from compas.geometry import Rotation
from compas.geometry import Translation
from compas.geometry import Vector
from compas.geometry import is_coplanar
from compas_viewer import Viewer
from compas_viewer.config import Config
from compas_viewer.scene import Collection


class Arrow:
    """A force arrow as a mesh: cylinder shaft, conical head.

    This is the compas_view2 ``Arrow`` shape the original visualization used, with its
    parameters: the widths are world units, the head takes ``head_portion`` of the
    length. The compas_viewer ``VectorObject`` is no substitute - it draws a fixed
    pixel-width line with a four-sided pyramid head that degenerates to a sliver for
    axis-parallel directions.
    """

    def __init__(self, position=[0, 0, 0], direction=[0, 0, 1], head_portion=0.2, head_width=0.07, body_width=0.02):
        super().__init__()
        self.position = Vector(*position)
        self.direction = Vector(*direction)
        self.head_portion = head_portion
        self.head_width = head_width
        self.body_width = body_width

    def mesh(self):
        length = self.direction.length
        if length == 0:
            return None
        zaxis = self.direction.unitized()
        body_length = length * (1 - self.head_portion)
        # widths scale with the length, exactly as the compas_view2 Arrow scaled its
        # unit shape by the direction vector: a short force renders as a thin pin, a
        # long one as a bold arrow, and relative proportions stay comparable
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

    def add_to_scene(self, viewer, facecolor: Color, opacity=0.999):
        mesh = self.mesh()
        if mesh is None:
            return
        # opacity just under 1 on purpose: it moves the arrow into the renderer's
        # transparency pass, which draws after the semi-transparent block faces - so an
        # arrow inside a block stays bold, as it was in compas_view2, instead of being
        # washed out by the faces compositing over it
        viewer.scene.add(mesh, facecolor=facecolor, show_lines=False, opacity=min(opacity, 0.999))


def draw_blocks(assembly, viewer, edge=True, tol=0.0):
    supports = []
    blocks = []
    supportedges = []
    blockedges = []
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, "block")
        if assembly.graph.node_attribute(node, "is_support"):
            supports.append(block)
        else:
            blocks.append(block)
        if not edge:
            continue
        for edge in block.edges():
            if tol != 0.0:
                fkeys = block.edge_faces(edge)
                ps = [
                    block.face_center(fkeys[0]),
                    block.face_center(fkeys[1]),
                    *block.edge_coordinates(edge),
                ]

                if is_coplanar(ps, tol=tol):
                    continue
            if assembly.graph.node_attribute(node, "is_support"):
                supportedges.append(Line(*block.edge_coordinates(edge)))
            else:
                blockedges.append(Line(*block.edge_coordinates(edge)))
    if len(blocks) != 0:
        viewer.scene.add(
            Collection(blocks),
            show_faces=True,
            show_lines=False,
            opacity=0.6,
            facecolor=Color(0.9, 0.9, 0.9),
        )
    if len(supports) != 0:
        viewer.scene.add(
            Collection(supports),
            show_faces=True,
            show_lines=False,
            opacity=0.5,
            facecolor=Color.from_hex("#f79d84"),
        )
    # edges share one thickness; each set takes the color of the meshes it outlines:
    # dark gray on the light gray blocks, salmon on the salmon supports
    if len(blockedges) != 0:
        viewer.scene.add(Collection(blockedges), linecolor=Color(0.3, 0.3, 0.3), linewidth=1.5)
    if len(supportedges) != 0:
        viewer.scene.add(Collection(supportedges), linecolor=Color.from_hex("#f79d84"), linewidth=1.5)


def draw_interfaces(assembly, viewer):
    interfaces = []
    faces = []
    for edge in assembly.graph.edges():
        interface = assembly.graph.edge_attribute(edge, "interface")
        if interface is not None:
            corners = np.array(interface.points)
            faces.append(Polyline(np.vstack((corners, corners[0]))))
            if assembly.graph.node_attribute(edge[0], "is_support") or assembly.graph.node_attribute(
                edge[1], "is_support"
            ):
                continue
            polygon = Polygon(interface.points)
            interfaces.append(Mesh.from_polygons([polygon]))
        if assembly.graph.edge_attribute(edge, "interfaces") is None:
            continue
        for subinterface in assembly.graph.edge_attribute(edge, "interfaces"):
            corners = np.array(subinterface.points)
            faces.append(Polyline(np.vstack((corners, corners[0]))))
            polygon = Polygon(subinterface.points)
            interfaces.append(Mesh.from_polygons([polygon]))

    if len(interfaces) != 0:
        viewer.scene.add(
            Collection(interfaces),
            show_lines=False,
            show_points=False,
            facecolor=(0.8, 0.8, 0.8),
        )
    if len(faces) != 0:
        viewer.scene.add(
            Collection(faces),
            linecolor=Color.from_hex("#fac05e"),
            linewidth=10,
            pointsize=10,
            show_points=True,
            pointcolor=(0, 0, 0),
        )


def draw_forces(assembly, viewer, scale=1.0, resultant=True, nodal=False):
    locs = []
    res_np = []
    res_nn = []
    fnp = []
    fnn = []
    ft = []
    for edge in assembly.graph.edges():
        interface = assembly.graph.edge_attribute(edge, "interface")
        if interface is None:
            break
        forces = interface.forces
        if forces is None:
            continue
        corners = np.array(interface.points)
        frame = interface.frame
        w, u, v = frame.zaxis, frame.xaxis, frame.yaxis
        if nodal:
            for i, corner in enumerate(corners):
                pt = Point(*corner)
                force = forces[i]["c_np"] - forces[i]["c_nn"]
                p1 = pt + w * force * 0.5 * scale
                p2 = pt - w * force * 0.5 * scale
                if force >= 0:
                    fnn.append(Line(p1, p2))
                else:
                    fnp.append(Line(p1, p2))
                ft_uv = (u * forces[i]["c_u"] + v * forces[i]["c_v"]) * 0.5 * scale
                p1 = pt + ft_uv
                p2 = pt - ft_uv
                ft.append(Line(p1, p2))
        if resultant:
            sum_n = sum(force["c_np"] - force["c_nn"] for force in forces)
            sum_u = sum(force["c_u"] for force in forces)
            sum_v = sum(force["c_v"] for force in forces)
            if sum_n == 0:
                continue
            resultant_pos = np.average(
                np.array(corners),
                axis=0,
                weights=[force["c_np"] - force["c_nn"] for force in forces],
            )
            locs.append(Point(*resultant_pos))
            # resultant
            resultant_f = (w * sum_n + u * sum_u + v * sum_v) * 0.5 * scale
            p1 = resultant_pos + resultant_f
            p2 = resultant_pos - resultant_f
            if sum_n >= 0:
                res_np.append(Line(p1, p2))
            else:
                res_nn.append(Line(p1, p2))
    if len(locs) != 0:
        viewer.scene.add(Collection(locs), pointsize=12, pointcolor=Color.from_hex("#386641"))
    if len(res_np) != 0:
        viewer.scene.add(Collection(res_np), linewidth=8, linecolor=Color(0, 0.3, 0))
    if len(res_nn) != 0:
        viewer.scene.add(Collection(res_nn), linewidth=8, linecolor=Color(0.8, 0, 0))
    if len(fnn) != 0:
        viewer.scene.add(Collection(fnn), linewidth=5, linecolor=Color.from_hex("#00468b"))
    if len(fnp) != 0:
        viewer.scene.add(Collection(fnp), linewidth=5, linecolor=Color(1, 0, 0))
    if len(ft) != 0:
        viewer.scene.add(Collection(ft), linewidth=5, linecolor=Color(1.0, 0.5, 0.0))


def draw_forcesline(assembly, viewer, scale=1.0, resultant=True, nodal=False):
    locs = []
    res_np = []
    res_nn = []
    fnp = []
    fnn = []
    ft = []
    # total_reaction = 0
    for edge in assembly.graph.edges():
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            forces = interface.forces
            if forces is None:
                continue
            corners = np.array(interface.points)
            frame = interface.frame
            w, u, v = frame.zaxis, frame.xaxis, frame.yaxis
            if nodal:
                for i, corner in enumerate(corners):
                    pt = Point(*corner)
                    force = forces[i]["c_np"] - forces[i]["c_nn"]
                    p1 = pt + w * force * 0.5 * scale
                    p2 = pt - w * force * 0.5 * scale
                    if force >= 0:
                        fnn.append(Line(p1, p2))
                    else:
                        fnp.append(Line(p1, p2))
                    ft_uv = (u * forces[i]["c_u"] + v * forces[i]["c_v"]) * 0.5 * scale
                    p1 = pt + ft_uv
                    p2 = pt - ft_uv
                    ft.append(Line(p1, p2))
            if resultant:
                sum_n = sum(force["c_np"] - force["c_nn"] for force in forces)
                # tension when the NET normal force is tensile - see the comment in
                # draw_forcesdirect: per-vertex tests false-positive on the solver's
                # relaxed-complementarity noise (vertex forces negative by ~1e-3 on
                # compressed interfaces), while the net sign colors the resultant by
                # what it actually depicts
                is_tension = sum_n < 0
                sum_u = sum(force["c_u"] for force in forces)
                sum_v = sum(force["c_v"] for force in forces)
                if sum_n == 0:
                    continue
                resultant_pos = np.average(
                    np.array(corners),
                    axis=0,
                    weights=[force["c_np"] - force["c_nn"] for force in forces],
                )
                locs.append(Point(*resultant_pos))
                # resultant
                resultant_f = (w * sum_n + u * sum_u + v * sum_v) * 0.5 * scale
                # print((w * sum_n + u * sum_u + v * sum_v).length * 100000, "edge: ", edge)

                # if assembly.graph.node_attribute(edge[0], "is_support") or assembly.graph.node_attribute(
                #     edge[1], "is_support"
                # ):
                #     print((w * sum_n + u * sum_u + v * sum_v).z)
                # total_reaction += abs((w * sum_n + u * sum_u + v * sum_v).z * 100000)

                p1 = resultant_pos + resultant_f
                p2 = resultant_pos - resultant_f

                if not is_tension:
                    res_np.append(Line(p1, p2))
                else:
                    res_nn.append(Line(p1, p2))
    if len(locs) != 0:
        viewer.scene.add(Collection(locs), pointsize=12, pointcolor=Color.from_hex("#386641"))
    if len(res_np) != 0:
        viewer.scene.add(Collection(res_np), linewidth=8, linecolor=Color(0, 0.3, 0))
    if len(res_nn) != 0:
        viewer.scene.add(Collection(res_nn), linewidth=8, linecolor=Color(0.8, 0, 0))
    if len(fnn) != 0:
        viewer.scene.add(Collection(fnn), linewidth=5, linecolor=Color.from_hex("#00468b"))
    if len(fnp) != 0:
        viewer.scene.add(Collection(fnp), linewidth=5, linecolor=Color(1, 0, 0))
    if len(ft) != 0:
        viewer.scene.add(Collection(ft), linewidth=5, linecolor=Color(1.0, 0.5, 0.0))
    # print("total reaction: ", total_reaction)


def draw_forcesdirect(assembly, viewer, scale=1.0, resultant=True, nodal=False):
    locs = []
    res_np = []
    res_nn = []
    fnp = []
    fnn = []
    ft = []
    for edge in assembly.graph.edges():
        thres = 1e-6
        if assembly.graph.node_attribute(edge[0], "is_support") and not assembly.graph.node_attribute(
            edge[1], "is_support"
        ):
            flip = False
        else:
            flip = True
        if assembly.graph.edge_attribute(edge, "interfaces") is None:
            continue
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            forces = interface.forces
            if forces is None:
                continue
            corners = np.array(interface.points)
            frame = interface.frame
            w, u, v = frame.zaxis, frame.xaxis, frame.yaxis
            if nodal:
                for i, corner in enumerate(corners):
                    pt = Point(*corner)
                    force = forces[i]["c_np"] - forces[i]["c_nn"]
                    if (w * force * scale).length == 0:
                        continue
                    if flip:
                        f = Arrow(pt, w * force * scale * -1)
                    else:
                        f = Arrow(pt, w * force * scale)
                    if force >= 0:
                        fnp.append(f)
                    else:
                        fnn.append(f)
                    ft_uv = (u * forces[i]["c_u"] + v * forces[i]["c_v"]) * scale
                    if ft_uv.length == 0:
                        continue
                    if flip:
                        f = Arrow(pt, ft_uv * -1)
                    else:
                        f = Arrow(pt, ft_uv)
                    ft.append(f)
            if resultant:
                sum_n = sum(force["c_np"] - force["c_nn"] for force in forces)
                # The resultant arrow depicts the NET interface force, so its color
                # follows the sign of the net normal component: tensile (negative)
                # renders red, compressive renders green - the same rule draw_forces
                # uses. The old test flagged the interface red as soon as ANY vertex
                # was more than 1e-5 in tension, but the solver relaxes
                # complementarity (f_n (d_n + eps) = 0 with eps ~ 1e-4) and stops at
                # IPOPT's "acceptable" tolerance, so vertex normal forces routinely
                # come out negative by up to ~1e-3 on interfaces that are in overall
                # compression - a false positive (the shelf example rendered red
                # although its published screenshot is green). Genuine tension - a
                # hanging or peeling interface - makes the net normal force itself
                # negative and still renders red.
                is_tension = sum_n < 0
                sum_u = sum(force["c_u"] for force in forces)
                sum_v = sum(force["c_v"] for force in forces)
                if abs(sum_n) <= thres:
                    resultant_pos = np.average(
                        np.array(corners),
                        axis=0,
                        weights=[sqrt(force["c_u"] ** 2 + force["c_v"] ** 2) for force in forces],
                    )
                    friction = True
                else:
                    resultant_pos = np.average(
                        np.array(corners),
                        axis=0,
                        weights=[force["c_np"] - force["c_nn"] for force in forces],
                    )
                    friction = False
                resultant_f = (w * sum_n + u * sum_u + v * sum_v) * scale
                if resultant_f.length >= thres:
                    locs.append(Point(*resultant_pos))
                if flip:
                    f = Arrow(resultant_pos, resultant_f * -1)
                else:
                    f = Arrow(resultant_pos, resultant_f)
                if friction:
                    f.add_to_scene(viewer, facecolor=(1.0, 0.5, 0.0))
                if not is_tension:
                    res_np.append(f)
                else:
                    res_nn.append(f)
    if len(locs) != 0:
        viewer.scene.add(Collection(locs), pointsize=12, pointcolor=Color.from_hex("#386641"))
    if len(res_np) != 0:
        for arrow in res_np:
            arrow.add_to_scene(viewer, facecolor=Color.from_hex("#386641"))
    if len(res_nn) != 0:
        for arrow in res_nn:
            arrow.add_to_scene(viewer, facecolor=Color(0.8, 0, 0))
    if len(fnp) != 0:
        for arrow in fnp:
            arrow.add_to_scene(viewer, facecolor=Color.from_hex("#00468b"), opacity=0.5)
    if len(fnn) != 0:
        for arrow in fnn:
            arrow.add_to_scene(viewer, facecolor=Color(1, 0, 0), opacity=0.5)
    if len(ft) != 0:
        for arrow in ft:
            arrow.add_to_scene(viewer, facecolor=Color(1.0, 0.5, 0.0), opacity=0.5)


def draw_displacements(assembly, viewer, dispscale=1.0, tol=0.0):
    blocks = []
    nodes = []
    for node in assembly.graph.nodes():
        if assembly.graph.node_attribute(node, "is_support"):
            continue
        block = assembly.graph.node_attribute(node, "block")
        displacement = assembly.graph.node_attribute(node, "displacement")
        if displacement is None:
            continue
        displacement = np.array(displacement) * dispscale
        vec = (
            np.array([1, 0, 0]) * displacement[3]
            + np.array([0, 1, 0]) * displacement[4]
            + np.array([0, 0, 1]) * displacement[5]
        ).tolist()
        R = Rotation.from_axis_angle_vector(vec, point=block.center())
        T = Translation.from_vector(displacement[0:3])
        new_block = block.transformed(R).transformed(T)
        nodes.append(Point(*new_block.center()))
        for edge in block.edges():
            if tol != 0.0:
                fkeys = block.edge_faces(edge)
                ps = [
                    block.face_center(fkeys[0]),
                    block.face_center(fkeys[1]),
                    *block.edge_coordinates(edge),
                ]
                if is_coplanar(ps, tol=tol):
                    continue
            blocks.append(Line(*new_block.edge_coordinates(edge)))
    if len(blocks) != 0:
        viewer.scene.add(Collection(blocks), linewidth=1, linecolor=Color(0.7, 0.7, 0.7))
    if len(nodes) != 0:
        viewer.scene.add(Collection(nodes), pointcolor=Color(0.7, 0.7, 0.7))


def draw_weights(assembly, viewer, scale=1.0, density=1.0):
    weights = []
    blocks = []
    supports = []
    # total_weights = 0
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, "block")
        if assembly.graph.node_attribute(node, "is_support"):
            supports.append(Point(*block.center()))
            continue
        d = block.attributes["density"] if "density" in block.attributes else density
        weights.append(
            Arrow(
                block.center(),
                [0, 0, -block.volume() * d * scale],
            )
        )
        # print("self-weight", -block.volume() * density)
        # total_weights += block.volume() * 2500 * 9.8
        blocks.append(Point(*block.center()))

    # print("total self-weight: ", total_weights)

    if len(supports) != 0:
        viewer.scene.add(Collection(supports), pointsize=20, pointcolor=Color.from_hex("#ee6352"))
    if len(blocks) != 0:
        viewer.scene.add(Collection(blocks), pointsize=30, pointcolor=Color.from_hex("#3284a0"))
    if len(weights) != 0:
        for weight in weights:
            weight.add_to_scene(viewer, facecolor=Color.from_hex("#59cd90"))


def cra_view(
    assembly,
    scale=1.0,
    density=1.0,
    dispscale=1.0,
    tol=1e-5,
    grid=False,
    resultant=True,
    nodal=False,
    edge=True,
    blocks=True,
    interfaces=True,
    forces=True,
    forcesdirect=True,
    forcesline=False,
    weights=True,
    displacements=True,
):
    """CRA Viewer, creating new viewer.

    Parameters
    ----------
    assembly : :class:`~compas_assembly.datastructures.Assembly`
        The rigid block assembly.
    scale : float, optional
        Force scale.
    density : float, optional
        Density of the block material.
    dispscale : float, optional
        virtual displacement scale.
    tol : float, optional
        Tolerance value to consider faces to be planar.
    grid : bool, optional
        Show view grid.
    resultant : bool, optional
        Plot resultant forces.
    nodal : bool, optional
        Plot nodal forces.
    edge : bool, optional
        Plot block edges.
    blocks : bool, optional
        Plot block.
    interfaces : bool, optional
        Plot interfaces.
    forces : bool, optional
        Plot forces.
    forcesdirect : bool, optional
        Plot forces as vectors.
    forcesline : bool, optional
        Plot forces as lines.
    weights : bool, optional
        Plot block self weight as vectors.
    displacements : bool, optional
        Plot virtual displacements.

    Returns
    -------
    None
    """

    viewer = Viewer(config=Config(vectorsize=0.15))

    if blocks:
        draw_blocks(assembly, viewer, edge, tol)
    if interfaces:
        draw_interfaces(assembly, viewer)
    if forces:
        draw_forces(assembly, viewer, scale, resultant, nodal)
    if forcesdirect:
        draw_forcesdirect(assembly, viewer, scale, resultant, nodal)
    if forcesline:
        draw_forcesline(assembly, viewer, scale, resultant, nodal)
    if weights:
        draw_weights(assembly, viewer, scale, density)
    if displacements:
        draw_displacements(assembly, viewer, dispscale, tol)

    viewer.show()


def cra_view_ex(
    viewer,
    assembly,
    scale=1.0,
    density=1.0,
    dispscale=1.0,
    tol=1e-5,
    resultant=True,
    nodal=False,
    edge=True,
    blocks=True,
    interfaces=True,
    forces=True,
    forcesdirect=True,
    forcesline=False,
    weights=True,
    displacements=True,
):
    """CRA Viewer using existing view.

    Parameters
    ----------
    viewer : compas_viewer.Viewer
        External viewer object.
    assembly : :class:`~compas_assembly.datastructures.Assembly`
        The rigid block assembly.
    scale : float, optional
        Force scale.
    density : float, optional
        Density of the block material.
    dispscale : float, optional
        virtual displacement scale.
    tol : float, optional
        Tolerance value to consider faces to be planar.
    resultant : bool, optional
        Plot resultant forces.
    nodal : bool, optional
        Plot nodal forces.
    edge : bool, optional
        Plot block edges.
    blocks : bool, optional
        Plot block.
    interfaces : bool, optional
        Plot interfaces.
    forces : bool, optional
        Plot forces.
    forcesdirect : bool, optional
        Plot forces as vectors.
    forcesline : bool, optional
        Plot forces as lines.
    weights : bool, optional
        Plot block self weight as vectors.
    displacements : bool, optional
        Plot virtual displacements.

    Returns
    -------
    None
    """

    if blocks:
        draw_blocks(assembly, viewer, edge, tol)
    if interfaces:
        draw_interfaces(assembly, viewer)
    if forces:
        draw_forces(assembly, viewer, scale, resultant, nodal)
    if forcesdirect:
        draw_forcesdirect(assembly, viewer, scale, resultant, nodal)
    if forcesline:
        draw_forcesline(assembly, viewer, scale, resultant, nodal)
    if weights:
        draw_weights(assembly, viewer, scale, density)
    if displacements:
        draw_displacements(assembly, viewer, dispscale, tol)
