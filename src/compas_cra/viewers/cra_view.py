"""CRA view style using compas_view2"""

from math import sqrt
import numpy as np

from compas.datastructures import Mesh
from compas.geometry import Polyline, Point, Line, Polygon
from compas.geometry import Rotation, Translation
from compas.geometry import is_coplanar
from compas.utilities import hex_to_rgb

try:
    from compas_view2 import app
    from compas_view2.collections import Collection
    from compas_view2.shapes import Arrow
except ImportError:
    pass


def hextorgb(hex):
    return tuple(i / 255 for i in hex_to_rgb(hex))


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
                fkeys = block.edge_faces(*edge)
                ps = [
                    block.face_center(fkeys[0]),
                    block.face_center(fkeys[1]),
                    *block.edge_coordinates(*edge),
                ]

                if is_coplanar(ps, tol=tol):
                    continue
            if assembly.graph.node_attribute(node, "is_support"):
                supportedges.append(Line(*block.edge_coordinates(*edge)))
            else:
                blockedges.append(Line(*block.edge_coordinates(*edge)))
    if len(blocks) != 0:
        viewer.add(
            Collection(blocks),
            show_faces=True,
            show_edges=False,
            opacity=0.6,
            facecolor=(0.9, 0.9, 0.9),
        )
    if len(supports) != 0:
        viewer.add(
            Collection(supports),
            show_faces=True,
            show_edges=False,
            opacity=0.5,
            facecolor=hextorgb("#f79d84"),
        )
    if len(blockedges) != 0:
        viewer.add(Collection(blockedges), linewidth=1.5)
    if len(supportedges) != 0:
        viewer.add(Collection(supportedges), linecolor=hextorgb("#f79d84"), linewidth=4)


def draw_interfaces(assembly, viewer):
    interfaces = []
    faces = []
    for edge in assembly.graph.edges():
        interface = assembly.graph.edge_attribute(edge, "interface")
        if interface is not None:
            corners = np.array(interface.points)
            faces.append(Polyline(np.vstack((corners, corners[0]))))
            if assembly.graph.node_attribute(
                edge[0], "is_support"
            ) or assembly.graph.node_attribute(edge[1], "is_support"):
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
        viewer.add(
            Collection(interfaces),
            show_edges=False,
            show_vertices=False,
            facecolor=(0.8, 0.8, 0.8),
        )
    if len(faces) != 0:
        viewer.add(
            Collection(faces),
            linecolor=hextorgb("#fac05e"),
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
        viewer.add(Collection(locs), size=12, color=hextorgb("#386641"))
    if len(res_np) != 0:
        viewer.add(Collection(res_np), linewidth=8, linecolor=(0, 0.3, 0))
    if len(res_nn) != 0:
        viewer.add(Collection(res_nn), linewidth=8, linecolor=(0.8, 0, 0))
    if len(fnn) != 0:
        viewer.add(Collection(fnn), linewidth=5, linecolor=hextorgb("#00468b"))
    if len(fnp) != 0:
        viewer.add(Collection(fnp), linewidth=5, linecolor=(1, 0, 0))
    if len(ft) != 0:
        viewer.add(Collection(ft), linewidth=5, linecolor=(1.0, 0.5, 0.0))


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
                is_tension = False
                for force in forces:
                    if force["c_np"] - force["c_nn"] <= -1e-5:
                        is_tension = True

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
                # print((w * sum_n + u * sum_u + v * sum_v).length * 100000, "edge: ", edge)

                # if assembly.graph.node_attribute(edge[0], 'is_support') or assembly.graph.node_attribute(edge[1], 'is_support'):
                #     print((w * sum_n + u * sum_u + v * sum_v).z)
                # total_reaction += abs((w * sum_n + u * sum_u + v * sum_v).z * 100000)

                p1 = resultant_pos + resultant_f
                p2 = resultant_pos - resultant_f

                if not is_tension:
                    res_np.append(Line(p1, p2))
                else:
                    res_nn.append(Line(p1, p2))
    if len(locs) != 0:
        viewer.add(Collection(locs), size=12, color=hextorgb("#386641"))
    if len(res_np) != 0:
        viewer.add(Collection(res_np), linewidth=8, linecolor=(0, 0.3, 0))
    if len(res_nn) != 0:
        viewer.add(Collection(res_nn), linewidth=8, linecolor=(0.8, 0, 0))
    if len(fnn) != 0:
        viewer.add(Collection(fnn), linewidth=5, linecolor=hextorgb("#00468b"))
    if len(fnp) != 0:
        viewer.add(Collection(fnp), linewidth=5, linecolor=(1, 0, 0))
    if len(ft) != 0:
        viewer.add(Collection(ft), linewidth=5, linecolor=(1.0, 0.5, 0.0))
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
        if assembly.graph.node_attribute(
            edge[0], "is_support"
        ) and not assembly.graph.node_attribute(edge[1], "is_support"):
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
                        f = Arrow(
                            pt,
                            w * force * scale * -1,
                            head_portion=0.2,
                            head_width=0.07,
                            body_width=0.02,
                        )
                    else:
                        f = Arrow(
                            pt,
                            w * force * scale,
                            head_portion=0.2,
                            head_width=0.07,
                            body_width=0.02,
                        )
                    if force >= 0:
                        fnp.append(f)
                    else:
                        fnn.append(f)
                    ft_uv = (u * forces[i]["c_u"] + v * forces[i]["c_v"]) * scale
                    if ft_uv.length == 0:
                        continue
                    if flip:
                        f = Arrow(
                            pt,
                            ft_uv * -1,
                            head_portion=0.2,
                            head_width=0.07,
                            body_width=0.02,
                        )
                    else:
                        f = Arrow(
                            pt,
                            ft_uv,
                            head_portion=0.2,
                            head_width=0.07,
                            body_width=0.02,
                        )
                    ft.append(f)
            if resultant:
                is_tension = False

                for force in forces:
                    if force["c_np"] - force["c_nn"] <= -1e-5:
                        is_tension = True

                sum_n = sum(force["c_np"] - force["c_nn"] for force in forces)
                sum_u = sum(force["c_u"] for force in forces)
                sum_v = sum(force["c_v"] for force in forces)
                if abs(sum_n) <= thres:
                    resultant_pos = np.average(
                        np.array(corners),
                        axis=0,
                        weights=[
                            sqrt(force["c_u"] ** 2 + force["c_v"] ** 2)
                            for force in forces
                        ],
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
                    f = Arrow(
                        resultant_pos,
                        resultant_f * -1,
                        head_portion=0.2,
                        head_width=0.07,
                        body_width=0.02,
                    )
                else:
                    f = Arrow(
                        resultant_pos,
                        resultant_f,
                        head_portion=0.2,
                        head_width=0.07,
                        body_width=0.02,
                    )
                if friction:
                    viewer.add(f, facecolor=(1.0, 0.5, 0.0), show_edges=False)
                if not is_tension:
                    res_np.append(f)
                else:
                    res_nn.append(f)
    if len(locs) != 0:
        viewer.add(Collection(locs), size=12, color=hextorgb("#386641"))
    if len(res_np) != 0:
        viewer.add(Collection(res_np), facecolor=hextorgb("#386641"), show_edges=False)
    if len(res_nn) != 0:
        viewer.add(Collection(res_nn), facecolor=(0.8, 0, 0), show_edges=False)
    if len(fnp) != 0:
        viewer.add(
            Collection(fnp),
            facecolor=hextorgb("#00468b"),
            show_edges=False,
            opacity=0.5,
        )
    if len(fnn) != 0:
        viewer.add(Collection(fnn), facecolor=(1, 0, 0), show_edges=False, opacity=0.5)
    if len(ft) != 0:
        viewer.add(
            Collection(ft), facecolor=(1.0, 0.5, 0.0), show_edges=False, opacity=0.5
        )


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
                fkeys = block.edge_faces(*edge)
                ps = [
                    block.face_center(fkeys[0]),
                    block.face_center(fkeys[1]),
                    *block.edge_coordinates(*edge),
                ]
                if is_coplanar(ps, tol=tol):
                    continue
            blocks.append(Line(*new_block.edge_coordinates(*edge)))
    if len(blocks) != 0:
        viewer.add(Collection(blocks), linewidth=1, linecolor=(0.7, 0.7, 0.7))
    if len(nodes) != 0:
        viewer.add(Collection(nodes), color=(0.7, 0.7, 0.7))


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
                head_portion=0.2,
                head_width=0.07,
                body_width=0.02,
            )
        )
        # print("self-weight", -block.volume() * density)
        # total_weights += block.volume() * 2500 * 9.8
        blocks.append(Point(*block.center()))

    # print("total self-weight: ", total_weights)

    if len(supports) != 0:
        viewer.add(Collection(supports), size=20, color=hextorgb("#ee6352"))
    if len(blocks) != 0:
        viewer.add(Collection(blocks), size=30, color=hextorgb("#3284a0"))
    if len(weights) != 0:
        viewer.add(Collection(weights), facecolor=hextorgb("#59cd90"), show_edges=False)


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

    viewer = app.App(width=1600, height=1000, viewmode="shaded", show_grid=grid)

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

    viewer.run()


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
    viewer : compas_view2.app.App
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


if __name__ == "__main__":
    pass
