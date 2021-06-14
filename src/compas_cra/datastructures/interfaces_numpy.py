#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Identify interfaces for CRA extended assembly data structures.
"""
from math import fabs

from numpy import array
from numpy import float64
from scipy.linalg import solve
from scipy.spatial import cKDTree
from shapely.geometry import Polygon

from compas.geometry import Frame
from compas.geometry import local_to_world_coordinates_numpy

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['assembly_interfaces_numpy']


def _find_nearest_neighbours(cloud, nmax):
    tree = cKDTree(cloud)
    nnbrs = [tree.query(root, nmax) for root in cloud]
    nnbrs = [(d.flatten().tolist(), n.flatten().tolist()) for d, n in nnbrs]
    return nnbrs


def assembly_interfaces_numpy(assembly,
                              nmax=10,
                              tmax=1e-6,
                              amin=1e-1,
                              concave=True):
    """Identify the interfaces between the blocks of an assembly.

    Parameters
    ----------
    assembly : compas_cra.datastructures.CRA_Assembly
        An assembly of discrete blocks.
    nmax : int, optional
        Maximum number of neighbours per block.
        Default is ``10``.
    tmax : float, optional
        Maximum deviation from the perfectly flat interface plane.
        Default is ``1e-6``.
    amin : float, optional
        Minimum area of a "face-face" interface.
        Default is ``1e-1``.
    concave : bool, optional
        Identify concave interfaces
        Default is ``True``

    References
    ----------

    Examples
    --------
    .. code-block:: python

        pass

    """
    node_index = assembly.key_index()
    index_node = assembly.index_key()

    blocks = assembly.nodes_attribute('block')

    nmax = min(nmax, len(blocks))

    block_cloud = assembly.nodes_attributes('xyz')
    block_nnbrs = _find_nearest_neighbours(block_cloud, nmax)

    # k:      key of the base block
    # i:      index of the base block
    # block:  base block
    # nbrs:   list of indices of the neighbouring blocks
    # frames: list of frames for each of the faces of the base block

    # f0:   key of the current base face
    # A:    uvw base frame of f0
    # o:    origin of the base frame of f0
    # xyz0: xyz coordinates of the nodes of f0
    # rst0: local coordinates of the nodes of f0, with respect to the frame f0
    # p0:   2D polygon of f0 in local coordinates

    # j:   index of the current neighbour
    # n:   key of the current neighbour
    # nbr: neighbour block
    # k_i: key index map for the nodes of the nbr block
    # xyz: xyz coorindates of all nodes of nbr
    # rst: local coordinates of all nodes of nbr, with respect to the frame f0

    # f1:   key of the current neighbour face
    # rst1: local coordinates of the nodes of f1, with respect to the frame f0
    # p1:   2D polygon of f1 in local coordinates

    count = 0
    for node in assembly.nodes():

        i = node_index[node]

        block = blocks[i]
        nbrs = block_nnbrs[i][1]

        frames = block.frames()

        for f0, (origin, uvw) in frames.items():
            A = array(uvw, dtype=float64)
            o = array(origin, dtype=float64).reshape((-1, 1))
            xyz0 = array(block.face_coordinates(f0),
                         dtype=float64).reshape((-1, 3)).T
            rst0 = solve(A.T, xyz0 - o).T.tolist()
            p0 = Polygon(rst0)

            for j in nbrs:
                n = index_node[j]

                if n == node:
                    continue

                if not concave:
                    if node in assembly.edge and n in assembly.edge[node]:
                        continue

                if n in assembly.edge and node in assembly.edge[n]:
                    continue

                if assembly.node_attribute(node, 'is_support') and \
                   assembly.node_attribute(n, 'is_support'):
                    continue

                nbr = blocks[j]
                k_i = {key: index
                       for index, key in enumerate(nbr.vertices())}
                xyz = array(nbr.vertices_attributes('xyz'),
                            dtype=float64).reshape((-1, 3)).T
                rst = solve(A.T, xyz - o).T.tolist()
                rst = {key: rst[k_i[key]] for key in nbr.vertices()}

                for f1 in nbr.faces():

                    rst1 = [rst[key] for key in nbr.face_vertices(f1)]

                    if any(fabs(t) > tmax for r, s, t in rst1):
                        continue
                    p1 = Polygon(rst1)

                    if p1.area < amin:
                        continue

                    if p0.intersects(p1):
                        intersection = p0.intersection(p1)
                        area = intersection.area
                        if not concave and count >= 1:
                            break

                        if area >= amin:
                            coords = [
                                [x, y, 0.0]
                                for x, y, z in
                                intersection.exterior.coords]

                            coords = local_to_world_coordinates_numpy(
                                Frame(o, A[0], A[1]), coords)

                            if concave:
                                assembly.add_to_interfaces(
                                    node, n,
                                    itype='face_face',
                                    isize=area,
                                    ipoints=coords.tolist()[:-1],
                                    iframe=Frame(origin, uvw[0], uvw[1]))
                            else:
                                assembly.add_interface(
                                    node, n,
                                    itype='face_face',
                                    isize=area,
                                    ipoints=coords.tolist()[:-1],
                                    iframe=Frame(origin, uvw[0], uvw[1]))

                            count += 1
    return assembly


if __name__ == "__main__":
    pass
