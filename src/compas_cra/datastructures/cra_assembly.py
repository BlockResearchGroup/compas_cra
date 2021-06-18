#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CRA assembly data structures
"""

from compas.geometry import Frame
from compas.geometry import Rotation
from compas.geometry import Pointcloud
from compas_assembly.datastructures import Assembly
from compas_assembly.datastructures import Interface

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['CRA_Assembly']


class CRA_Assembly(Assembly):
    """Extended data structure for concave assemblies."""

    def __init__(self):
        super(CRA_Assembly, self).__init__()
        self.attributes.update({'name': 'CRA_Assembly'})
        self.default_node_attributes.update({
            'block': None,
            'displacement': [0, 0, 0, 0, 0, 0]
        })
        self.default_edge_attributes.update({
            'interface': None,
            'interfaces': []
        })

    def add_to_interfaces(self, u, v, itype, isize, ipoints, iframe):
        interface = Interface(itype=itype, isize=isize, ipoints=ipoints,
                              iframe=iframe)
        self.add_interface_to_interfaces(u, v, interface)

    def add_interface_to_interfaces(self, u, v, interface):
        if not self.has_edge(u, v):
            self.add_edge(u, v, interfaces=[interface])
        else:
            interfaces = self.edge_attribute((u, v), "interfaces")
            interfaces.append(interface)
            self.edge_attribute((u, v), "interfaces", interfaces)

    def add_interfaces_from_meshes(self, meshes, u, v):
        for mesh in meshes:
            for f in mesh.faces():
                pt = mesh.face_coordinates(f)
                interface = Interface(itype='face_face',
                                      isize=mesh.face_area(f),
                                      ipoints=pt,
                                      iframe=Frame.from_points(pt[0],
                                                               pt[1], pt[2]))
                self.add_interface_to_interfaces(u, v, interface)

    def set_boundary_conditions(self, keys):
        for key in keys:
            self.set_boundary_condition(key)

    def set_boundary_condition(self, key):
        self.node_attribute(key, "is_support", True)

    def is_block_support(self, key):
        return self.node_attribute(key, "is_support")

    def rotate_assembly(self, o, axis, rad):
        R = Rotation().from_axis_and_angle(axis, angle=rad, point=o)
        self.transform(R)
        for edge in self.edges():
            for interface in self.edge_attribute(edge, "interfaces"):
                interface.points = [list(c) for c in
                                    Pointcloud(interface.points).transformed(
                                        R)]
                interface.frame.transform(R)

    def get_weight_total(self, density=1):
        weight = 0
        for node in self.nodes():
            block = self.node_attribute(node, 'block')
            weight += block.volume() * density
        return weight

    def get_weight_mean(self, density=1):
        n = self.number_of_nodes()
        w = self.get_weight_total(density)
        return w / n


if __name__ == '__main__':
    pass
