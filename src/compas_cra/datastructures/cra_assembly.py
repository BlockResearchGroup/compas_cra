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

    # TODO: make density as default attribute for the assembly class member

    def __init__(self):
        super(CRA_Assembly, self).__init__()
        self.attributes.update({'name': 'CRA_Assembly'})
        self.graph.default_node_attributes.update({
            'block': None,
            'displacement': [0, 0, 0, 0, 0, 0]
        })
        self.graph.default_edge_attributes.update({
            'interface': None,
            'interfaces': []
        })

    def add_to_interfaces(self, u, v, type, size, points, frame):
        """Add interface from attributes to edge (u, v) interfaces"""
        interface = Interface(type=type, size=size, points=points,
                              frame=frame)
        self.add_interface_to_interfaces(u, v, interface)

    def add_interface_to_interfaces(self, u, v, interface):
        """Add interface to edge (u, v) interfaces"""
        if not self.graph.has_edge(u, v):
            self.graph.add_edge(u, v, interfaces=[interface])
        else:
            interfaces = self.graph.edge_attribute((u, v), "interfaces")
            interfaces.append(interface)
            self.graph.edge_attribute((u, v), "interfaces", interfaces)

    def add_interfaces_from_meshes(self, meshes, u, v):
        """Add interfaces from meshes to edge (u, v) interfaces"""
        for mesh in meshes:
            for f in mesh.faces():
                pt = mesh.face_coordinates(f)
                interface = Interface(itype='face_face',
                                      isize=mesh.face_area(f),
                                      ipoints=pt,
                                      iframe=Frame.from_points(pt[0], pt[1], pt[2]))
                self.add_interface_to_interfaces(u, v, interface)

    def set_boundary_conditions(self, keys):
        """Set blocks as boundary conditions"""
        for key in keys:
            self.set_boundary_condition(key)

    def set_boundary_condition(self, key):
        """Set block as boundary condition"""
        self.graph.node_attribute(key, "is_support", True)

    def is_block_support(self, key):
        """Check if the block is a support"""
        return self.graph.node_attribute(key, "is_support")

    def rotate_assembly(self, o, axis, rad):
        """Rotate the entire assembly"""
        R = Rotation().from_axis_and_angle(axis, angle=rad, point=o)
        self.transform(R)
        for edge in self.edges():
            for interface in self.graph.edge_attribute(edge, "interfaces"):
                interface.points = [list(c) for c in
                                    Pointcloud(interface.points).transformed(
                                        R)]
                interface.frame.transform(R)

    def move_block(self, key, vector=(0, 0, 0)):
        """Move block with vector"""
        from compas.geometry import Translation

        self.graph.node_attribute(key, "block").transform(
            Translation.from_vector(vector))

    def get_weight_total(self, density=1):
        """Get total assembly weight"""
        weight = 0
        for node in self.nodes():
            block = self.graph.node_attribute(node, 'block')
            weight += block.volume() * density
        return weight

    def get_weight_mean(self, density=1):
        """Get assembly mean weight"""
        n = self.graph.number_of_nodes()
        w = self.get_weight_total(density)
        return w / n


if __name__ == '__main__':
    pass
