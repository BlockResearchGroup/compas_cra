"""CRA assembly data structures"""

from compas.geometry import Frame
from compas.geometry import Rotation
from compas.geometry import Pointcloud
from compas_assembly.datastructures import Assembly
from compas_assembly.datastructures import Interface
from compas_assembly.datastructures import Block


class CRA_Assembly(Assembly):
    """Extended data structure for concave assemblies."""

    def __init__(self):
        # keep this old coding style for Rhino IronPython. super().__init__() doesn't work
        super(CRA_Assembly, self).__init__()

        self.attributes.update({"name": "CRA_Assembly"})
        self.graph.default_node_attributes.update(
            {"block": None, "displacement": [0, 0, 0, 0, 0, 0]}
        )
        self.graph.default_edge_attributes.update({"interface": None, "interfaces": []})

    def add_blocks_from_rhinomeshes(self, guids):
        """Add multiple blocks from their representation as as Rhino meshes.

        Parameters
        ----------
        guids : list of str
            A list of GUIDs identifying the meshes representing the blocks of the assembly.

        Returns
        -------
        list
            The keys of the added blocks.

        """
        keys = []
        for guid in guids:
            block = Block.from_rhinomesh(guid)
            key = self.add_block(block)
            keys.append(key)
        return keys

    def add_to_interfaces(self, u, v, type, size, points, frame):
        """Add interface to edge (u, v) interfaces.

        Parameters
        ----------
        u : int
            block_j id.
        v : int
            block_k id.
        type : str
            Interface type.
        size : float
            Interface area.
        points : int
            Interface vertices.
        frame : int
            Local coordinate.

        Returns
        -------
        None

        """
        interface = Interface(type=type, size=size, points=points, frame=frame)
        self.add_interface_to_interfaces(u, v, interface)

    def add_interface_to_interfaces(self, u, v, interface):
        """Add interface to edge (u, v) interfaces.

        Parameters
        ----------
        u : int
            block_j id.
        v : int
            block_k id.
        interface : class:`compas_assembly.datastructures.Interface`
            Interface.

        Returns
        -------
        None

        """
        if not self.graph.has_edge(u, v):
            self.graph.add_edge(u, v, interfaces=[interface])
        else:
            interfaces = self.graph.edge_attribute((u, v), "interfaces")
            interfaces.append(interface)
            self.graph.edge_attribute((u, v), "interfaces", interfaces)

    def add_interfaces_from_meshes(self, meshes, u, v):
        """Add interfaces from meshes to edge (u, v) interfaces.

        Parameters
        ----------
        meshes : list of compas.datastructures.Mesh
            Meshes.
        u : int
            block_j id.
        v : int
            block_k id.

        Returns
        -------
        None

        """
        for mesh in meshes:
            for f in mesh.faces():
                pt = mesh.face_coordinates(f)
                interface = Interface(
                    type="face_face",
                    size=mesh.face_area(f),
                    points=pt,
                    frame=Frame.from_points(pt[0], pt[1], pt[2]),
                )
                self.add_interface_to_interfaces(u, v, interface)

    def set_boundary_conditions(self, keys):
        """Set blocks as boundary conditions.

        Parameters
        ----------
        keys : list of int
            Assembly node keys.

        Returns
        -------
        None

        """
        for key in keys:
            self.set_boundary_condition(key)

    def set_boundary_condition(self, key):
        """Set block as boundary condition.

        Parameters
        ----------
        key : int
            Assembly node key.

        Returns
        -------
        None

        """
        self.graph.node_attribute(key, "is_support", True)

    def is_block_support(self, key):
        """Check if the block is a support.

        Parameters
        ----------
        key : int
            Assembly node key.

        Returns
        -------
        None

        """
        return self.graph.node_attribute(key, "is_support")

    def rotate_assembly(self, o, axis, angle, is_rad=False):
        """Rotate the entire assembly.

        Parameters
        ----------
        o : list
            Rotation origin.
        axis : list
            Rotation axis.
        angle : float
            Rotation angle.
        is_rad : bool, optional
            True: angle is radian. False: angle is degree.

        Returns
        -------
        None

        """
        from math import pi

        rad = angle
        if not is_rad:
            rad = angle * pi / 180

        R = Rotation().from_axis_and_angle(axis, angle=rad, point=o)
        self.transform(R)
        for edge in self.edges():
            for interface in self.graph.edge_attribute(edge, "interfaces"):
                interface.points = [
                    list(c) for c in Pointcloud(interface.points).transformed(R)
                ]
                interface.frame.transform(R)

    def move_block(self, key, vector=(0, 0, 0)):
        """Move block with vector.

        Parameters
        ----------
        key : int
            Assembly node key.
        vector : list, optional
            Translation vector.

        Returns
        -------
        None

        """
        from compas.geometry import Translation

        self.graph.node_attribute(key, "block").transform(
            Translation.from_vector(vector)
        )

    def get_weight_total(self, density=1):
        """Get total assembly weight.

        Parameters
        ----------
        density : float, optional
            Material density.

        Returns
        -------
        None

        """
        weight = 0
        for node in self.nodes():
            block = self.graph.node_attribute(node, "block")
            weight += block.volume() * density
        return weight

    def get_weight_mean(self, density=1):
        """Get assembly mean weight.

        Parameters
        ----------
        density : float, optional
            Material density.

        Returns
        -------
        None

        """
        n = self.graph.number_of_nodes()
        w = self.get_weight_total(density)
        return w / n


if __name__ == "__main__":
    pass
