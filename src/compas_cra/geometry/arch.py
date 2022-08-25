from math import radians
from compas.geometry import Rotation
from compas.geometry import add_vectors
from compas.geometry import subtract_vectors
from compas.geometry import transform_points
from compas.geometry import angle_vectors
from compas.datastructures import Mesh
from compas_assembly.datastructures import Block
from compas_cra.datastructures import CRA_Assembly


class Arch(object):
    """Create voussoir geometry for a semi-circular arch with given height and span.

    Parameters
    ----------
    height : float
        The distance between the base of the arch and the highest point of the intrados.
    span : float
        The distance between opposite intrados points at the base.
    thickness : float
        The distance between intrados and extrados.
    depth : float
        The depth of the arch.
    n : int
        Number of blocks
    """

    def __init__(
        self, height, span, thickness, depth, num_blocks=None, extra_support=False
    ):
        super().__init__()
        self.height = height
        self.span = span
        self.thickness = thickness
        self.depth = depth
        self.num_blocks = num_blocks
        self.extra_support = extra_support

    def assembly(self):
        """Create assembly.

        Returns
        -------
        assembly : compas_assembly.datastructures.Assembly
            Arch as assembly

        """

        assembly = CRA_Assembly()
        for mesh in self.blocks():
            assembly.add_block(mesh.copy(cls=Block))

        if self.extra_support is False:
            assembly.graph.node_attribute(0, "is_support", True)
            assembly.graph.node_attribute(self.num_blocks - 1, "is_support", True)
        else:
            assembly.graph.node_attribute(self.num_blocks, "is_support", True)
            assembly.graph.node_attribute(self.num_blocks + 1, "is_support", True)

        return assembly

    def blocks(self):
        """Create blocks.

        Returns
        -------
        list
            A list of blocks defined as simple meshes.

        """
        if self.height > self.span / 2:
            raise Exception("Not a semicircular arch.")

        radius = self.height / 2 + self.span**2 / (8 * self.height)
        top = [0.0, 0.0, self.height]
        left = [-self.span / 2, 0.0, 0.0]
        center = [0.0, 0.0, self.height - radius]
        vector = subtract_vectors(left, center)
        springing = angle_vectors(vector, [-1.0, 0.0, 0.0])
        sector = radians(180) - 2 * springing
        angle = sector / self.num_blocks

        a = top
        b = add_vectors(top, [0, self.depth, 0])
        c = add_vectors(top, [0, self.depth, self.thickness])
        d = add_vectors(top, [0, 0, self.thickness])

        R = Rotation.from_axis_and_angle([0, 1.0, 0], 0.5 * sector, center)
        bottom = transform_points([a, b, c, d], R)

        faces = [
            [0, 1, 2, 3],
            [7, 6, 5, 4],
            [3, 7, 4, 0],
            [6, 2, 1, 5],
            [7, 3, 2, 6],
            [5, 1, 0, 4],
        ]
        faces_inverse = [list(reversed(f)) for f in faces]

        R = Rotation.from_axis_and_angle([0, 1.0, 0], 0.5 * sector, center)
        bottom = transform_points([a, b, c, d], R)
        blocks = []
        for i in range(self.num_blocks):
            R = Rotation.from_axis_and_angle([0, 1.0, 0], -angle, center)
            top = transform_points(bottom, R)
            vertices = bottom + top
            faces = [
                [0, 1, 2, 3],
                [7, 6, 5, 4],
                [3, 7, 4, 0],
                [6, 2, 1, 5],
                [7, 3, 2, 6],
                [5, 1, 0, 4],
            ]
            mesh = Mesh.from_vertices_and_faces(vertices, faces)
            blocks.append(mesh)
            bottom = top

        if self.extra_support:
            a = [0.0, -self.thickness / 2, self.height - self.thickness / 2]
            b = add_vectors(a, [0, self.depth + self.thickness, 0])
            c = add_vectors(a, [0, self.depth + self.thickness, 2 * self.thickness])
            d = add_vectors(a, [0, 0, 2 * self.thickness])

            R_right = Rotation.from_axis_and_angle([0, 1.0, 0], 0.5 * sector, center)
            R_left = Rotation.from_axis_and_angle([0, 1.0, 0], -0.5 * sector, center)
            R_rightb = Rotation.from_axis_and_angle([0, 1.0, 0], 0.51 * sector, center)
            R_leftb = Rotation.from_axis_and_angle([0, 1.0, 0], -0.51 * sector, center)

            top_right_support = transform_points([a, b, c, d], R_right)
            bottom_right_support = transform_points([a, b, c, d], R_rightb)
            vertices_right_support = top_right_support + bottom_right_support
            block = Mesh.from_vertices_and_faces(vertices_right_support, faces_inverse)
            blocks.append(block)

            top_left_support = transform_points([a, b, c, d], R_left)
            bottom_left_support = transform_points([a, b, c, d], R_leftb)
            vertices_left_support = top_left_support + bottom_left_support
            block = Mesh.from_vertices_and_faces(vertices_left_support, faces)
            blocks.append(block)

        return blocks
