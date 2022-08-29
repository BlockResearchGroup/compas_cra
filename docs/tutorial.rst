********************************************************************************
Tutorial
********************************************************************************


How to use CRA for your analysis
================================


.. _Creating_Geometry:

1. Creating geometries
----------------------

We create two blocks: one as support and another one as free block.

.. code-block:: python

    from compas.geometry import Box, Frame, Translation

    support = Box(Frame.worldXY(), 4, 2, 1)  # supporting block
    free1 = Box(
        Frame.worldXY().transformed(
            Translation.from_vector([0, 0, 1])
            * Rotation.from_axis_and_angle([0, 0, 1], 0.2)
        ), 1, 3, 1
    )  # block to analyse

2. CRA Assembly data structure
------------------------------

Add them to assembly data structure.

.. code-block:: python

    from compas_assembly.datastructures import Block
    from compas_cra.datastructures import CRA_Assembly

    assembly = CRA_Assembly()
    assembly.add_block(Block.from_shape(support), node=0)
    assembly.add_block(Block.from_shape(free1), node=1)

.. figure:: /_images/tutorial_cubes_1.png
    :figclass: figure
    :class: figure-img img-fluid


.. _BC:

3. Boundary conditions
----------------------

Set boundary conditions.

.. code-block:: python

    assembly.set_boundary_conditions([0])

.. figure:: /_images/tutorial_cubes_2.png
    :figclass: figure
    :class: figure-img img-fluid

4. Identifying interfaces
-------------------------
Then we identify planar interfaces between blocks automatically.

.. code-block:: python

    from compas_cra.algorithms import assembly_interfaces_numpy
    assembly_interfaces_numpy(assembly)

.. figure:: /_images/tutorial_cubes_3.png
    :figclass: figure
    :class: figure-img img-fluid

5. Solving equilibrium
----------------------

:mod:`compas_cra` provide three solvers:

- RBE Solve: :mod:`compas_cra.equilibrium.rbe_solve`.
- CRA Solve: :mod:`compas_cra.equilibrium.cra_solve`.
- CRA Penalty Solve: :mod:`compas_cra.equilibrium.cra_penalty_solve`.

.. code-block:: python

    from compas_cra.equilibrium import cra_solve
    cra_solve(assembly, verbose=True, timer=True)

6. Visualisation
----------------

.. code-block:: python

    from compas_cra.viewers import cra_view
    cra_view(assembly, resultant=False, nodal=True, grid=True)

.. figure:: /_images/tutorial_cubes_4.png
    :figclass: figure
    :class: figure-img img-fluid

The complete tutorial script can be downloaded from
`scripts/tutorial_cubes.py <https://github.com/BlockResearchGroup/compas_cra/blob/main/scripts/tutorial_cubes.py>`_

To see more how to construct assembly and solve equilibrium, please check :ref:`Examples`.


Optional: Export geometry from CAD software (Rhino)
=========================================

For the step :ref:`Creating_Geometry`, we can also input geometry from CAD software.
Here we use Rhino as an example.


Export mesh blocks as Assembly json file
----------------------------------------

Use this script at `scripts/mesh_to_assembly_json.py <https://github.com/BlockResearchGroup/compas_cra/blob/main/scripts/mesh_to_assembly_json.py>`_
to select Rhino mesh blocks and export to assembly data structure as a ``.json`` file.

.. literalinclude:: ../scripts/mesh_to_assembly_json.py
    :language: python

**Note**: The selection sequence is important because it represents the node indices.

Then we can load the ``.json`` file from local path.

.. code-block:: python

    import os
    import compas
    import compas_cra

    from compas_cra.datastructures import CRA_Assembly

    FILE_I = os.path.join(compas_cra.DATA, "XXX.json")  # or your own path
    assembly = compas.json_load(FILE_I)
    assembly = assembly.copy(cls=CRA_Assembly)

After loading the ``.json`` file and convert it to :mod:`compas_cra.datastructures.CRA_Assembly`,
we can follow the previous step :ref:`BC` for the analysis.

Export mesh blocks and interfaces as Assembly json file
-------------------------------------------------------

Currently, we do not implement automatic interface detection algorithm for blocks with curve/free-form interfaces,
so they have to be discretised manually as planar faces or triangles.

Use this script at `scripts/mesh_to_assembly_interfaces_json.py <https://github.com/BlockResearchGroup/compas_cra/blob/main/scripts/mesh_to_assembly_interfaces_json.py>`_
to select mesh blocks with interfaces and export to assembly data structure and store it as json file.

.. literalinclude:: ../scripts/mesh_to_assembly_interfaces_json.py
    :language: python

**Note**:

- The selection sequence is important because it represents the node indices.
- The interface normal directions are important, it must point from **assign interface from** to **assign interface to**. For example, in Figure 1, **assign interface from: 1** and **assign interface to: 0** for the interface from 1 to 0.


.. list-table:: Figure 1
   :class: borderless

   * - .. image:: /_images/tutorial_interface1.png
          :width: 100 %
     - .. image:: /_images/tutorial_interface2.png
          :width: 100 %

More Rhino files and precomputed :code:`.json` files are located at `data <https://github.com/BlockResearchGroup/compas_cra/blob/main/data>`_ folder.

**Note**:

Every time a new file is opened in Rhino, be sure to restart Rhino or reset the Python Script Engine before running scripts.
