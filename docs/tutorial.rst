********************************************************************************
Tutorial
********************************************************************************


Export geometry from CAD software (Rhino)
=========================================

Rhino Installation
-------------------

Please make sure install :mod:`compas_cra` first, see :ref:`Installation`.

:mod:`compas_cra` is developed independent of the functionality of CAD software.
However, CAD software is still useful to create complicated geometrical objects.
For a more detailed information on how to install COMPAS and its packages for Rhino,
please refer to `Working in Rhino <https://compas.dev/compas/latest/gettingstarted/rhino.html>`_ page of the COMPAS documentation.

In order to install :mod:`compas_cra` for Rhino, do

::

    $ python -m compas_rhino.uninstall
    $ python -m compas_rhino.install
    $ python -m compas_rhino.install -p compas_cra

Every time a new file is opened in Rhino, be sure to restart Rhino or reset the Python Script Engine before running scripts.


Export mesh blocks as Assembly json file
----------------------------------------

Use this script at `scripts/mesh_to_assembly_json.py <https://github.com/BlockResearchGroup/compas_cra/blob/main/scripts/mesh_to_assembly_json.py>`_
to select mesh blocks and export to assembly data structure and store it as json file.

.. literalinclude:: ../scripts/mesh_to_assembly_json.py
    :language: python

**Note**: The selection sequence is important because it represents the node indices.


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

Using CRA Solvers
=========================================

:mod:`compas_cra` provide three solvers:

- RBE Solve: :mod:`compas_cra.equilibrium.rbe_solve`.
- CRA Solve: :mod:`compas_cra.equilibrium.cra_solve`.
- CRA Penalty Solve: :mod:`compas_cra.equilibrium.cra_penalty_solve`.

To see more how to construct assembly and solve equilibrium, please check :ref:`Examples`.
