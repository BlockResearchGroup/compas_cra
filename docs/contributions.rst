********************************************************************************
How to contribute
********************************************************************************

.. highlight:: bash

Contributions are welcome and very much appreciated!

Code contributions
==================

We accept code contributions through pull requests.
In short, this is how that works.

1. Fork `the repository <https://github.com/petrasvestartas/compas_cra>`_ and clone the fork.
2. Create a virtual environment using your tool of choice (e.g. ``virtualenv``, ``conda``, etc).
3. Build the solver once — it is compiled into the package, so an install from source
   needs a staged IPOPT tree. See :ref:`Installation` for the compilers and libraries
   this needs on your platform, and for the Docker route if you would rather not
   install them:

::

    $ packaging/build_ipopt.sh

4. Install the package and its development dependencies:

::

    $ pip install -e ".[dev]"


5. Make sure all tests pass:

::

    $ invoke test


6. Start making your changes to the **main** branch (or branch off of it).
7. Make sure all tests still pass:

::

    $ invoke test


8. Add yourself **Contributors** section to ``AUTHORS.md``.
9. Commit your changes and push your branch to GitHub.
10. Create a `pull request <https://help.github.com/articles/about-pull-requests/>`_ through the GitHub website.

During development, use `pyinvoke <http://docs.pyinvoke.org/>`_ tasks on the
command line to ease recurring operations:

* ``invoke clean``: Clean all generated artifacts.
* ``invoke check``: Run various code and documentation style checks.
* ``invoke docs``: Generate documentation.
* ``invoke test``: Run all tests and checks in one swift command.
* ``invoke``: Show available tasks.

We highly recommend you have `Black <https://black.readthedocs.io/en/stable/index.html>`_ installed in your favourite editor/IDE.
For setting up Black in your editor, please refer to their documentation `Editor integration <https://black.readthedocs.io/en/stable/integrations/editors.html>`_.

Bug reports
===========

When `reporting a bug <https://github.com/petrasvestartas/compas_cra/issues>`_
please include:

* Operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Feature requests and feedback
=============================

The best way to send feedback is to file an issue on
`Github <https://github.com/petrasvestartas/compas_cra/issues>`_.
If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.


