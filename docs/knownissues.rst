********************************************************************************
Known Issues
********************************************************************************

- pyomo cannot find ipopt location: ``pyomo.common.errors.ApplicationError: No executable found for solver 'ipopt'``.
  Please refer this thread for solution: https://groups.google.com/g/open-dsopf/c/wYPbZp-HLCw?pli=1

- :mod:`compas_cra` uses `IPOPT <https://coin-or.github.io/Ipopt/>`_ solver, so it might not work for PC with AMD processor.
