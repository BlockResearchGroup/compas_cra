# Known issues

- `compas_cra` uses the [IPOPT](https://coin-or.github.io/Ipopt/) solver with
  the MUMPS linear solver. Ill-conditioned assemblies can end in a non-optimal
  termination condition, which is reported as a `ValueError` by the solver
  functions.

- The solver is a compiled extension inside the package (`compas_cra._native`),
  shipped in the platform wheels. A source checkout has no solver until it is built:
  run `packaging/build_ipopt.sh` and then `pip install .` — see
  `native/README.md`. Without it the solver functions raise
  `no NLP backend available`.
