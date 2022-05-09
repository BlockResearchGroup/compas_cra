def test_mosek():
    import mosek
    import pyomo.environ as pyo

    model = pyo.ConcreteModel()
    model.x = pyo.Var([1, 2], domain=pyo.NonNegativeReals)
    model.OBJ = pyo.Objective(expr=2 * model.x[1] + 3 * model.x[2])
    model.Constraint1 = pyo.Constraint(expr=3 * model.x[1] + 4 * model.x[2] >= 1)

    with pyo.SolverFactory("mosek") as solver:
        # options - MOSEK parameters dictionary, using strings as keys (optional)
        # tee - write log output if True (optional)
        # soltype - accepts three values : bas, itr and itg for basic,
        # interior point and integer solution, respectively. (optional)
        result = solver.solve(
            model,
            options={'dparam.optimizer_max_time': 100.0,
                     'iparam.intpnt_solve_form': int(mosek.solveform.dual)},
            tee=True, soltype='itr')

        assert result.Solver._active

