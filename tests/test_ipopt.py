import pyomo.environ as pyo
import pytest

from compas_cra._ipopt import executable
from compas_cra._ipopt import ipopt_version
from compas_cra.equilibrium._solver import ipopt_solver


def test_ipopt_executable_is_available():
    exe = executable()
    assert exe is not None, "no ipopt executable found"
    assert exe.is_file()
    assert "Ipopt" in (ipopt_version() or "")


def test_ipopt():
    model = pyo.ConcreteModel()
    model.x = pyo.Var([1, 2], domain=pyo.NonNegativeReals)
    model.OBJ = pyo.Objective(expr=2 * model.x[1] + 3 * model.x[2])
    model.Constraint1 = pyo.Constraint(expr=3 * model.x[1] + 4 * model.x[2] >= 1)
    with ipopt_solver() as solver:
        result = solver.solve(model, tee=True)
        assert result.Solver._active
        assert result.solver.termination_condition == pyo.TerminationCondition.optimal
        assert pyo.value(model.OBJ) == pytest.approx(2 / 3, rel=1e-6)
