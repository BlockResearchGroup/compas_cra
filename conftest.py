import pytest
import compas
import compas_cra
import math
import numpy


def pytest_ignore_collect(path):
    if "rhino" in str(path):
        return True

    if "blender" in str(path):
        return True

    if "ghpython" in str(path):
        return True


@pytest.fixture(autouse=True)
def add_compas(doctest_namespace):
    doctest_namespace["compas"] = compas


@pytest.fixture(autouse=True)
def add_compas_cra(doctest_namespace):
    doctest_namespace["compas_cra"] = compas_cra


@pytest.fixture(autouse=True)
def add_math(doctest_namespace):
    doctest_namespace["math"] = math


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace["np"] = numpy
