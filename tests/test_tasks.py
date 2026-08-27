"""The platform logic of the invoke tasks, mocked per-OS.

The Windows path is exercised for real on Windows dev machines; the macOS and Linux
branches cannot run here, but their environment construction is pure logic - these
tests pin it so a refactor cannot silently break a platform nobody is sitting at.
Skipped wherever the dev tooling is not installed (the CI test job installs only the
wheel and pytest).
"""

import os
import sys

import pytest

pytest.importorskip("invoke", reason="dev tooling not installed")
pytest.importorskip("compas_invocations2", reason="dev tooling not installed")

# tasks.py lives at the repository root, which the bare `pytest` binary does not put
# on sys.path (python -m pytest would)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tasks  # noqa: E402


def test_native_build_env_macos(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Darwin")
    monkeypatch.delenv("MACOSX_DEPLOYMENT_TARGET", raising=False)
    env = tasks._native_build_env()
    assert env == {"MACOSX_DEPLOYMENT_TARGET": "11.0"}


def test_native_build_env_macos_respects_override(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Darwin")
    monkeypatch.setenv("MACOSX_DEPLOYMENT_TARGET", "12.0")
    assert tasks._native_build_env()["MACOSX_DEPLOYMENT_TARGET"] == "12.0"


def test_native_build_env_linux_is_empty(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Linux")
    assert tasks._native_build_env() == {}


def test_native_build_env_windows(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Windows")
    monkeypatch.setattr(tasks, "_msys2_root", lambda: "C:\msys64")
    env = tasks._native_build_env()
    assert env["CMAKE_GENERATOR"] == "Ninja"
    assert env["PATH"].startswith(os.path.join("C:\msys64", "ucrt64", "bin"))


def test_verbosity_default_and_override(monkeypatch):
    monkeypatch.delenv("VERBOSITY", raising=False)
    assert tasks._verbosity() == "3"
    monkeypatch.setenv("VERBOSITY", "1")
    assert tasks._verbosity() == "1"
