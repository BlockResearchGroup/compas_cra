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
    monkeypatch.setattr(tasks, "_msys2_root", lambda: r"C:\msys64")
    env = tasks._native_build_env()
    assert env["CMAKE_GENERATOR"] == "Ninja"
    assert env["PATH"].startswith(os.path.join(r"C:\msys64", "ucrt64", "bin"))


def test_verbosity_default_and_override(monkeypatch):
    monkeypatch.delenv("VERBOSITY", raising=False)
    assert tasks._verbosity() == "3"
    monkeypatch.setenv("VERBOSITY", "1")
    assert tasks._verbosity() == "1"


def test_brew_prefers_path(monkeypatch):
    monkeypatch.setattr(tasks.shutil, "which", lambda name: "/opt/homebrew/bin/brew")
    assert tasks._brew() == "/opt/homebrew/bin/brew"


def test_brew_falls_back_to_prefixes(monkeypatch):
    """A brew installed moments ago is on disk before it is ever on PATH."""
    monkeypatch.setattr(tasks.shutil, "which", lambda name: None)
    monkeypatch.setattr(tasks.os, "access", lambda path, mode: path == "/usr/local/bin/brew")
    assert tasks._brew() == "/usr/local/bin/brew"


def test_brew_returns_none_when_absent(monkeypatch):
    monkeypatch.setattr(tasks.shutil, "which", lambda name: None)
    monkeypatch.setattr(tasks.os, "access", lambda path, mode: False)
    assert tasks._brew() is None


def test_use_brew_prepends_once(monkeypatch):
    monkeypatch.setenv("PATH", "/usr/bin")
    tasks._use_brew("/opt/homebrew/bin/brew")
    assert os.environ["PATH"].split(os.pathsep)[0] == "/opt/homebrew/bin"
    tasks._use_brew("/opt/homebrew/bin/brew")
    assert os.environ["PATH"].count("/opt/homebrew/bin") == 1


def test_install_homebrew_refused_by_env(monkeypatch):
    """The escape hatch fails instead of installing anything."""
    monkeypatch.setenv("COMPAS_CRA_NO_BREW_INSTALL", "1")
    with pytest.raises(tasks.Exit):
        tasks._install_homebrew(None)


def test_patch_venv_activation_is_a_noop_on_linux(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Linux")
    tasks._patch_venv_activation()  # must not raise, must not look for a venv


def test_append_once_is_idempotent(tmp_path):
    script = tmp_path / "activate"
    script.write_text("# original\n", encoding="utf-8")
    tasks._append_once(str(script), [tasks.VENV_PATCH_MARK_MACOS, "export FOO=1"])
    tasks._append_once(str(script), [tasks.VENV_PATCH_MARK_MACOS, "export FOO=1"])
    assert script.read_text(encoding="utf-8").count("export FOO=1") == 1


def test_append_once_skips_a_missing_file(tmp_path):
    tasks._append_once(str(tmp_path / "nope"), ["export FOO=1"])


class _FakeResult:
    def __init__(self, ok):
        self.ok = ok


class _FakeCtx:
    """Records the commands a task would run, and answers each with a fixed result."""

    def __init__(self, results):
        self.results = list(results)
        self.commands = []

    def run(self, command, **kwargs):
        self.commands.append(command)
        return _FakeResult(self.results.pop(0))


def test_refuse_root_allows_a_normal_user(monkeypatch):
    monkeypatch.setattr(tasks.os, "geteuid", lambda: 501, raising=False)
    tasks._refuse_root()


def test_refuse_root_rejects_sudo(monkeypatch):
    monkeypatch.setattr(tasks.os, "geteuid", lambda: 0, raising=False)
    with pytest.raises(tasks.Exit):
        tasks._refuse_root()


def test_install_homebrew_pairs_sudo_with_the_installer(monkeypatch):
    """One ctx.run, so `sudo -v` and the installer share a pty - and a credential."""
    monkeypatch.delenv("COMPAS_CRA_NO_BREW_INSTALL", raising=False)
    monkeypatch.setattr(tasks, "_brew", lambda: "/opt/homebrew/bin/brew")
    ctx = _FakeCtx([True])
    assert tasks._install_homebrew(ctx) == "/opt/homebrew/bin/brew"
    assert len(ctx.commands) == 1
    assert ctx.commands[0].startswith("sudo -v && ")
    assert tasks.HOMEBREW_INSTALL_URL in ctx.commands[0]


def test_install_homebrew_retries_interactively(monkeypatch):
    """A failed non-interactive run falls back to letting the installer prompt."""
    monkeypatch.delenv("COMPAS_CRA_NO_BREW_INSTALL", raising=False)
    brews = iter([None, "/opt/homebrew/bin/brew"])
    monkeypatch.setattr(tasks, "_brew", lambda: next(brews))
    ctx = _FakeCtx([False, True])
    assert tasks._install_homebrew(ctx) == "/opt/homebrew/bin/brew"
    assert len(ctx.commands) == 2
    assert not ctx.commands[1].startswith("sudo -v")


def test_install_homebrew_exits_when_no_brew_appears(monkeypatch):
    monkeypatch.delenv("COMPAS_CRA_NO_BREW_INSTALL", raising=False)
    monkeypatch.setattr(tasks, "_brew", lambda: None)
    with pytest.raises(tasks.Exit):
        tasks._install_homebrew(_FakeCtx([False, False]))


def test_ipopt_work_dir_stays_in_tree_without_spaces():
    base = os.path.join(os.sep, "home", "dev", "compas_cra")
    assert tasks._ipopt_work_dir(base) == os.path.join(base, "build", "ipopt")
    assert tasks._ipopt_prefix_override(base) is None


def test_ipopt_work_dir_redirects_a_spaced_path(monkeypatch, capsys):
    """Autotools cannot build under whitespace, so the build leaves the tree."""
    monkeypatch.setattr(tasks.os.path, "expanduser", lambda p: p.replace("~", os.path.join(os.sep, "home", "dev")))
    base = os.path.join(os.sep, "home", "dev", "untitled folder", "compas_cra")
    expected = os.path.join(os.sep, "home", "dev", ".compas_cra", "ipopt")
    assert tasks._ipopt_work_dir(base) == expected
    assert tasks._ipopt_prefix_override(base) == os.path.join(expected, "stage")


def test_ipopt_work_dir_exits_when_home_is_spaced_too(monkeypatch):
    monkeypatch.setattr(
        tasks.os.path, "expanduser", lambda p: p.replace("~", os.path.join(os.sep, "home", "spaced name"))
    )
    base = os.path.join(os.sep, "srv", "untitled folder", "compas_cra")
    with pytest.raises(tasks.Exit):
        tasks._ipopt_work_dir(base)


def test_msys_path():
    assert tasks._msys_path(r"C:\Users\dev\.compas_cra\ipopt") == "/c/Users/dev/.compas_cra/ipopt"


def test_codesign_extension_signs_each_so(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(tasks.sysconfig, "get_path", lambda name: os.path.join(os.sep, "venv", "site-packages"))
    sos = [os.path.join(os.sep, "venv", "site-packages", "compas_cra", "_native", "_core.so")]
    monkeypatch.setattr(tasks.glob, "glob", lambda pattern, recursive: list(sos))
    ctx = _FakeCtx([True])
    tasks._codesign_extension(ctx)
    assert ctx.commands == ['codesign --force --sign - "%s"' % sos[0]]


def test_codesign_extension_is_a_noop_off_macos(monkeypatch):
    monkeypatch.setattr(tasks.platform, "system", lambda: "Windows")
    tasks._codesign_extension(None)  # must not touch ctx at all
