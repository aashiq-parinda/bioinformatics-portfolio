"""Root pytest configuration and dynamic multi-project module loader."""

import sys
from pathlib import Path

import pytest


def _set_project_path(file_path: Path) -> None:
    for parent in file_path.parents:
        if parent.parent.name == "projects":
            proj_str = str(parent)
            if proj_str in sys.path:
                sys.path.remove(proj_str)
            sys.path.insert(0, proj_str)
            # Evict previous 'src' modules from cache
            sys.modules.pop("src", None)
            for k in list(sys.modules.keys()):
                if k.startswith("src."):
                    del sys.modules[k]
            break


def pytest_pycollect_makemodule(module_path, parent):
    _set_project_path(Path(module_path))
    return None


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    _set_project_path(Path(item.path))
