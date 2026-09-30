"""Modules run with ``python -m`` execute their ``if __name__ == "__main__"`` block
the moment it is reached, so anything defined below it does not exist yet.
multisample_sv_calling once defined save_svComposites_to_json below its main
block, and --tmp-dir-path died with a NameError when run as a module."""

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "svirlpool"


def _is_main_guard(node: ast.stmt) -> bool:
    return (
        isinstance(node, ast.If)
        and isinstance(node.test, ast.Compare)
        and isinstance(node.test.left, ast.Name)
        and node.test.left.id == "__name__"
    )


@pytest.mark.parametrize(
    "module", sorted(SRC.rglob("*.py")), ids=lambda p: str(p.relative_to(SRC))
)
def test_the_main_block_is_the_last_statement(module: Path):
    body = ast.parse(module.read_text()).body
    guards = [i for i, node in enumerate(body) if _is_main_guard(node)]
    if guards:
        after = [getattr(n, "name", type(n).__name__) for n in body[guards[-1] + 1 :]]
        assert not after, f"defined after the main block: {after}"
