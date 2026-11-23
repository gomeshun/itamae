"""Fetch the official mathlib artifacts needed by the named Lean modules."""
import json
import subprocess
import sys


def main() -> None:
    targets = sys.argv[1:] or ["HzExample"]
    graph = subprocess.check_output(
        ["lake", "-q", "--json", "query", *(f"+{name}:transImports" for name in targets)],
        text=True,
    )
    imports = {mod for row in graph.splitlines() for mod in json.loads(row)}
    modules = sorted(mod for mod in imports if mod.startswith("Mathlib."))
    if not modules:
        raise SystemExit("Lake returned no mathlib dependencies for the requested modules")
    subprocess.run(["lake", "exe", "cache", "get", *modules], check=True)


if __name__ == "__main__":
    main()
