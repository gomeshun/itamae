"""Replay exact scalar or spectrum endpoints with the Lean kernel."""

from fractions import Fraction
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile


def literal(x: Fraction) -> str:
    return f"({x.numerator}/{x.denominator} : ℚ)"


def _budget(name: str, default: int) -> int:
    text = os.environ.get(name, str(default))
    if not re.fullmatch(r"[0-9]+", text) or int(text) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(text)


def check_interval(
    source: str, lo: Fraction, hi: Fraction, project: Path, *, max_memory_mb: int | None = None
) -> None:
    """Prove equality for a closed Lean definition named result, with a bounded checker."""
    timeout = _budget("PHYSICSLEAN_REPLAY_TIMEOUT", 120)
    heartbeats = _budget("PHYSICSLEAN_REPLAY_HEARTBEATS", 8_000_000)
    if max_memory_mb is None:
        max_memory_mb = _budget("PHYSICSLEAN_REPLAY_MEMORY_MB", 4096)
    proof = f"""
set_option maxRecDepth 100000 in
set_option maxHeartbeats {heartbeats} in
theorem concrete : result = some ⟨{literal(lo)}, {literal(hi)}⟩ := by decide +kernel
#print axioms concrete
"""
    check_proof(source + proof, project, max_memory_mb=max_memory_mb, timeout_seconds=timeout)


def check_proof(source, project, *, max_memory_mb, timeout_seconds=120):
    """Check an explicit proof file and audit its final theorem named concrete."""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "Replay.lean"
        path.write_text(source, encoding="utf-8")
        command = ["lake", "env", "lean", f"-M{max_memory_mb}", "-j1", str(path)]
        # Lake spawns Lean as a child. Stop both on timeout so a costly proof
        # cannot continue consuming memory after its caller has returned.
        with subprocess.Popen(
            command,
            cwd=project,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            start_new_session=os.name == "posix",
        ) as check:
            try:
                output, _ = check.communicate(timeout=timeout_seconds)
            except subprocess.TimeoutExpired as exc:
                if os.name == "posix":
                    os.killpg(check.pid, signal.SIGKILL)
                elif os.name == "nt":
                    subprocess.run(
                        ["taskkill", "/F", "/T", "/PID", str(check.pid)],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                else:
                    check.kill()
                check.communicate()
                raise RuntimeError(
                    f"Concrete kernel replay exceeded {timeout_seconds} seconds"
                ) from exc
    if check.returncode:
        raise RuntimeError(f"Concrete kernel replay failed\n{output[-3000:]}")
    lines = output.splitlines()
    if "'concrete' does not depend on any axioms" in lines:
        return
    allowed = {"propext", "Classical.choice", "Quot.sound"}
    for line in lines:
        audit = re.fullmatch(r"'concrete' depends on axioms: \[(.*?)\]", line)
        if audit and set(audit[1].split(", ")) <= allowed:
            return
    raise RuntimeError(f"Unexpected concrete proof audit\n{output[-1500:]}")


def check_table(args, raw_rows, layout, h, lo, hi, project) -> None:
    """Replay the request using rationals in the table's original column order and units."""
    tokens = [str(arg) for arg in args]
    if layout not in {"wmap7", "camb", "canonical"} or any(
        not re.fullmatch(r"[a-zA-Z0-9./,\-]+", token) for token in tokens
    ):
        raise ValueError("invalid spectrum replay request")
    # Separate declarations avoid a very deeply nested list literal when
    # replaying an entire CAMB table rather than a few measured knots.
    definitions = []
    groups = []
    for start in range(0, len(raw_rows), 32):
        name = f"fileRows{start}"
        groups.append(name)
        knots = ", ".join(f"⟨{literal(x)}, {literal(y)}⟩" for x, y in raw_rows[start : start + 32])
        definitions.append(f"private def {name} : List SashimiTable.Knot := [{knots}]")
    rows = "\n".join(definitions)
    arguments = ", ".join(map(json.dumps, tokens))
    source = f"""import SashimiTableRequest
open FloatLib.Numerics
set_option maxRecDepth 100000
{rows}
private def fileRows : List SashimiTable.Knot := [{", ".join(groups)}].flatten
def result : Option RationalInterval := do
  let canonical ← SashimiTableRequest.canonicalRows "{layout}" {literal(h)} fileRows
  let enclosure ← (SashimiTableRequest.run [{arguments}] canonical).toOption
  enclosure
"""
    check_interval(source, lo, hi, project)
