"""Temporary helper to update packaging docs for the sashimi-itamae rename."""

from __future__ import annotations

import os
from pathlib import Path

PIN = os.environ["PINNED_RENAME_SHA"]


def replace_once(path: str, old: str, new: str) -> None:
    target = Path(path)
    text = target.read_text()
    if old not in text:
        raise RuntimeError(f"expected text not found in {path}: {old[:80]!r}")
    target.write_text(text.replace(old, new, 1))


replace_once(
    "notebooks/usage_walkthrough.ipynb",
    'version(\\"itamae\\")',
    'version(\\"sashimi-itamae\\")',
)

replace_once(
    "README.md",
    """The PyPI distribution rename is a release blocker tracked in
[#3](https://github.com/gomeshun/itamae/issues/3); the Python import namespace
remains `itamae` during development.""",
    """The public Python **distribution** is named `sashimi-itamae` to avoid the
unrelated PyPI project named `itamae`. The Python **import namespace** remains
`itamae`. This packaging decision is tracked in
[#3](https://github.com/gomeshun/itamae/issues/3).""",
)

replace_once(
    "README.md",
    """The SASHIMI ITAMAE core is not released on PyPI. That distribution name is
occupied by an unrelated project, so do not install it by an unqualified index
name. Use the reviewed source revision for reproducible development:

```bash
uv pip install \"itamae[full] @ git+https://github.com/gomeshun/itamae.git@1c5b1ad67725671fd4dc1a2306d2731ba337e563\"
```""",
    f"""The SASHIMI ITAMAE core is not yet released on PyPI. Its collision-free
distribution name is `sashimi-itamae`, while code continues to use
`import itamae`. Until publication, install an exact reviewed VCS revision:

```bash
uv pip install \"sashimi-itamae[full] @ git+https://github.com/gomeshun/itamae.git@{PIN}\"
python -c \"import itamae; print(itamae.__version__)\"
```

Do **not** install bare `itamae` from PyPI; that name belongs to an unrelated
project. Downstream SASHIMI package metadata should depend on `sashimi-itamae`,
but Python imports do not change.""",
)

replace_once(
    "PLAN.md",
    """ITAMAE is intended for publication under a collision-free distribution name
selected in [GOV-02 / #3](https://github.com/gomeshun/itamae/issues/3). The current
metadata still says `itamae`, but that PyPI name belongs to an unrelated
project. The Python import namespace can remain `itamae`. No new distribution
identity is selected by this document.""",
    """ITAMAE uses the collision-free public distribution name `sashimi-itamae`,
selected in [GOV-02 / #3](https://github.com/gomeshun/itamae/issues/3). The
Python import namespace remains `itamae`. The unrelated PyPI distribution named
`itamae` must never be used as the SASHIMI core dependency.""",
)

replace_once(
    "PLAN.md",
    'name = "itamae"',
    'name = "sashimi-itamae"',
)

replace_once(
    "PLAN.md",
    """Install the reviewed development source explicitly until the naming issue is
resolved. Do not use `uv add itamae` to obtain the SASHIMI core from PyPI.

```bash
uv pip install \"itamae[full] @ git+https://github.com/gomeshun/itamae.git@1c5b1ad67725671fd4dc1a2306d2731ba337e563\"
```""",
    f"""Install the reviewed development source explicitly until `sashimi-itamae` is
published. Do not use `uv add itamae`; the import name and distribution name are
deliberately different. Downstream requirements use `sashimi-itamae`, while
Python code continues to import `itamae`.

```bash
uv pip install \"sashimi-itamae[full] @ git+https://github.com/gomeshun/itamae.git@{PIN}\"
```""",
)
