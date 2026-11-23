#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
example=$(mktemp -d)
trap 'rm -rf "$example"' EXIT
python3 scripts/hz_cert.py --h0 70 --omega-m 3/10 --z 1 --digits 8 \
  --certificate "$example/hz_z1.json" --lean "$example/HzExample.lean"
python3 scripts/hz_cert.py --verify "$example/hz_z1.json"
diff -u generated/hz_z1.json "$example/hz_z1.json"
diff -u generated/HzExample.lean "$example/HzExample.lean"
python3 -m unittest discover -s tests -v
lake build
