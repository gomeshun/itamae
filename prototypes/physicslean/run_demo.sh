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
lake build
python3 -m unittest discover -s tests -v
for test_file in tests/*.lean; do
  lake env lean "$test_file"
done
# Check the shipped concrete witnesses, including the per-redshift equalities.
for proof in results/eps-eh/transfer.lean results/eps-eh/growth.lean; do
  {
    echo 'import HzAudit'
    cat "$proof"
    echo
    echo 'run_cmd PhysicsLeanAudit.currentWitnesses'
  } > "$example/CheckedWitness.lean"
  lake env lean "$example/CheckedWitness.lean"
done
