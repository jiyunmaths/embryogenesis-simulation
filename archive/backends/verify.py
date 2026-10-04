"""Verify retired prototype copies without importing or launching a backend."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent.parent
manifest = json.loads((root/'manifest.json').read_text())
for entry in manifest['files']:
    for filename in (repo/entry['original'], root/entry['snapshot']):
        if hashlib.sha256(filename.read_bytes()).hexdigest() != entry['sha256']:
            raise SystemExit('Archive/original source changed: '+str(filename))
print(f"Verified {len(manifest['files'])} archival source snapshots and retained originals.")
