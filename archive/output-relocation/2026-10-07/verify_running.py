"""Observe the resumed worker without importing or changing scientific code."""

import hashlib
import io
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np


evidence = Path(__file__).resolve().parent
resume = json.loads((evidence / "resume.json").read_text())
manifest_bytes = (evidence / "manifest.json").read_bytes()
assert hashlib.sha256(manifest_bytes).hexdigest() == resume["manifest_sha256"]
manifest = json.loads(manifest_bytes)
previous = resume["resumed_checkpoint"]
checkpoint = Path(previous["path"])
checkpoint_bytes = checkpoint.read_bytes()
with np.load(io.BytesIO(checkpoint_bytes), allow_pickle=False) as data:
    metadata = json.loads(str(data["metadata"]))
    saved = json.loads(str(data["long_experiment"]))
    ids = data["ids"].tolist()

history_relative = str(
    checkpoint.with_name("history.json").relative_to(Path(manifest["source_root"]))
)
prefix_bytes = (
    json.dumps(saved["history"][: previous["history_rows"]], indent=2, allow_nan=False)
    + "\n"
).encode()
prefix_sha256 = hashlib.sha256(prefix_bytes).hexdigest()
assert prefix_sha256 == manifest["files"][history_relative]["sha256"]
assert saved["job"] == previous["job"]
assert saved["protocol_hash"] == resume["protocol_sha256"]
assert ids == previous["ids"]
assert metadata["time"] > previous["time"]
assert metadata["step_number"] > previous["step_number"]

pid = resume["pid"]
process = subprocess.check_output(
    ["ps", "-p", str(pid), "-o", "pid,etimes,state,args"], text=True
)
assert "embryo.phase_carry_redistribution_response" in process
assert resume["manifest_sha256"] in process
gpu_output = subprocess.check_output(
    [
        "nvidia-smi",
        "--query-compute-apps=pid,gpu_uuid,process_name,used_gpu_memory",
        "--format=csv,noheader",
    ],
    text=True,
)
gpu_rows = [
    line for line in gpu_output.splitlines()
    if line.split(",", 1)[0].strip() == str(pid)
]
assert len(gpu_rows) == 1 and resume["device_uuid"] in gpu_rows[0]
observed = datetime.now(timezone.utc)
record = {
    "passed": True,
    "observed_utc": observed.isoformat(),
    "observed_local": observed.astimezone(ZoneInfo("America/New_York")).isoformat(),
    "pid": pid,
    "process_alive": True,
    "process_observation": process,
    "qualified_gpu_observation": gpu_rows[0],
    "device_uuid": resume["device_uuid"],
    "manifest_sha256": resume["manifest_sha256"],
    "protocol_sha256": saved["protocol_hash"],
    "status": json.loads((checkpoint.parent.parent / "status.json").read_text()),
    "original_checkpoint": previous,
    "observed_checkpoint": {
        "path": str(checkpoint),
        "sha256": hashlib.sha256(checkpoint_bytes).hexdigest(),
        "time": metadata["time"],
        "step_number": metadata["step_number"],
        "elapsed": metadata["time"] - saved["job"]["start"],
        "history_rows": len(saved["history"]),
        "ids": ids,
    },
    "original_history_prefix_preserved": True,
    "original_history_prefix_rows": previous["history_rows"],
    "original_history_prefix_sha256": prefix_sha256,
    "checkpoint_advanced": True,
    "observer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope": "Execution and checkpoint continuity only; no scientific outcome promoted.",
}
(evidence / "running-verification.json").write_text(
    json.dumps(record, indent=2, allow_nan=False) + "\n"
)
print(json.dumps({
    "passed": record["passed"],
    "pid": pid,
    "old_time": previous["time"],
    "new_time": metadata["time"],
    "original_history_prefix_preserved": True,
    "gpu_uuid": resume["device_uuid"],
}))
