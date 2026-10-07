"""Read-only archive selection: completed studies outside live dependency roots."""
from __future__ import annotations

import collections
import hashlib
import json
import re
from pathlib import Path

RECORD = Path(__file__).resolve().parent
ROOT = RECORD.parents[2]
ACTIVE = ROOT / "outputs/phase-carry-delayed-response"
PATTERN = re.compile(rb"outputs/[A-Za-z0-9_.+*/\-]+")
KEEP_READABLE = {"protocol.json", "status.json", "summary.json", "result.json",
                 "analysis.json", "assessment.json", "preflight.json"}


def main():
    queue = [ACTIVE / "protocol.json"]
    seen, references = set(), set()
    while queue:
        path = queue.pop()
        if path in seen:
            continue
        seen.add(path)
        with path.open("rb") as stream:
            tail = b""
            while block := stream.read(1024 * 1024):
                joined = tail + block
                for match in PATTERN.finditer(joined):
                    relative = match.group(0).decode().rstrip(".")
                    references.add(relative)
                    target = ROOT / relative
                    if target.suffix == ".json" and target.is_file() and target not in seen:
                        queue.append(target)
                tail = joined[-4096:]
    protected = sorted({"/".join(p.split("/")[:2]) for p in references}
                       | {str(ACTIVE.relative_to(ROOT))})
    completed, texts, groups = [], [], collections.defaultdict(list)
    for folder in (ROOT / "outputs").iterdir():
        relative = str(folder.relative_to(ROOT))
        status = folder / "status.json"
        if not folder.is_dir() or relative in protected or not status.is_file():
            continue
        if json.loads(status.read_text()).get("state") != "completed":
            continue
        completed.append(relative)
        for path in folder.rglob("*"):
            if path.is_symlink() or not path.is_file():
                continue
            info = path.stat()
            if info.st_size <= 1024 * 1024 or info.st_nlink != 1:
                continue
            if path.suffix in (".json", ".html") and path.name not in KEEP_READABLE:
                texts.append(str(path.relative_to(ROOT)))
            elif path.suffix == ".npz":
                groups[info.st_size].append(path)
    duplicates = []
    for size, paths in groups.items():
        if len(paths) < 2:
            continue
        hashes = collections.defaultdict(list)
        for path in paths:
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            hashes[digest].append(str(path.relative_to(ROOT)))
        for digest, names in hashes.items():
            if len(names) > 1:
                duplicates.append(dict(sha256=digest, bytes_each=size, paths=sorted(names)))
    report = dict(active=str(ACTIVE.relative_to(ROOT)),
                  protocol_sha256=hashlib.sha256((ACTIVE / "protocol.json").read_bytes()).hexdigest(),
                  closure_json=sorted(str(p.relative_to(ROOT)) for p in seen),
                  references=sorted(references), protected_roots=protected,
                  completed_unprotected_roots=sorted(completed),
                  text_candidates=sorted(texts), duplicate_npz_groups=duplicates)
    destination = RECORD / "audit.json"
    if destination.exists():
        raise FileExistsError(destination)
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(report=str(destination.relative_to(ROOT)),
                          protected_roots=len(protected), json_files_scanned=len(seen),
                          text_candidates=len(texts), duplicate_npz_groups=len(duplicates)), indent=2))


if __name__ == "__main__":
    main()
