"""Exercise actual restore/refusal behavior with temporary text and binary files."""
import json
import os
import random
import stat
import subprocess
import tempfile
from pathlib import Path

from manage_archive import RECORD, ROOT, inside_root, restore_file, sha256, verify_blob, zstd

checks = []
with tempfile.TemporaryDirectory(prefix=".compression-restore-check-", dir=ROOT / "outputs") as work:
    work = Path(work)
    for name, contents in (("history.json", b'{"example": [1.0, 2.0, 3.0]}\n' * 4096),
                           ("checkpoint.npz", random.Random(7).randbytes(256 * 1024))):
        original = work / name
        original.write_bytes(contents)
        original.chmod(0o640)
        os.utime(original, ns=(1_700_000_000_123456789, 1_700_000_001_987654321))
        info = original.stat()
        blob = work / (name + ".zst")
        with blob.open("wb") as stream:
            subprocess.run([zstd(), "-q", "-3", "-c", str(original)], stdout=stream, check=True)
        record = dict(path=str(original.relative_to(ROOT)), blob=str(blob.relative_to(ROOT)),
                      sha256=sha256(original), compressed_sha256=sha256(blob), bytes=info.st_size,
                      atime_ns=info.st_atime_ns, mtime_ns=info.st_mtime_ns, mode=info.st_mode)
        verify_blob(record)
        original.unlink()
        assert restore_file(record)
        assert original.read_bytes() == contents and sha256(original) == record["sha256"]
        assert original.stat().st_mtime_ns == info.st_mtime_ns
        assert stat.S_IMODE(original.stat().st_mode) == 0o640
        assert original.stat().st_nlink == 1
        assert not restore_file(record)  # Idempotent; matching originals are kept.
        original.write_bytes(b"subsequent edits")
        try:
            restore_file(record)
        except FileExistsError:
            pass
        else:
            raise AssertionError("Restore overwrote a changed file")
        assert original.read_bytes() == b"subsequent edits"
        original.unlink()
        blob.write_bytes(blob.read_bytes() + b"corruption")
        try:
            restore_file(record)
        except ValueError:
            pass
        else:
            raise AssertionError("Corrupt blob was accepted")
        assert not original.exists()
        checks.append(dict(format=name, exact_byte_restore=True, mode_and_mtime_preserved=True,
                           independent_restored_inode=True, idempotent=True,
                           changed_original_refused=True, corrupt_archive_refused=True))
for unsafe in ("/tmp/file", "outputs/../embryo/model.py", "embryo/model.py"):
    try:
        inside_root(unsafe)
    except ValueError:
        pass
    else:
        raise AssertionError("Unsafe restore path accepted")
result = dict(passed=True, checks=checks, unsafe_paths_refused=True)
(RECORD / "restore-validation.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
