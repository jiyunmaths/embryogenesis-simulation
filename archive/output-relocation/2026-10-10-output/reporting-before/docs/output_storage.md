# Output storage on the data drive

Simulation outputs now live at:

```text
/mnt/data/Documents/Python/embryogenesis-simulation/outputs
```

The repository's original `outputs` directory is a symbolic link to that location. Keep `/mnt/data` mounted when using the simulations. Code, Git history, manuscript and the separate `output/presentations` folder remain in the original repository.

The move preserves file contents, permissions, modification times and the fourteen existing output symlinks. All 8,253 regular files were compared with SHA-256 before switching locations. The [migration manifest](../archive/output-relocation/2026-10-07/manifest.json) records original and destination metadata. It describes the cutover snapshot; running jobs can subsequently update their own checkpoints and logs.

## Continue or assess an existing recorded protocol

Some historical checks compare canonical absolute evidence names. A normal filesystem symlink changes those names when `Path.resolve()` is called. The external launcher preserves their original logical names while file access follows the new drive. It changes no scientific source, equations, parameters or recorded decision.

For example, **after the running redistribution-response study finishes**, assess it with:

```bash
python scripts/with_output_path_alias.py \
  --record archive/output-relocation/2026-10-07/manifest.json \
  --record-sha256 edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02 \
  -m embryo.phase_carry_redistribution_response assess
```

Use the same launcher prefix for other module commands that continue or assess a protocol recorded before relocation. Ordinary file access, viewers and the dashboard retain their original output paths through the link.

## Restore compressed historical files

Use the relocation-aware archive wrapper rather than calling the frozen managers directly:

```bash
python scripts/relocated_output_archive.py 2026-10-06 list
python scripts/relocated_output_archive.py 2026-10-06 restore outputs/domain-conservative-extended
python scripts/relocated_output_archive.py 2026-10-07-redistribution-launch verify
```

The wrapper permits only the approved storage link and checks relocated file hashes and metadata before translating device/inode identities. Original manifests and managers remain unchanged; restored copies remain independent files. See [historical compression instructions](output_compression.md).

## Operation evidence

- [Copy verification](../archive/output-relocation/2026-10-07/verification.log) and [cutover](../archive/output-relocation/2026-10-07/cutover.json).
- [Scientific compatibility](../archive/output-relocation/2026-10-07/scientific-compatibility.json) and [archive compatibility](../archive/output-relocation/2026-10-07/archive-compatibility.json).
- [Original-copy retirement](../archive/output-relocation/2026-10-07/retirement.json), [worker resume](../archive/output-relocation/2026-10-07/resume.json), [checkpoint advancement](../archive/output-relocation/2026-10-07/running-verification.json), and [operation report](../archive/output-relocation/2026-10-07/report.md).

The previous worker was already absent when inspected. The relaunched worker advanced beyond its retained checkpoint on the same qualified GPU, with its original history preserved. The move does not establish any new scientific outcome.
