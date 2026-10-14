# Output storage on the data drive

Large simulation artifacts now use these data-drive folders:

```text
/mnt/data/Documents/Python/embryogenesis-simulation/outputs
/mnt/data/Documents/Python/embryogenesis-simulation/output/studies
/mnt/data/Documents/Python/embryogenesis-simulation/output/kernel-cache
/mnt/data/Documents/Python/embryogenesis-simulation/output/pdf
```

The corresponding repository paths are symbolic links. Keep `/mnt/data` mounted when using the simulations. Code, Git history and manuscript remain in the original repository. The small `output/figures` and `output/presentations` folders remain physical local directories to preserve their six Git-tracked artifacts; verified copies also exist under the data-drive `output` folder.

The move preserves file contents, permissions, modification times and the fourteen existing output symlinks. All 8,253 regular files were compared with SHA-256 before switching locations. The [migration manifest](../archive/output-relocation/2026-10-07/manifest.json) records original and destination metadata. It describes the cutover snapshot; running jobs can subsequently update their own checkpoints and logs.

## New study storage move: 10 October 2026

The singular `output` tree had accumulated another **37.38 GiB of file contents**. All **14,960 regular files**, three symlinks and 5,500 directories passed source/destination SHA-256 and stable metadata checks before cutover. The verified local backup was removed only after all 33 active-study source pins, 38 input pins and 1,774 preserved evidence pins passed, including exact reuse provenance. Restoring the small tracked folders leaves approximately **37.39 GiB net local space recovered**. [Operation report](../archive/output-relocation/2026-10-10-output/report.md), [full cutover manifest](../archive/output-relocation/2026-10-10-output/manifest.json), [final layout](../archive/output-relocation/2026-10-10-output/layout.json).

The matched mechanical-control coordinator was gracefully interrupted for a stable copy and resumed on the same designated GPU as PID **470605**. Four completed new trajectories and the accepted feedback-off baseline are verified and reused. Only the unfinished `polarity_space88` path restarts from its unchanged initial source. Its complete interrupted attempt is retained under `output/studies/development-control-qualification/storage-interrupted-attempts/2026-10-10/polarity_space88`. An independent saved-data check already reproduces all 54 pre-interruption observations and all three original checkpoint arrays exactly; the resume worker also checks that full prefix when the repeat finishes. Model, source/kernel bytes, protocol and numerical limits remain unchanged; storage movement does not establish a scientific result.

Some protocols recorded before the move expect the original canonical study filenames. Their ordinary filesystem access works through the link, but commands that validate canonical names require `scripts/with_study_output_path_alias.py`, with `--record archive/output-relocation/2026-10-10-output/studies-manifest.json` and `--record-sha256 8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6` before the usual `-m MODULE ACTION ...` arguments. This launcher uses the existing tested path-alias implementation for the verified study subtree only. The older plural `outputs` launcher and its manifest remain unchanged. [Study subtree manifest](../archive/output-relocation/2026-10-10-output/studies-manifest.json), [resumption record](../archive/output-relocation/2026-10-10-output/resume.json), [path-alias checks](../archive/output-relocation/2026-10-10-output/path-alias-final-tests.log).

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
