# Verification and reproduction

## Confirm this publication package's original file bytes

From the extracted research root:

```bash
python3 07_Reproduce_And_Verify/verify_release.py
```

This verifies the exact original file hashes after renaming, checks the extracted raw paired records, sums the condition summaries, and checks source publication figure digests. It does **not** claim to rerun physical experiments or the native computational program.

## Reconstruct the original historical release layout, if required

```bash
python3 07_Reproduce_And_Verify/restore_original_layout.py /tmp/original-molecular-publication
cd /tmp/original-molecular-publication
python3 reproduce/VERIFY_PUBLICATION.py
```

The first step copies the originals back into their historical paths in a *new*, unused output directory, without altering this organized package. It gives the old relative-path-dependent scripts their original environment. Those historic scripts and manifests are not directly runnable from the reorganized folder hierarchy.

## Execute the original computational model separately

Original compiled binaries, C/C++ sources, Python tools and their original tests are preserved inside `08_Executable_Source_Archive/`. The full native software requires compatible compiler/runtime dependencies, described by its own historical documentation. It is separate from simple publication-data verification.
