# INCA Tests

This test suite runs INCA with fixed input templates and compares outputs against reference files.

## Run

From repo root:

```bash
python3 tests/run_tests.py
```

`make test` does the same. The methane test (`tests/data/ch4.*`) ships with the
repository. The H3 and Li tests need large files (~260 MB) that are not in the
repository; download them from Zenodo
(<https://doi.org/10.5281/zenodo.23184299>) and put the files
(`H3_FCI_cc-pVTZ.{fchk,dm2p}`, `Li_FCI_cc-pVTZ.{fchk,dm2p}`) into `tests/data/`. If they are
absent those tests are reported as SKIPPED, not failed.

The harness uses paths relative to `tests/` and will create `tests/work/` on each run.

## Current Input Format

All active tests use NAMELIST input blocks:

- `&files` for input filenames
- `&intracule_job` for intracule options
- `&cube_job` for cube options

See `tests/inputs/` for templates.

## References

Golden references live in `tests/refs/`. Missing references are treated as failures. Update them only when you intentionally change numerical results.

## Notes

- Legacy inputs in `tests/inputs/old_inputs/` are not used by the current code.
