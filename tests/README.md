# INCA Tests

This test suite runs INCA with fixed input templates and compares outputs against reference files.

## Run

From repo root:

```bash
python3 tests/run_tests.py
```

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
- Scripts `test.sh`, `test2.sh`, and `tests/parallel_test/*` still reference `roda.exe` and the legacy input format. They are not part of the current automated test path.
