# Studio background admission — Linux qualification

The [report](qualification.json) records the final installed wheel hash, all 61
Studio module hashes, versions and results. The two fresh environments matched
the wheel and source modules byte for byte. They used no editable Studio install
or `PYTHONPATH` overlay. Kernel version 0.20.0 is unchanged.

- [Studio](studio-tests.log.gz): 188 tests, 17 expected skips, no failures.
- [Separate Pinocchio worker](pinocchio-tests.log.gz): all 11 original tests and
  six applied-load derivative tests pass. These cover the 16 numerical cases
  skipped in the CAD/GUI environment; tree admission runs in both environments.
- [Keyboard CAD](keyboard.log.gz): the remaining opt-in GUI recipe passes at
  100%/200% scaling with both 420/280 px inspector widths.
- [Repeated Qt lifecycle checks](qt-stress.log.gz): the five background-admission
  tests repeated ten times in one process, 50 tests without failure. The window
  and CAD-dialog tests each exercise both implementations.
- Normal application startup, native dynamics, the OCCT worker self-check and
  a nonblack rendered scene pass. The application was installed Python code,
  not a newly packaged standalone binary.

Reproduce the full suite with `ci/studio.sh`, setting `PY` to the installed
CAD/GUI environment and `VINKULUM_PINOCCHIO_PYTHON` to the separate installed
worker. Run `ci/studio_keyboard.sh` with the same `PY` for dedicated X11/Openbox
keyboard sessions. The targeted cases are in
[`test_background_admission.py`](../../../apps/studio/tests/test_background_admission.py);
load its `BackgroundAdmission` class into one `unittest.TestSuite` ten times to
repeat the stress run, with `VINKULUM_3D_TESTS=1` under Xvfb.

Controlled barriers test event-loop progress, lease ownership and cancellation
publication boundaries. They do not measure worst-case JSON latency or certify
arbitrary concurrent operations. The
[background-admission contract](../../STUDIO_BACKGROUND_ADMISSION.md) lists
the remaining GUI work and the Python GIL limitation.
