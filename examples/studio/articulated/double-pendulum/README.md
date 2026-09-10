# Inspect two-link operators without installing Pinocchio

This is a captured calculation from **Studio 0.6.0.dev1 / Pinocchio 4.1.0**,
generated on Linux using `ci/studio_pinocchio_recipe.py` and the separate
installed worker. It supplies the actual matrix values shown in the
[workspace screenshot](../../../../docs/assets/studio-pinocchio.png).

Start Studio **0.6.0a1** (Linux preview or current source), open **Run → Articulated
operators · Pinocchio… → Open result…**, and select this folder's `result.json`.
No installed Pinocchio engine is needed to inspect it. The
[Linux alpha](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.6.0a1-linux)
ships these files in `Examples/pinocchio`. Earlier 0.5.0 binaries predate this workspace.

The folder contains the captured Studio project, requested state and result.
Keep the three JSON files together and unchanged: reopening checks their
fingerprints and the operator contract. The two revolute coordinates are
`q = [0.2, -0.3] rad`, velocity `[0.1, 0.4] rad/s`, requested acceleration
`[0.7, 0.2] rad/s²`, and actuator effort `[0.9, 0.1] N m`. Requested acceleration
and actuator effort are independent inputs for inverse and forward dynamics;
the displayed state is not an integrated trajectory.

To recompute with the qualified worker, from the repository root:

```sh
.venv-pinocchio/bin/python -m vinkulum_studio.pinocchio_backend \
  examples/studio/articulated/double-pendulum/project.vinkulum.json \
  --state examples/studio/articulated/double-pendulum/requested-state.json \
  --output /tmp/vinkulum-two-link-new-run
```

Use a new output directory. The result records actual library and adapter
fingerprints; a recomputation after source changes can have different provenance.
The scientific status is `NotAssessed`. Read the
[operator meanings, reference equations and limits](../../../../docs/PINOCCHIO_OPERATORS.md).
Original example data are distributed under Vinkulum's Apache-2.0 licence.
