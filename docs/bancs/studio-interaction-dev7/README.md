# Initial local interaction diagnostic

These two Linux/X11 observations use identical actions and the same dependency
versions, with an installed dev6 wheel before and editable dev7 work after.
They measure command synchronous return and the next event-loop callback.
They are exploratory observations, not a benchmark campaign or Mac evidence.
The final installed-wheel and native packaging checks must accompany a release.

The first candidate measured here predates the visible Back action and Workspaces
button. No claim that these samples describe the final source hash is made.
See [interaction scope](../../STUDIO_INTERACTION.md).

The subsequent installed-wheel observations are `installed-before.json` and
`installed-after.json`. They retain platform, package versions, all application
module SHA-256 values, OpenGL capabilities and timing boundaries. The recipe is
[`ci/studio_interaction_recipe.py`](../../../ci/studio_interaction_recipe.py).
The plate command blocked synchronously for 286.2 ms before and
99.2 ms after in that pair. Window size is 1280×800, on a 1600×1100
Xvfb screen; viewport dimensions are recorded in each file. The screenshot is the actual installed dev7
application. No first-frame deadline or Mac speedup follows from these samples.

The retained pairs used the default Linux Qt application style. The current
recipe explicitly selects Fusion and en_GB, as the ordinary Studio entry point
does; native workflow observations use that final recipe.
