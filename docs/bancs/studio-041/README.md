# Studio 0.4.1 — English workbench qualification

Local Linux x86-64 qualification with Python 3.14.7, kernel 0.19.0,
PySide6 6.11.2, VTK 9.7.0, OCCT 8.0.1.0 and adapted build123d
0.11.1+vinkulum.occt8. The [test log](studio-tests.txt) records **49 passing
Studio tests** under a real Qt/X11/OpenGL context in Xvfb.

```sh
PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh
```

The added language compatibility test checks English joint and law choices
against their preserved serialised identifiers. Existing tests check editing,
undo/redo, worker failure, capture of simulation inputs, physical reference
cases, result comparison, CAD mass properties and number precision across
Inspector sizes and themes.

The [CAD recipe report](recipe.json) comes from four real worker operations:
box, cylinder, subtraction and fillet. The [README screenshot](../../assets/studio-cad.png)
was captured at a 1440 × 950 logical window size with Qt scale factor 2. A
1280 × 844 capture at scale factor 1 was also inspected locally. These checks
cover the exercised workflows; they do not qualify all platforms or establish
a general scientific error bound.

The earlier [0.4.0 record](../studio-cad-040/README.md) contains the OCCT 8
adaptation's upstream test subset. Its results remain historical evidence for
that version. Each standalone 0.4.1 build must additionally pass checks on its
own extracted executable; consult the report attached to that distribution.
