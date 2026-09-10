# Vinkulum dans FreeCAD

FreeCAD est désormais l’interface de référence pour les nouveaux développements.
L’atelier Vinkulum utilise ses outils PartDesign, son arbre de document et son
panneau de tâches. Le moteur mécanique reste dans un processus séparé.

Cette première version couvre **une pièce rigide et un pivot fixé au monde**.
Le panneau définit la masse volumique, la position du pivot en millimètres, son
axe global X/Y/Z, la durée et le pas. Le calcul conserve la capture géométrique,
les propriétés physiques en SI et les échantillons natifs du mouvement. La lecture
anime une copie dédiée ; la pièce paramétrique reste éditable dans FreeCAD.

## Installer l’atelier

Construire l’archive depuis ce dépôt avec un Python standard :

```sh
python apps/freecad/package.py /tmp/Vinkulum-FreeCAD-0.1.0a1.zip
```

Dans la console Python de FreeCAD, `App.getUserAppDataDir()` donne le dossier
utilisateur de cette installation. Extraire le dossier `Vinkulum` de l’archive
dans son sous-dossier `Mod`, puis redémarrer FreeCAD. Choisir **Vinkulum** dans
la liste des ateliers. L’archive contient le code de l’atelier, son icône, sa
licence et les empreintes des fichiers.

Le premier paquet nécessite un environnement moteur séparé, préparé avec les
[dépendances CAO de Vinkulum](../../docs/STUDIO_CAD.md). Dans
**Vinkulum → Configurer le moteur…**, sélectionner son exécutable Python et le
dossier de conservation des calculs. Aucun lancement de l’interface Studio
n’est nécessaire. Pour une installation automatisée,
`VINKULUM_FREECAD_PYTHON` fournit la valeur initiale du chemin du moteur.

## Premier calcul

1. **Vinkulum → Exemple : pendule paramétrique** crée une esquisse et une
   extrusion PartDesign de 800 mm dans un nouveau document FreeCAD.
2. Sélectionner le corps et choisir **Analyser la pièce sélectionnée…**.
   Les valeurs initiales correspondent à une pièce de 20 × 30 mm, de masse
   volumique 7800 kg/m³, avec un pivot à l’origine sur l’axe global Y.
3. **Calculer**, puis **Lire le mouvement**. Le curseur sélectionne directement
   les échantillons calculés ; les commandes de vue restent celles de FreeCAD.
4. Fermer la tâche pour reprendre la modélisation. Modifier la longueur dans
   FreeCAD et relancer une capture pour calculer cette nouvelle pièce.

Chaque résultat devient une entrée dans l’arbre FreeCAD. Enregistrer le document
avec la commande habituelle, puis sélectionner cette entrée après réouverture
pour retrouver son mouvement sans relancer le moteur. La géométrie capturée est
conservée dans le `.FCStd`, et reconstruite à partir de cette capture à la reprise.
Les trajectoires et fichiers du calcul restent dans le dossier indiqué par la
propriété `VinkulumDirectory` du résultat : ce dossier doit rester accessible à
son emplacement. Le `.FCStd` seul ne constitue pas encore un paquet de transfert
complet vers une autre machine.

Un seul calcul Vinkulum est admis à la fois dans FreeCAD. Annuler, supprimer la
pièce source ou fermer son document arrête le processus avant de notifier la fin.
Un résultat antérieur reste disponible après annulation. Si la géométrie a changé
pendant le calcul, son résultat est refusé et une nouvelle capture est nécessaire.

## Qualification et suite

La [qualification conservée de l’atelier](../../docs/bancs/freecad-workbench-2026/README.md)
contient l’archive installable, trois calculs réels, les documents FreeCAD,
les références indépendantes et les dix scénarios d’arrêt du moteur.

Les recettes [qualify_workbench.FCMacro](qualify_workbench.FCMacro) et
[qualify_lifecycle.FCMacro](qualify_lifecycle.FCMacro) s’exécutent dans une session
FreeCAD dédiée sur l’archive installée. Elles couvrent le calcul et la lecture,
l’enregistrement/réouverture et la durée de vie des processus. La
[référence physique initiale](../../docs/bancs/freecad-bridge-2026/README.md)
compare les propriétés et le mouvement à un pendule physique indépendant.

Le premier environnement ciblé est FreeCAD 1.1.3 Linux x86-64, Python 3.11 et
PySide6. La conversion des assemblages FreeCAD, les liaisons attachées aux faces,
les corps multiples, les analyses FEM dans FreeCAD et la livraison macOS/Windows
restent à qualifier. La [direction d’intégration](../../docs/FREECAD_INTEGRATION.md)
décrit ces étapes.

## Contrat et expérience de transport initiale

This prototype tests using FreeCAD's existing parametric modelling interface
with Vinkulum's mechanics engine in a separate process. It edits a PartDesign
pad, captures its solid and physical properties, runs an explicit revolute
mechanism, and displays the returned native poses on a separate copy in FreeCAD.
The parametric source is unchanged by playback.

The [retained qualification](../../docs/bancs/freecad-bridge-2026/README.md)
contains four real desktop runs, editable `.FCStd` sources, STEP captures,
reopenable Vinkulum projects, native trajectory archives, screenshots and an
independent physical-pendulum reference. That initial feasibility experiment
established the transport now used by the installable workbench above.

## Process and geometry contract

```mermaid
flowchart LR
    F[FreeCAD parametric source] --> C[Captured STEP and SI properties]
    C --> V[Separate Vinkulum Python process]
    V --> R[Validated native trajectory]
    R --> P[Captured-shape playback in FreeCAD]
```

The tested official **FreeCAD 1.1.3** AppImage uses Python 3.11.14 and OCCT 7.8.1.
The Vinkulum worker uses Python 3.14.7, Studio 0.6.0a2.dev6, kernel 0.20.0 and
OCCT 8.0.1.0. Each process keeps its own Python and native libraries. The bridge
removes inherited `PYTHONHOME` and library/plugin overrides before launching the
explicitly selected Vinkulum interpreter; the AppImage sets those for FreeCAD's
own runtime. FreeCAD already supports external Python extensions through its
[workbench interface](https://github.com/FreeCAD/FreeCAD-documentation/blob/main/wiki/Workbench_creation.md).

[bridge.py](bridge.py) runs inside FreeCAD and captures one valid solid in a
fresh directory. STEP bytes, a request UUID, body identity, local BREP and global
placement identify that capture. Volume, centre and full centroidal inertia are
also recorded from FreeCAD. Its millimetre quantities become SI using `10⁻⁹`
for volume and `ρ × 10⁻¹⁵` for the inertia integral in mm⁵.

[worker.py](worker.py) runs in the existing Studio CAD environment. The production
OCCT 8 importer reads the STEP and recomputes all physical properties. Volume,
mass, centre and every inertia component must agree with the FreeCAD capture
within scale-dependent absolute budgets (`εV`, `εm`, `εL`, `εmL²`, `ε = 10⁻⁸`;
L is the imported bounding-box diagonal). The normal mechanical adapter creates
the explicitly captured world-frame revolute joint and validates the native
trajectory archive. The initial experiment used the world origin and Y axis;
the workbench exposes the pivot position and global X/Y/Z axes.
The selected allocation is two threads; CAD import and
mechanics run sequentially in this worker.

FreeCAD receives JSON poses in metres and row-major body-to-world matrices.
The bridge checks capture identity, increasing times, finite values, proper
rotations and the initial pose. It also rechecks the source geometry and its
global placement. Moving an enclosing Part can leave the local BREP unchanged;
that inherited placement is therefore part of the signature. Playback applies
the rigid transform about the captured centre to a separate `Part::Feature`.
It does not write the motion into the original pad or its design placement.
Replacing the input file also invalidates admission, even if a result carries
the hash of that replacement: identity is checked against the original capture.

The QProcess has a 60-second timeout and bounded diagnostic capture. The workbench
admits one Vinkulum job at a time per FreeCAD application and waits for the child
to stop before reporting completion, cancellation or failure. Its qualification
covers cancellation during startup, execution and after result preparation,
source deletion, document/application closure, timeout and stale geometry.
FreeCAD's ordinary save/cancel decision remains available when quitting.
Scheduling a broader set of concurrent solver services remains future work.

## Reproduce

Install the [Studio CAD environment](../../docs/STUDIO_CAD.md), then use a
dedicated FreeCAD process. The qualification creates and closes its own document
and exits that process; it refuses to start if documents are already open.
Use a new output directory:

```sh
export VINKULUM_FREECAD_PROTOTYPE="$PWD/apps/freecad"
export VINKULUM_FREECAD_PYTHON=/absolute/path/to/studio-environment/bin/python
export VINKULUM_FREECAD_CAPTURE=/tmp/freecad-vinkulum-capture
FreeCAD "$VINKULUM_FREECAD_PROTOTYPE/qualify.FCMacro"
"$VINKULUM_FREECAD_PYTHON" apps/freecad/verify.py \
  "$VINKULUM_FREECAD_CAPTURE" --plot
```

`FreeCAD` above denotes the executable for a normal FreeCAD installation. On
headless Linux, prefix it with `xvfb-run -a -s '-screen 0 1600x1100x24'`.
The qualification used the extracted official AppImage's `AppRun`, with
`XDG_CONFIG_HOME`, `XDG_DATA_HOME` and `XDG_CACHE_HOME` pointing to separate
temporary directories. The release URL and checked checksum are in the record.

[qualify.FCMacro](qualify.FCMacro) changes the pad length from 800 to 1200 mm,
with native time steps of 5 and 2.5 ms for each length. It checks all displayed
centres against the returned poses, confirms that the source stays unchanged,
and rejects changed geometry, a moved parent, a corrupted rotation and a
replaced request. A Qt
timer records continued event-loop delivery during each external calculation.

[verify.py](verify.py) needs numpy, and matplotlib only for `--plot`. It imports
neither FreeCAD, Studio nor a geometry/physics engine. It checks input/result
hashes, exact agreement between the native archive and the replayed data, the
analytic box mass/inertia, pivot kinematics and convergence against an independent
scalar pendulum equation. It can also verify the extracted retained archive.

The normal Studio suite includes [two external-worker regressions](../studio/tests/test_freecad_bridge.py)
using the retained FreeCAD STEP. They require the Studio CAD environment but no
running FreeCAD: one checks an actual native trajectory against the independent
pendulum equation; the other changes the captured mass and checks rejection
before simulation. Run them directly with:

```sh
"$VINKULUM_FREECAD_PYTHON" -m unittest discover \
  -s apps/studio/tests -p test_freecad_bridge.py -v
```

## Scope of the result

The positive experiment is one uniform rectangular rigid body, one explicitly
defined revolute joint, two lengths and two seconds of motion from rest at 20°.
It establishes this parametric-edit → physical-capture → mechanics → native-pose
display path across the two runtimes. It does not establish general assembly
constraint conversion, persistent face attachment, arbitrary mechanisms,
multi-window scheduling, FEM or macOS/Windows packaging. STEP transfers evaluated
geometry; the original FreeCAD feature graph remains in `.FCStd`.

All bridge code, models, reference calculations and capture data are original
Vinkulum contributions under [Apache-2.0](../../LICENSE). The unchanged official
FreeCAD runtime remains available from its upstream project with its own licence.
