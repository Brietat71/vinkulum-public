# Studio CAD 0.4.0 — qualification locale

Environnement : Linux x86-64, glibc 2.39, CPython 3.14.7, Qt/PySide6 6.11.2,
VTK 9.7.0, noyau Vinkulum 0.19.0. CAD : OCCT 8.0.1.0 via OCP,
build123d 0.11.1+vinkulum.occt8, ocpsvg 0.6.0+vinkulum.occt8.

- `studio-tests.txt` : 48 tests Studio, CAD, document, processus et rendu, tous
  réussis sous Xvfb. Le contrôle de focus utilise une fenêtre réellement active.
- `build123d-tests.txt` : 96 tests amont ciblés, tous réussis ; ce sous-ensemble
  ne qualifie pas toute l'API build123d.
- `recipe.json` : boîte, cylindre, soustraction et congé exécutés par de vrais
  workers CAD. Masse, volume, nombre de triangles et temps mesurés de chaque
  opération sont conservés.
- [Capture réelle](../../assets/studio-cad.png) et
  [document de la pièce](../../../examples/studio/platine-percee.vinkulum.json).

Les captures ont été contrôlées à 1280×844 en échelle 1 et à 1440×950 en échelle
2. Les dimensions de la fenêtre sont en pixels logiques. L'inspecteur adapte
la précision d'affichage à la largeur disponible, préserve les valeurs exactes
non modifiées et les restitue au focus et dans l'infobulle. Le test modifie
aussi une valeur de très faible amplitude puis vérifie son application exacte.

Les mesures de la recette sont des observations locales, sans protocole de
benchmark statistique : l'import à froid prend environ 1,7 s par worker ;
l'opération et sa capture prennent environ 8 ms pour la boîte et 153 ms pour le
congé de cette pièce. Cela désigne le démarrage des workers comme une piste
d'optimisation ; cela ne prouve pas un avantage face à un autre logiciel.

Reproduire depuis l'environnement décrit dans [Studio CAD](../../STUDIO_CAD.md) :

```bash
PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh
QT_SCALE_FACTOR=1 python ci/studio_cad_recipe.py /tmp/studio-cad-small --width 1280 --height 844
QT_SCALE_FACTOR=2 python ci/studio_cad_recipe.py /tmp/studio-cad-hidpi
```

Les captures nécessitent un bureau OpenGL ou Xvfb. La construction autonome
est indépendante : `ci/linux_bundle.sh` extrait et teste l'archive produite
avant de la livrer. Son rapport identifie la plateforme et le binaire testés.
Une qualification Linux ne vaut pas qualification macOS.
