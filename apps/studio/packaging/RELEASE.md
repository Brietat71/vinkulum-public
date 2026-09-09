Vinkulum Studio 0.4.0 — conception CAD et mécanismes 3D.

Cette version ajoute OCCT 8.0.1 et build123d 0.11.1 adapté : primitives,
extrusions, opérations booléennes, congés sur toutes les arêtes et échange STEP
pour une pièce solide. Masse, centre de masse et inertie sont calculés sur le
BREP ; les projets conservent la géométrie exacte et une représentation de rendu.
Une opération CAD s'exécute dans un processus séparé et son échec conserve le
document. Les esquisses contraintes interactives et l'arbre de fonctions
régénérable ne sont pas encore disponibles.

L'interface propose les ateliers Modéliser / Simuler / Examiner, un inspecteur
plus dense, des nombres compacts conservant leur valeur exacte, les champs
X/Y/Z et les outils de sélection / déplacement / rotation. L'examen conserve
les entrées des calculs comparés et les grilles temporelles natives.

Le noyau reste Vinkulum 0.19.0. Python 3.14, Qt/PySide6 et VTK sont embarqués
avec le moteur CAD. Les versions et conditions d'installation sont précisées
dans les notices de la plateforme. Chaque distribution doit passer les tests
de son propre exécutable extrait : opération CAD, aller-retour STEP, double
pendule, rendu OpenGL et identification des versions. Le rapport joint décrit
la plateforme effectivement testée ; ce n'est pas une certification générale.

La première qualification de cette version est locale sous Linux x86-64,
Ubuntu 24.04 / glibc 2.39 / X11. La recette Apple Silicon est préparée ; une
qualification Linux ne vaut pas qualification macOS. Si un DMG est joint à
cette release, consulter son rapport macOS, sa signature et ses instructions
spécifiques. La signature prévue est ad hoc, sans notarisation Apple.

Pinocchio et les futurs connecteurs de solveurs sont décrits dans la feuille
de route ; ils ne sont pas embarqués dans cette version.

Code original Apache-2.0 ; licences et adaptations tierces incluses. Sources
exactes : commit ou tag associé à la distribution et provenance de construction.
