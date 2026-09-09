# Vinkulum Studio 0.2.0 — éditeur de mécanismes rigides 3D

Application locale PySide6 / VTK : création de corps et de liaisons, déplacement
à la souris, propriétés numériques, lois de mouvement et charges temporelles.
Le noyau Vinkulum 0.19.0 calcule dans un processus séparé ; Studio affiche ses
positions, orientations, vitesses et coordonnées de liaison.

## Installer et lancer

Depuis la racine du dépôt, avec Python **3.14**, Rust/Cargo, un compilateur C++17
et un éditeur de liens système :

```sh
python3.14 -m venv .venv-studio
. .venv-studio/bin/activate
python -m pip install .
python -m pip install './apps/studio[test]'
vinkulum-studio
```

Une roue du noyau **0.19.0 compatible avec Python et la plateforme** peut remplacer
la compilation `pip install .`. Une roue Linux ne fonctionne pas sur macOS.
PySide6 **6.11.2** et VTK **9.7.0** sont des dépendances séparées avec leurs propres
licences ; le code original de Studio est sous Apache-2.0. Aucun installateur
embarquant Qt/VTK n'est livré. Après installation, Studio fonctionne sans réseau.

Sous Linux, Qt/X11 exige notamment `libxcb-cursor0`, `libxcb-icccm4`,
`libxcb-keysyms1`, `libxcb-image0`, `libxcb-render-util0`, `libxcb-util1` et un
pilote OpenGL. Les tests automatisés utilisent également `xvfb` et `xauth`.
La qualification macOS ARM64 sur un bureau réel reste à effectuer.

## Construire et calculer

1. Choisir **Nouveau** ou un exemple : pendule G0, double pendule, bielle-manivelle,
   corps soumis à une force temporelle.
2. Ajouter une boîte, un cylindre ou une sphère. Sélectionner un corps dans
   l'arbre ou la scène, saisir ses dimensions, sa masse et sa pose dans les
   propriétés, puis appliquer. L'inertie homogène est calculée ; une matrice
   explicite peut être fournie.
3. Déplacer ou orienter le corps avec le manipulateur. Une opération terminée
   crée une entrée d'historique ; annuler/rétablir restaure le document.
4. Ajouter un pivot, une rotule, une glissière ou un encastrement entre deux
   corps, ou entre un corps et le sol. Le point et l'axe mondiaux initialisent
   deux repères locaux indépendants, ensuite éditables.
5. Configurer éventuellement la loi du pivot/de la glissière, ou une charge
   appliquée à un corps. Le dialogue propose constante, linéaire et table.
6. Corriger les diagnostics, régler durée et pas, puis **Calculer la conception**.
7. En mode **Résultat**, lire l'animation, déplacer le curseur temporel, choisir
   une courbe et exporter le CSV avec sa provenance. Revenir à **Conception**
   pour modifier le document.

La caméra permet l'orbite, le zoom, le cadrage et les vues orthographiques.
Les unités sont SI : m, kg, s, N, N·m, kg·m². Les champs d'orientation utilisent
les degrés avec la convention `Rz × Ry × Rx` ; le document stocke des matrices.
Les lois angulaires utilisent les radians. Le cylindre a son axe local Z.

**Déplacer un corps ne résout pas les contraintes.** Les ancrages incohérents
restent visibles et bloquent le lancement ; l'initialisation native refuse de
déplacer silencieusement les corps. Elle peut initialiser les vitesses imposées
par les lois de mouvement. La grille ne représente pas un contact avec le sol.

Les composantes de force et de moment sont exprimées dans le repère mondial.
Le point d'application est local au corps : son bras de levier tourne avec lui.
Les flèches orange (force) et violettes (moment) donnent une direction, leur
longueur est symbolique. Une table interpole linéairement et prolonge les
valeurs aux extrémités. Aucune expression Python n'est exécutée.

## Documents et résultats

Les objets ont des UUID persistants. Le JSON utilise
`format: vinkulum-studio-project`, `schema_version: 1` et un projet immuable.
Les fichiers G0 restent importables. Supprimer un corps conserve les références
cassées des liaisons/charges afin de pouvoir les diagnostiquer et les réparer.
La sauvegarde et les exports remplacent le fichier après écriture complète et
`fsync` du temporaire ; ils ne garantissent pas la résistance universelle à une
perte d'alimentation.

Chaque calcul capture son propre projet. Le document reste éditable pendant
l'exécution ; les propriétés et la scène du résultat restent celles du projet
capturé, en lecture seule. Une erreur, un crash ou une annulation conserve le
résultat précédent. Un seul worker est autorisé ; la fermeture le termine.

Le résultat temporaire est un NPZ borné et validé sans pickle : empreinte,
identité du projet, en-têtes avant allocation, dimensions, finitude, chronologie
et matrices de rotation. La provenance contient les versions et l'empreinte du
binaire natif. Le CSV exporte les échantillons natifs, avec l'état initial,
les unités et le projet capturé. Les pivots sont affichés en angle principal
`[-π, π]`, pas en compteur de tours. La lecture choisit les échantillons sans
interpolation dynamique. Les résultats restent en mémoire pour la session ;
le JSON sauvegarde la conception, le CSV permet de conserver les données.

## Domaine et qualification

32 corps, 64 liaisons, 128 charges, 20 000 pas demandés et 100 000 couples
corps/échantillon au maximum ; archive de résultat limitée à 64 Mio. Ces budgets
ne garantissent ni la convergence ni une limite de mémoire imposée par l'OS.
Le worker n'est pas un bac à sable de sécurité. Contact, flexibles, import CAO,
URDF, collaboration et synthèse automatique de mécanismes ne sont pas inclus.

Le statut scientifique d'une trajectoire reste **`NotAssessed`** : les contrôles
d'intégrité et les références physiques testées ne constituent pas une borne
d'erreur pour tout modèle saisi. Voir la [recette et ses preuves](../../docs/STUDIO_3D.md).

```sh
PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh
```

Sous Linux, cette commande utilise un vrai contexte Qt/X11/OpenGL sous Xvfb.
Sous macOS, la lancer depuis une session de bureau avec les droits d'accès à
l'écran. Une recette humaine reste nécessaire pour le confort, les raccourcis
et les gestes propres aux périphériques de chaque plateforme.
