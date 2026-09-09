# Corrections de robustesse et protocole de validation — 6 septembre 2026

Cette note accompagne l'intégration de `ccd-capsule-balayee` dans `main`.
L'analyse a été faite sur `a2bf879`, après les commits `162652b` (CCD sur
primitives), `9517cd7` (gardes et SVD bornée) et `a2bf879` (contact des capsules
et distance segment/triangle). Ces trois commits font aussi partie de
l'intégration ; leur formulation et leurs limites sont décrites dans le
[journal](JOURNAL.md).

## Défauts reproduits et corrections

| Déclencheur | Avant | Maintenant | Régression |
|---|---|---|---|
| Chute libre sur 0,01 s avec un pas demandé de 1 s | −0,2483 mm au lieu de −0,4905 mm : l'accélération initiale était limitée avec le pas nominal avant sa réduction à la durée restante | La limite utilise `min(h, t_end − t)` ; position et vitesse retrouvent la solution analytique | `verification.robustesse`, au démarrage et après une reprise, avec trois pas nominaux |
| Loi constante vide, loi linéaire tronquée, table de longueur impaire | Indexation Rust hors bornes, `PanicException` | `ValueError` à la déclaration ; longueur et finitude des paramètres contrôlées | `verification.robustesse` |
| `bloque_t=[3]` ou `bloque_r=[3]` | Liaison acceptée, puis panique à l'évaluation | `IndexError` à la déclaration ; axes autorisés : 0, 1, 2 | `verification.robustesse` |
| Orientation `diag(-1, 1, 1)` | Réflexion acceptée comme rotation | `corps` impose aussi `det(R)=+1`, en plus de l'orthogonalité | `verification.robustesse` |
| `pose_etat` avec temps NaN, données non finies ou orientation dégénérée | État invalide accepté ; un temps NaN pouvait donner une simulation vide | Validation de l'état entier avant toute mutation ; une erreur laisse le modèle intact | `verification.robustesse`, dont une erreur sur le dernier corps après une modification proposée du premier |
| Sphère d'appariement à rayon négatif ou infini, ou centre NaN | Géométrie acceptée sans erreur | `ValueError` ; centre fini et rayon fini strictement positif | `verification.robustesse`, avec déclaration valide en contrôle positif |
| Pas infini | Accepté par le contrôle `h > 0` | Refus dans l'API Python et dans `Modele::simule` | `verification.robustesse` |

Le correctif du pas initial conserve le mécanisme de limitation des
accélérations initiales sur les modèles raides. Il corrige son utilisation
d'un pas plus grand que la durée effectivement disponible ; il ne garantit
pas l'exactitude d'un pas arbitrairement grand sur tout modèle.

### Reproduire la chute libre

Depuis le dépôt et avec le paquet installé dans le venv :

```python
from vinkulum import Noyau

n = Noyau(g=[0.0, 0.0, -9.81])
n.corps("chute", 1.0, [1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0])
n.simule(0.01, 1.0)
assert abs(n.etat()[1][0][2] - (-0.0004905)) < 1e-15
```

### Compatibilité de l'API

Les appels valides conservent leurs signatures. Les entrées suivantes sont
désormais refusées plus tôt :

- Une loi `constante` exige exactement une valeur finie, une loi `lineaire`
  exactement deux. Une `table` exige au moins une paire `(temps, valeur)`
  finie ; ses temps restent strictement croissants.
- `pose_etat` accepte encore les petites perturbations d'orientation utilisées
  par Floquet et les ré-orthonormalise. Les matrices réfléchies ou dégénérées
  sont refusées. Ses tableaux de corps doivent tous avoir la taille du modèle.
- Dans `pose_etat`, `v_i=[]` laisse les inflows existants inchangés ; un tableau
  non vide doit fournir exactement une valeur finie par inflow. Les valeurs
  excédentaires ne sont plus ignorées.
- Les nouvelles erreurs de domaine sont des `ValueError`, et les indices
  d'axes hors bornes des `IndexError`. Elles sont interceptables par
  `except Exception`, contrairement aux anciennes paniques PyO3.

## Validation aérodynamique : publier toutes les hypothèses

Le score Maryland retenait, pour chaque rotor et chaque modèle de sillage,
le résultat le plus proche de la mesure parmi `N_crit=5` et `N_crit=9`.
Il mesurait donc le meilleur accord obtenu dans cette plage de paramètres.

Chaque combinaison rotor × modèle de sillage × `N_crit` a maintenant sa propre
ligne et son propre verdict. Les tolérances et les paramètres physiques sont
conservés. Le contrôle `verification.validation_scores` vérifie que changer
la mesure de référence change les verdicts sans sélectionner ni modifier les
prédictions publiées, et qu'une erreur technique reste visible.

La campagne complète historique, calculée sous Python 3.13.12 avant le
relèvement du minimum, a été recalculée en **1 959 s** : **44/62 lignes dans
leur tolérance**, contre 44/54 dans la précédente présentation. Les huit
lignes supplémentaires sont les hypothèses auparavant écartées de l'affichage
principal ; ce changement de dénominateur ne traduit pas une régression du
solveur. Les seize variantes Maryland figurent dans
[le rapport généré](VALIDATION.md). Les lignes comptent des grandeurs et des
variantes de modèles, pas 62 essais indépendants.

Exemple : la plaque cambrée sans nappe donne −2,92 % pour `N_crit=5` et
−11,87 % pour `N_crit=9`. Le second écart ne disparaît plus derrière le premier.
Ces résultats explicitent une sensibilité aux hypothèses ; ils ne constituent
pas une validation indépendante d'un paramètre ajusté sur les mêmes essais.

Le processus rend un code non nul en cas d'erreur d'exécution de la campagne.
Un écart physique hors tolérance reste un résultat publié et ne suffit pas,
à lui seul, à faire échouer la commande. Un code de sortie nul ne signifie
donc pas que toutes les prédictions sont dans leur tolérance.

## Installation et commandes de vérification

La CI utilise `PY` si fourni, sinon le Python du venv actif, sinon
`.venv/bin/python` dans le dépôt. Elle exige un venv ; `MATURIN` permet de
choisir un exécutable particulier. Elle utilise `uv` s'il est disponible,
sinon `pip` dans ce venv. L'extra `verification` est installé avec le noyau et
les sorties complètes des contrôles sont conservées à l'écran.

**Python 3.14 ou plus est obligatoire**, conformément à la consigne projet.
`requires-python` empêche l'installation sur une version antérieure, l'import
depuis les sources la refuse aussi, et la CI la refuse avant de compiler.
`.python-version` sélectionne
3.14 ; le workflow conservé dans `ci/` cible 3.14 et 3.14t. Un venv 3.13 doit
être recréé, pas simplement réactivé.

Avec Rust et `uv` installés, depuis la racine du dépôt :

```bash
uv venv --python 3.14 .venv
source .venv/bin/activate
uv pip install maturin
ci/local.sh
ci/installe_hook.sh
```

La création du venv n'est nécessaire que s'il n'existe pas encore ou utilise
une version de Python devenue incompatible. Le venv local 3.13 a été conservé
dans `.venv-python313-backup` ; le nouveau `.venv` utilise Python 3.14.7.
`ci/installe_hook.sh` configure `core.hooksPath=ci/hooks` **pour ce clone**.
Git ne transmet pas cette configuration lors d'un clone : chaque nouveau
clone doit activer son hook. Les contrôles GitHub Actions restent désactivés.
Le hook contrôle l'arbre de travail local avant le push ; il ne remplace pas
un contrôle exécuté sur chaque commit distant. `git push --no-verify` peut le
contourner.

```bash
# Contrôles supplémentaires de mécanique et de contact
ci/local.sh --bancs

# Campagne physique complète, distincte de la CI de chaque push (~33 min ici)
maturin develop --uv --release --extras validation
python -m vinkulum.validation

# Régénérer la référence d'API après une modification de son contenu
python -m vinkulum.doc --ecrire
```

La campagne réécrit `docs/VALIDATION.md` depuis les calculs ; ce rapport et
`docs/API.md` ne sont pas édités à la main. La campagne complète n'est pas
exécutée par `ci/local.sh`, y compris avec `--bancs`.

## Vérifications effectuées et portée

### Première campagne, avant le minimum Python 3.14

- `cargo fmt --check`, Clippy avec `-D warnings` et les 13 tests Rust : succès.
- CI locale complète : succès ; vérification Python en 76 s, campagne
  statique 12/12, référence d'API de 178 entrées à jour.
- Dernière garde sur les sphères : roue reconstruite, Clippy et les **48 cas
  de robustesse** vérifiés depuis cette roue avant son installation locale.
- Campagne physique complète : aucune erreur d'exécution ; les 18 écarts
  hors tolérance restent publiés dans le rapport.

Environnement du premier recalcul (historique, avant le minimum 3.14) :
Linux, Python 3.13.12, Rust 1.98.1, NumPy 2.5.2,
SciPy 1.18.1, NeuralFoil 0.3.3 et maturin 1.15.0. Les durées sont des mesures
sur cette machine, pas des seuils de réussite.

Le rapport physique complet est conservé comme trace de ce calcul sous
3.13.12. Il ne constitue pas une exécution complète de cette campagne sous
3.14. La prise en charge actuelle exige Python 3.14 ou plus, quels que soient
les interpréteurs cités dans les anciennes mesures du journal.

### Validation pour l'intégration sous Python 3.14.7

Le venv a été recréé avec `/usr/bin/python3.14`, puis les dépendances et
l'extension ont été installées pour cet interpréteur. NumPy 2.5.2, SciPy
1.18.1, NeuralFoil 0.3.3 et maturin 1.15.0 sont conservés.

| Contrôle rejoué sous 3.14.7 | Résultat |
|---|---|
| `ci/local.sh` | Succès : formatage, Clippy, 13 tests Rust, reconstruction et installation de l'extension, vérification complète et fraîcheur de l'API |
| `verification` dans cette CI | 71 s ; 48 refus typés, campagne statique 12/12 |
| `contact.capsules_primitives()` | 12 cas, pire écart géométrique 10⁻¹² ; chute et rebond de la tige cohérents entre boîte et maillage |
| `contact.tir_paroi` avec et sans CCD | Sans CCD : x final +0,030 m ; avec CCD : −0,054678 m pour boîte/cylindre, −0,036103 m pour maillage |
| Roue construite pour 3.14 | Tag `cp314-cp314`, `Requires-Python: >=3.14`, seule extension embarquée : `_vinkulum.cpython-314-x86_64-linux-gnu.so` |
| Calcul NeuralFoil NACA 4504, Re 24 000, `N_crit=5` | 53 points finis, angles strictement croissants |
| Tentatives sous le venv 3.13 sauvegardé | Import refusé ; `ci/local.sh` et `ci/installe_neuf.sh` refusent avant compilation |

Ces vérifications portent sur Python 3.14 standard. La cible 3.14t reste
déclarée dans le workflow préparé ; elle n'a pas été rejouée lors de cette
intégration. Le rapport physique historique de 1 959 s n'a pas été recalculé
intégralement sous 3.14 ; les commandes ci-dessus permettent de le faire.

Ces corrections ne ferment pas les écarts de modèle aérodynamique. Elles
n'étendent pas non plus le domaine de l'adjoint temporel, ni la couverture
géométrique du CCD des capsules. Ces limites restent celles du
[journal](JOURNAL.md) et des docstrings concernées.
