# Réduction matérielle par ports — API 0.11.0

`vinkulum.reduction_ports` construit une réponse fréquentielle réduite d'un
modèle linéaire conservatif défini par son énergie, `K = D.T @ D`, et sa masse
physique `M`. Le noyau peut fournir directement D pour la partie matérielle
de ses poutres. La préparation conserve les déformations élémentaires ;
elle n'a pas besoin de retrouver leur facteur dans une raideur déjà assemblée.

Python 3.14 ou ultérieur est requis. L'extension native ne dépend pas de
SciPy ; la couche de réduction l'utilise via l'extra `reduction`. Le projet
est privé et n'est pas publié sur un index de paquets. Depuis les sources
de la version 0.11.0 :

```sh
python -m pip install '.[reduction]'
```

Une roue locale peut également être installée en ajoutant `[reduction]`
à son chemin : `python -m pip install '/chemin/vers/vinkulum-…whl[reduction]'`,
en remplaçant ce chemin par le nom complet de la roue disponible.

## Une console native, de la création à la réponse

Cet exemple est autonome. Quatre poutres intégrées de 0,5 m relient cinq
corps, chacun portant 0,5 kg et une inertie diagonale de 0,01 kg·m². Les
rigidités sont `EA=10 000 N`, `GA=2 000 N`, `GJ=200 N·m²`, `EI=100 N·m²`.
La force de 1 N agit selon y au dernier corps. Les fréquences sont des
**pulsations en rad/s**, jamais des valeurs en Hz.

```python
import numpy as np
from vinkulum import Noyau
from vinkulum.reduction_ports import reduire_poutres

n = Noyau([0., 0., 0.])
nb, pas = 4, 0.5
for j in range(nb + 1):
    n.corps(f"section{j}", 0.5, np.diag([0.01]*3).ravel().tolist(),
            [j*pas, 0., 0.])
n.liaison("racine", None, 0)
EA, GA, GJ, EI = 10_000., 2_000., 200., 100.
for j in range(nb):
    n.poutre(f"poutre{j}", j, j+1, EA, GA, GJ, EI,
             formulation="integree")

L = nb*pas
metrique = np.diag([EA/L, EI/L**3, EI/L**3, GJ/L, EI/L, EI/L])
libres = np.arange(6, 6*(nb+1))   # supprimer les six DDL de la racine
ports = np.arange(6*nb, 6*(nb+1)) # conserver les six DDL du dernier corps
r = reduire_poutres(n, libres, ports, metrique, omega_max=0.5,
                    tolerance=1e-10, max_blocs=1, max_directions=18)

force = np.array([0., 1., 0., 0., 0., 0.])  # N puis N.m
bilan = r.adapter([0., 0.25, 0.5], force,
                  tolerance_relative=1e-8, max_etapes=8)
rep = r.reponse(0.5, force)
print(bilan["statut"])
print("uy au port [m] :", rep["deplacement_ports_physique"][1])
print("borne relative deformation :", rep["deformation"]["borne_relative"])
print("certificat spectral / champ :",
      r.certificat_spectral["certification_machine"], rep["certification_machine"])
```

Exécuté avec la roue 0.11.0 installée hors dépôt,
cet exemple donne
`tolerance_aux_frequences_demandees`, une flèche d'environ
`0.02781517265484532 m` et une borne relative en déformation proche de
`6.24e-14`. Le minorant spectral automatique vaut environ `607.93 s^-2` ;
la base finale contient 12 directions intérieures. Les derniers chiffres
des résultats flottants peuvent varier selon l'environnement.

La dernière ligne affiche `True False` : le minorant spectral possède un
encadrement des arrondis ; les bornes sur le champ restent calculées en
doubles. Ces deux statuts ont des portées différentes.

## Coordonnées, métrique et travail des forces

Chaque corps du noyau porte six DDL physiques, dans l'ordre
`(tx, ty, tz, rx, ry, rz)`. Les translations et les perturbations de rotation
sont exprimées dans le repère monde ; les rotations sont perturbées à gauche.
Pour le corps d'indice `j`, les indices globaux vont de `6*j` à `6*j+5`.

Dans `depuis_noyau` et `reduire_poutres`, `libres` et `interface` désignent
ces **indices globaux du noyau**. Leur ordre est conservé. Tous les indices
d'interface doivent appartenir à `libres` ; les autres DDL libres constituent
l'intérieur. Dans le constructeur matriciel, les indices désignent les
colonnes de D et doivent former une partition complète, disjointe et non vide.

La métrique de port H est une petite matrice SPD, dans l'ordre d'interface.
Elle fixe l'échelle physique de la tolérance. Pour sa Cholesky `H=L L.T`,
la normalisation est `W=L^-T` et le déplacement de port s'écrit `u_S=W y` :

\[
u_S^T H u_S=\|y\|_2^2,\qquad
f_S^T u_S=(W^T f_S)^T y.
\]

`reponse`, `reponses` et `adapter` prennent **f_S en unités physiques** et appliquent
eux-mêmes `W.T`. Pour un port à six DDL, il s'agit de trois forces en N
puis trois moments en N·m. Ne pas leur fournir une force déjà normalisée.
La métrique diagonale de l'exemple fixe des échelles de raideur ; elle n'est
pas annoncée comme le Schur statique exact de la console.

`r.schur(omega)` renvoie le Schur **normalisé**, approchant `W.T S(omega) W`.
`r.borne_uniforme` estime son erreur en norme opérateur 2 sur toute la bande
`0 <= omega <= omega_max`. `tolerance` vise cette erreur de matrice normalisée,
et ne représente pas directement une erreur relative sur un déplacement.
Changer H change donc la signification de cette tolérance absolue.

## Construire et lire une réduction

Les deux points d'entrée sont :

```python
from vinkulum.reduction_ports import ReductionMaterielle

# Alias de reduire_poutres(...).
r = ReductionMaterielle.depuis_noyau(noyau, libres, interface, metrique,
                                    omega_max, **options)

# Énergie et masse fournies directement, sans extraction du noyau.
r = ReductionMaterielle(d, m, interieur, interface, metrique,
                        omega_max, **options)
```

D et M doivent être réels et finis, M symétrique et positive pour définir
les normes physiques. Les couplages de masse `M_II`, `M_IS` et `M_SS` sont
conservés. D est rectangulaire ; les modes rigides sont permis avant
sélection, mais `D_I.T D_I` doit devenir défini positif lorsque les ports
sont bloqués. Les charges du modèle réduit sont appliquées aux ports ;
ni amortissement ni charge intérieure ne sont inclus dans cette API.

| Option | Rôle |
|---|---|
| `lambda_min=None` | Calcul automatique du minorant intérieur, avec certificat spectral dirigé. |
| `lambda_min=...` | Minorant positif fourni par l'appelant, en s⁻², visant les données D/M retenues. Il reste une hypothèse ; `certificat_spectral` vaut `None`. |
| `methode="qr"` | Facteur par rotations sur D ; nécessaire au certificat automatique. |
| `methode="lu_energie"` | LU équilibrée avec corrections évaluées par D ; exige un minorant fourni. Aucun basculement automatique entre méthodes. |
| `tolerance=1e-8` | Tolérance absolue du Schur normalisé. |
| `max_blocs=12`, `max_directions=128` | Limites de construction de la base intérieure. |
| `budget_qr=10_000_000` | Budget du calcul QR. |
| `precision_spectrale=80`, `budget_spectral=10_000_000` | Précision Decimal et budget spectral ; la sélection flottante et le calcul dirigé ont chacun ce budget d'opérations. |
| `tolerance_contraintes=1e-12` | Pour l'extraction native seulement : tolérance sur les valeurs de phi, dans leurs unités natives. |

Il faut `omega_max**2 < lambda_min`. Le minorant par trace peut être
conservateur : un refus de bande ne prouve pas l'existence d'une fréquence
propre dans cette bande.

Pour `rep = r.reponse(omega, force)`, les sorties principales sont :

| Clé | Signification |
|---|---|
| `deplacement_ports_physique` | Déplacement du port, dans l'ordre d'interface. |
| `champ_physique` | Champ reconstruit sur tous les DDL retenus, dans l'ordre de `libres` ou des colonnes de D. |
| `coordonnees_physiques` | Correspondance de ce champ avec les indices physiques. |
| `coordonnees_ports_normalisees`, `normalisation_ports` | y et une copie de W. |
| `borne_ports_normalises` | Majorant de l'erreur sur y en norme 2, soit l'erreur physique en norme H. Ce n'est pas une longueur en mètres. |
| `masse`, `deformation` | Dictionnaires donnant `norme_candidate`, `borne_absolue` et `borne_relative`. |
| `marge` | Plus petite valeur singulière du Schur réduit diminuée de l'enveloppe d'erreur. |
| `certification_machine` | `False` pour le contrôle de réponse. |

Les normes du champ sont `sqrt(u.T M u)` et `||D u||₂`. Une borne relative
est calculée par `b / (norme_candidate - b)` seulement si le dénominateur
est positif. Elle vise la norme de la solution exacte du modèle retenu.
Une force nulle produit un champ nul et des bornes relatives `None`.

`r.reconstruire(omega, y)` est également disponible, mais prend les
**coordonnées de port normalisées**, contrairement à `reponse` qui prend
une force physique. L'extraction et les réponses retournent des copies ;
modifier ensuite l'état du noyau ne met pas à jour la réduction déjà créée.

## Plusieurs charges à la même fréquence

Après l'exemple de console ci-dessus :

```python
forces = np.zeros((6, 3))
forces[1, 0] = 1.       # force y de 1 N
forces[2, 1] = 1.       # force z de 1 N
forces[:, 2] = forces[:, 0] - forces[:, 1]
reponses = r.reponses(0.5, forces)
assert len(reponses) == 3
np.testing.assert_allclose(
    reponses[2]["champ_physique"],
    reponses[0]["champ_physique"] - reponses[1]["champ_physique"],
    atol=1e-13, rtol=1e-12)
```

`forces` est une matrice réelle finie `(nombre_ports, nombre_charges)`,
avec au moins une colonne. Chaque dictionnaire suit le contrat de
`reponse`, sans partager ses tableaux modifiables avec les autres résultats.
La colonne nulle reste permise ; les charges complexes sont refusées.

Le Schur, sa factorisation, sa marge et les relèvements locaux sont
préparés une fois par pulsation. Les normes sont calculées sur chaque
champ reconstruit, avec `M @ champ` et `D @ champ`. Une petite combinaison
peut disparaître par arrondi dans `y.T @ Gram @ y` : cette expression
n'est pas utilisée pour la norme du champ demandé.

`reponse` et `reponses` partagent automatiquement **la dernière fréquence**.
La mémoire de cet état dépend des DDL et ports locaux, sans croître avec
le nombre de fréquences parcourues. L'enrichissement invalide l'état ; le
facteur statique et son certificat spectral sont conservés.
Plusieurs lectures simultanées d'une réduction inchangée gardent chacune
leur état de fréquence. L'appelant doit synchroniser tout enrichissement
ou changement de modèle avec ces lectures. Modifier directement D, M ou
les facteurs internes ne constitue pas une mise à jour prise en charge :
construire une nouvelle réduction pour de nouvelles données mécaniques.

Pour les assemblages, `ControleChamp.reponses(omega, forces_globales)`
partage ces opérations dans un lot. Un contrôle devient périmé après
modification des applications, masses ou raideurs externes, ou après
l'enrichissement d'un composant ; reconstruire alors `ControleChamp`.
Ce niveau n'a pas le cache automatique inter-appels de `ReductionMaterielle`.

## Adapter et reconnaître les refus

`r.adapter(frequences, force, tolerance_relative=..., max_etapes=...)`
enrichit les directions intérieures, puis réévalue les bornes des deux
normes physiques. Le facteur statique et son certificat spectral sont
réutilisés. Les résultats sont dans `reponses`, avec le même contrat
physique que `reponse`, et les étapes dans `historique`.

Les statuts sont `tolerance_aux_frequences_demandees`, `budget_etapes` et
`stagnation_ou_budget_directions`. Seul le premier signifie que les deux
bornes relatives satisfont la demande à tous les points fournis. Cette
acceptation ponctuelle n'établit pas une tolérance relative entre les points.
Le statut de construction `r.statut`, `tolerance_estimee` ou
`tolerance_non_atteinte`, concerne séparément le Schur normalisé.

La coercivité de l'intérieur ne supprime pas les résonances du système
complet. Si la marge de réponse n'est pas positive, les bornes de port et
de champ sont `None` ; un système réduit singulier peut aussi faire échouer
la résolution. Un budget épuisé ou une déflation de directions ne vaut pas
acceptation. Les entrées invalides, fréquences hors bande, forces complexes
et cas non démontrés par le certificat sont refusés explicitement.

## Contraintes et portée exacte du certificat

`depuis_noyau` accepte des contraintes fixes holonomes satisfaites à la
tolérance de phi demandée, puis exige **`G[:, libres] == 0` coefficient par
coefficient**. Aucun seuil ne supprime une petite entrée de G. Une liaison
couplant plusieurs DDL libres nécessite une application admissible construite
par l'appelant, puis le constructeur matriciel avec D et M transformés.
Les contacts non lisses, liaisons non holonomes et contraintes pilotées sont
refusés par cet adaptateur. La tolérance de phi n'est pas un certificat des
arrondis du système de contraintes.

Le certificat automatique démontre seulement

\[
0<\lambda_*\leq\lambda_{\min}(D_I^T D_I,M_{II}),
\]

pour les **valeurs binary64 d'entrée interprétées exactement**. Il utilise
une inverse sélectionnée et des opérations Decimal dirigées pour majorer
la trace et l'écart énergétique du facteur QR aux données D originales.
Ce n'est ni un certificat d'un continuum mécanique ni un certificat du
calcul flottant complet de la réponse. Les diagnostics du certificat
comprennent `lambda_min`, `eta_superieur`, `trace_superieure`, les versions
Decimal de ces bornes et les compteurs d'opérations.

Ce chemin automatique exige `M_II` exactement symétrique, SPD par blocs
connexes d'au plus six DDL, éventuellement dispersés par permutation. Les
doublons des matrices creuses sont refusés avant leur sommation flottante.
Le remplissage de la sélection et les budgets peuvent entraîner un refus,
même si une minoration existe mathématiquement. Une masse plus générale
peut être utilisée avec un minorant fourni ; sa validité et celle des
hypothèses de positivité restent alors à la charge du modèle appelant.

L'extraction native livre la contribution **matérielle des seules poutres**.
Elle reste disponible avec des déformations pondérées non nulles ; cela
n'autorise pas à remplacer la tangente précontrainte par DᵀD. La réduction
n'inclut pas les termes géométriques, réactions, couples, contacts,
superéléments, charges ou effets de mouvement absents de cette énergie.
`r.origine["donnees_natives"]` conserve z, l'état et les diagnostics de
périmètre. L'appelant choisit explicitement ce modèle matériel ; aucun
diagnostic vide ne certifie qu'il reproduit la dynamique complète du noyau.
Voir [le contrat natif et ses identités d'énergie](FACTEURS_ENERGIE_NOYAU.md).

## Coût et composition

La préparation comprend l'extraction éventuelle, le facteur statique, le
certificat automatique, les directions de réduction et l'audit uniforme.
`r.preparation_s` les comptabilise. `schur(omega)` interroge les petites
matrices réduites ; `reponse` coûte davantage, car elle reconstruit le champ
et évalue ses deux normes physiques. L'adaptation ajoute des directions et
des audits ; son bilan contient son propre `preparation_s`.

Comparer des méthodes exige donc de séparer préparation, interrogation du
Schur et réponse avec reconstruction, pour une même erreur demandée et
les mêmes sorties. Aucun gain universel sur un solveur multicorps n'est
déduit de cette API.

Pour plusieurs sous-structures, `AssemblagePorts` prend une application
linéaire **de déplacement physique** par composant, puis effectue la
normalisation locale. `ControleChamp` et `adapter_champ` sont disponibles
pour ce niveau d'assemblage. Les masses et raideurs externes doivent définir
des énergies positives pour ces normes. Ces objets restent des modèles
fréquentiels linéaires ; ils ne créent ni superélément dynamique natif ni
intégrateur réduit non linéaire.

Les contrats exécutables sont dans
[test_reduction_ports.py](../python/vinkulum/test_reduction_ports.py),
[test_reponses_groupees.py](../python/vinkulum/test_reponses_groupees.py) et
[test_inverse_selectionnee.py](../ci/test_inverse_selectionnee.py).
