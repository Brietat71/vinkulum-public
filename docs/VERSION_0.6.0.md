# Vinkulum 0.6.0 — spectre complet et diagnostics fiables

Ce lot étend les analyses mécaniques et corrige les contrôles d'utilisation.
L'intégrateur, les formulations de contact et les modèles aérodynamiques
conservent leur périmètre. Aucune dépendance d'exécution n'est ajoutée.

## Spectre et bilan

`Noyau.spectre(t=None)` fournit toutes les racines de la linéarisation
mécanique locale, comme couples `(réel, imaginaire)` en s⁻¹. Les multiplicités
et les deux membres des paires conjuguées sont conservés ; le classement est
par partie réelle décroissante, puis partie imaginaire croissante. Le calcul
utilise la cinématique tournante existante et résout la masse réduite par
Cholesky, sans l'inverser explicitement.

Sur `x'' + 5x' + 4x = 0`, les racines sont −1 et −4. L'ancienne API
`modes_complexes()` rend une liste vide puisque rien n'oscille ; `spectre()`
rend les deux racines. Sur `x'' − 4x = 0`, les racines −2 et +2 révèlent la
divergence. `modes_complexes(combien=6, t=None)` conserve ses triplets
`(fréquence Hz, ζ, σ)`, son filtre oscillant et son tri par fréquence.

`Noyau.bilan_stabilite(t=None, tol=1e-8)` rend un dictionnaire :

| Champ | Signification |
|---|---|
| `t` | Date analysée ; défaut : horloge courante, sans mutation |
| `racines` | Couples du spectre complet |
| `abscisse_spectrale` | Maximum des parties réelles ; `None` sans spectre |
| `tol` | Tolérance absolue finie et ≥0, en s⁻¹ |
| `statut` | Un des cinq résultats ci-dessous |
| `limites` | Liste des limites d'interprétation et motifs de hors domaine |

- `croissance_detectee` : au moins une partie réelle > `tol`.
- `decroissance_spectrale` : toutes les parties réelles < `-tol`.
- `marginal_ou_indetermine` : aucun des deux cas précédents, notamment modes
  rigides, oscillations non amorties et racines proches de l'axe imaginaire.
- `sans_ddl` : espace admissible vide, y compris modèle vide.
- `hors_domaine` : contacts, candidats à l'appariement ou appariement
  configuré, pales ou inflows présents. Le bilan rend `racines=[]` et
  `abscisse_spectrale=None` ; `spectre()` refuse ces cas par `RuntimeError`.

Une erreur numérique produit toujours une exception, jamais un résultat
« sans ddl » ou un statut rassurant. Les dates non finies sont refusées,
y compris hors domaine. Une tolérance invalide produit `ValueError`.

**Le bilan concerne une linéarisation locale**, sans vérifier que l'état est
un équilibre. Il ne certifie pas la stabilité globale d'une trajectoire.
Pour un système périodique, employer Floquet. Le spectre reste dense, les
tangentes mécaniques conservent leurs approximations documentées, et la
tolérance est un seuil de lecture, pas une borne d'erreur sur les racines.
Les racines multiples, notamment à l'amortissement critique, sont sensibles
aux petites erreurs des tangentes ; elles ne sont ni fusionnées ni arrondies
artificiellement. Les APIs historiques conservent leurs domaines antérieurs.

## Diagnostics et migration

`controle(t=None)` et `phi_dot(t=None)` emploient le même diagnostic
`Φ̇ = G·u + ∂Φ/∂t`, évalué à état fixé. Le contrôle omettait le terme temporel
et pouvait signaler une vitesse incompatible sur une liaison pilotée bien
suivie. Les évaluations aux dates perturbées sont vérifiées. Aux cassures
des lois tabulées, la différence centrée reste une approximation : aucune
dérivée classique n'y est garantie.

Le défaut de `controle()` devient **la date courante**, au lieu de zéro.
Pour contrôler explicitement à zéro, utiliser `controle(t=0)`.
Le retour reste une liste de textes ; la formulation de certains messages
change. Les diagnostics laissent inchangés l'horloge, les poses, le schéma,
les multiplicateurs et les compteurs.

Les corps reliés par superéléments, couples, contacts, candidats à
l'appariement, pales ou efforts sont pris en compte dans la recherche des
corps isolés. Ce diagnostic décrit les interactions déclarées ; il ne prouve
pas qu'elles bloquent effectivement un mouvement ou exercent un effort.
Le contrôle des inerties utilise une décomposition bornée et normalisée.

## Étude de convergence

Les formats de `convergence.etude()` et `verdict()` sont conservés. Sont
désormais refusés : entrées non finies, pas non positifs ou non représentables,
facteurs non positifs/non croissants, moins de deux facteurs, réutilisation
du même noyau, dates initiales différentes, observables vides, non finies ou
de dimensions variables. Les résultats d'observables sont copiés pour éviter
qu'un tampon réutilisé efface les différences. Une simulation qui n'atteint
pas la date finale demandée produit `RuntimeError`.

L'ordre n'est calculé que pour trois pas successifs de rapports égaux
(tolérance relative `1e-12` sur leurs logarithmes), avec écarts non nuls.
Des facteurs comme `(1, 2, 5)` restent utilisables pour les écarts, avec
`ordre=None`. Deux facteurs donnent également un ordre indéterminé.

Le verdict mesure l'accord du premier calcul avec le plus fin, sur
l'observable finale choisie et avec une tolérance absolue. Il ne garantit
ni une erreur absolue par rapport à la solution exacte, ni la validité du
modèle physique. La fonction de construction doit restituer le même problème ;
les contrôles d'identité et de date ne démontrent pas cette équivalence.
Les noyaux construits sont conservés jusqu'à la fin de l'étude pour détecter
leur réutilisation ; la mémoire augmente donc avec le nombre de raffinements.

## Validation et coût

Les 35 régressions Python passent : 24 historiques et 11 nouvelles méthodes
de test, plusieurs paramétrées. Elles couvrent les oscillateurs analytiques
sous-amorti, suramorti, critique et divergent, les modes rigides, la toupie
tournante, les diagnostics temporels et la convergence contre un oscillateur
exact. Le test de débordement en sous-processus couvre aussi les nouvelles
APIs. La CI locale étendue (`ci/local.sh --bancs`) a réussi : formatage,
Clippy sans avertissement, **20/20 tests Rust**, vérification physique
complète en **73 s**, **12/12 cas statiques convergents**, référence d'API à
jour (**183 entrées**), bancs étendus et suite dédiée aux contacts/CCD.
La garde supplémentaire sur les pas raffinés arrondis à la même valeur a
ensuite été contrôlée par les 11 tests d'analyse et la vérification d'API.

La référence 0.5.2 reproduit les défauts : contrôle d'une liaison suivie
annonçant à tort une erreur de pose de 0,1 m et de vitesse de 0,5 m/s,
réutilisation d'un même noyau acceptée dans l'étude de convergence, ordre
annoncé à 1,388 pour les facteurs irréguliers `(1, 2, 5)`. Ces cas sont
couverts par les nouvelles régressions.

Le [banc spectral](../ci/banc_spectre.py) compare une console de 30 corps
contre l'extension v0.5.2 conservée avant reconstruction. Il mesure
`modes_complexes(combien=12)` : tangentes, masse, base et Schur, construction
exclue. Huit processus isolés suivent l'ordre ABBA à un puis quatre threads
Rayon ; chacun réalise un échauffement puis sept répétitions. BLAS/OpenMP
restent à un thread. Les empreintes des extensions et les 56 mesures sont
dans le [relevé brut](bancs/spectre-0.6.0.json).

| Threads Rayon | Médiane v0.5.2 → v0.6.0 | Variation | MAD avant / après |
|---|---:|---:|---:|
| 1 | 79,307 → 79,181 ms | −0,16 % | 0,472 / 0,691 ms |
| 4 | 78,981 → 79,766 ms | +0,99 % | 0,489 / 0,511 ms |

Les sorties historiques sont identiques sur ce cas (écart maximal nul,
comparaison à `rtol=atol=1e-8`). Le petit ralentissement à quatre threads est
publié ; aucun gain de vitesse n'est revendiqué. Ce lot vise la fiabilité
des analyses et ce seul cas ne décrit pas les performances générales.

```bash
python -m unittest vinkulum.test_noyau vinkulum.test_analyses -v
ci/local.sh --bancs
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ci/banc_spectre.py
```
