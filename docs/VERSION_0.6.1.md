# Vinkulum 0.6.1 — campagnes et gates parallèles

Les petits modèles passent peu de temps dans chaque opération numérique.
Les distribuer **entre processus indépendants** permet d'exploiter plusieurs
cœurs sans multiplier les synchronisations à l'intérieur d'un petit Newton.
Ce correctif change l'exécution des campagnes ; les modèles, leurs paramètres,
assertions physiques et résultats publics sont conservés.

## Exécution

Les mêmes listes de cas pilotent les modes séquentiel et parallèle :

| Commande | Groupes indépendants | Mode par défaut |
|---|---:|---|
| `python -m vinkulum.verification` | 41 | Parallèle |
| `python -m vinkulum.bancs rapide` | 46 | Parallèle |
| `python -m vinkulum.contact` | 9 | Parallèle |
| `python -m vinkulum.bancs` | 46, avec tous les raffinements | Séquentiel pour les mesures travail–précision |

La CI locale et le hook pre-push appellent ces campagnes et bénéficient
automatiquement de la parallélisation. Le modèle GitHub Actions active aussi
la concurrence sur la campagne complète, qui conserve tous ses raffinements.
Les builds, la vérification d'API et les phases successives de CI restent
ordonnés selon leurs dépendances.

Chaque cas dispose d'un processus Python neuf. Ses imports scientifiques se
font après le réglage de son affinité CPU sous Linux. Le budget est réparti
entre les processus, Rayon est borné à la part attribuée et BLAS/OpenMP à un
thread. Les sous-processus d'un test héritent de cette affinité.
`RAYON_NUM_THREADS`, si cette variable est fixée, peut borner davantage
le nombre de threads. Un cas composite déjà dans un worker exécute ses
sous-cas sur place : aucun pool
supplémentaire n'est lancé par l'ordonnanceur.

Quatre contrôles fondés sur des temps restent **exclusifs**, après les autres
cas : `verification.echelle`, `verification.refus_gardes`,
`contact.appariement` et `contact.maillages`. Ils récupèrent le
budget CPU complet. Les autres durées mesurées sous concurrence décrivent
la validation, et ne doivent pas servir à comparer la vitesse isolée des
modèles. Une archive de bancs précise `mesures_concurrentes`.

## Réglages et journaux

| Variable | Défaut et effet |
|---|---|
| `VINKULUM_BANCS_JOBS` | Jusqu'à 8 processus, bornés par le nombre de cas et le budget CPU ; `1` pour le séquentiel |
| `VINKULUM_BANCS_CPUS` | CPU accessibles au processus ; peut réduire le budget total |
| `VINKULUM_BANCS_TIMEOUT` | 900 secondes par cas, finies et strictement positives |
| `VINKULUM_BANCS_LOGS` | Répertoire parent des journaux ; défaut : répertoire temporaire du système |

```bash
VINKULUM_BANCS_JOBS=8 VINKULUM_BANCS_CPUS=16 ci/local.sh --bancs
VINKULUM_BANCS_JOBS=1 python -m vinkulum.verification
VINKULUM_BANCS_LOGS=target/campagnes python -m vinkulum.contact
```

Les CPU demandés sont plafonnés à ceux accessibles, sans supposer que le
nombre de CPU de la machine est celui attribué au processus. Hors Linux,
les limites de threads restent appliquées mais l'affinité matérielle n'est
pas imposée. Le budget CPU ne constitue pas une limite mémoire : réduire
le nombre de processus pour des campagnes qui consomment beaucoup de RAM.

Le répertoire de chaque campagne contient un journal par cas et un
`rapport.json` : cas prévus, cas terminés, statuts, durées, CPU attribués et
motifs d'échec. Les messages de fin de cas s'affichent immédiatement, sans
entrelacer les journaux. Les journaux sont conservés pour le diagnostic.
Les fichiers de résultats internes sont propres à la campagne.

Un cas échoué, tué, dépassant son délai ou ne produisant pas son résultat
fait échouer la campagne et donc la gate. Les autres cas déjà programmés
sont également collectés. Les expirations et interruptions arrêtent le
processus concerné et, sous POSIX, son groupe de sous-processus.

## Validation et mesures

Les tests de l'ordonnanceur vérifient la même couverture en séquentiel et en
parallèle, la conservation des résultats, le chevauchement des cas éligibles,
l'exclusivité des mesures, les budgets CPU/threads, les journaux et la
propagation des échecs et délais dépassés. Ils font partie de la gate.

La vérification précédente v0.6.0 avait réussi en 73 secondes. Le nouveau
découpage ajoute les trois tests de l'ordonnanceur et le coût d'un processus
neuf par groupe. Sa comparaison contrôlée utilise le même corpus de 41
groupes, sur 16 CPU accessibles, sans compilation ni autre campagne en parallèle :

| Exécution | Durée de la campagne | Groupes réussis |
|---|---:|---:|
| 1 processus, budget 16 CPU | 89,33 s | 41/41 |
| 8 processus, budget total 16 CPU | 24,72 s | 41/41 |

Soit **3,61×** entre ces deux exécutions du même découpage. La comparaison
avec les 73 s historiques inclut le changement d'isolation et les nouveaux
tests. Il s'agit d'un passage par configuration sur cette machine ; le gain
porte sur la campagne, pas sur le temps d'intégration d'un modèle isolé.
Le [relevé brut](bancs/campagnes-0.6.1.json) contient les cas, leurs statuts,
durées, affectations CPU et l'empreinte de l'ordonnanceur.

La CI étendue a réussi : formatage, Clippy sans avertissement, 20 tests Rust,
38 tests Python, référence d'API à jour (183 entrées), vérification 41/41,
bancs rapides 46/46 et contacts complets 9/9. Les deux dernières campagnes
ont pris respectivement 171,70 s et 144,20 s. Quelques cas longs dominent
encore leur fin ; ces durées ne constituent pas une mesure de gain sans
un passage séquentiel correspondant.

La revue finale a ajouté l'exclusivité de `refus_gardes` et `maillages` :
le premier contrôle une part de temps interne au solveur, le second
l'exposant du coût des requêtes BVH. La vérification complète est rejouée
avec ce réglage dans la comparaison ci-dessus ; le contrôle des maillages
est aussi rejoué isolément. Le hook pre-push exécute la CI de livraison.
