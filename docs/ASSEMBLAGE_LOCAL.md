# Assemblage local des analyses et de la statique

Le noyau calcule désormais les tangentes d'analyse à partir des éléments,
puis assemble directement le système statique en triplets. Sur le chemin
creux, il ne construit plus une grande matrice dense pour la convertir.

Sur Princeton à 480 intervalles, le contrôle alterné mesure **11,31 s →
0,943 s**, soit **11,99×**, avec un pic mémoire médian de **191,4 → 66,5 Mio**.
La position du bout est identique entre versions et le résidu libre reste
inférieur à `3e-9 N`, en mode strict à `tol=1e-8`.

Il s'agit d'une comparaison entre deux états de Vinkulum. Le témoin est le
commit `49bdaaa2f249dc008a8e93d6747b3752d8d1644e`, avec son extension conservée
avant modification. Ces mesures ne constituent pas un nouveau classement
face à MBDyn ou Simpack.

## Calcul effectué

À état interne et multiplicateurs fixés, les matrices sont
`K = -∂f/∂q + ∂(Gᵀλ)/∂q` et `C = -∂f/∂u`. Les rotations sont perturbées
à gauche dans le repère monde. Les dérivées du gyroscopique, des couples,
des pales, des contacts et des superéléments sont assemblées localement.
Les liaisons apportent la dérivée de leurs réactions. Les poutres gardent
leurs blocs mixtes exacts et leur bloc entre rotations par différences
finies centrées.

Les contributions partageant un coefficient sont sommées dans leur ordre
de dispersion, avant l'équilibrage diagonal. Aucun seuil n'élimine les
petits termes : seul zéro est omis. L'analyse symbolique de la factorisation
est réutilisée uniquement lorsque le motif est identique. Le contrôle du
résidu linéaire et les replis pour les systèmes singuliers sont conservés.

Le coût de différenciation dépend des interactions locales, au lieu de
réévaluer le modèle entier pour chacune des `6N` perturbations. Les matrices
publiques de `k_c_m_z()` restent denses ; leur construction et leurs copies
vers Python peuvent donc encore dominer le coût de cette API.

## Trois défauts numériques exposés

- **Roulement chargé.** La dérivée de `Gᵀλ` utilisait un point matériel sur
  la roue, alors que le résidu utilisait le point de contact géométrique.
  À `h=0,05 s`, l'écart du bloc dynamique au contrôle indépendant passe
  d'environ `5,55e-3` à `4,93e-12` sur le témoin. La correction concerne
  aussi le jacobien d'intégration.
- **Inflow spatial.** Les translations des pales modifient la vitesse
  induite, même lorsque les états d'inflow sont figés. Ces colonnes
  manquaient dans le jacobien dynamique. Les trois témoins — harmoniques,
  profil radial et carte — passent d'écarts de `4 à 6e-6` à moins de
  `7e-14`. En inflow uniforme, ces colonnes sont nulles et ne sont pas
  propagées dans l'arithmétique duale.
- **Petit amortissement sous précharge.** Pour
  `τ = 10^12(1-θ) - ω`, l'amortissement vaut exactement `1`. La différence
  de deux grands couples arrondis donnait `0` ; la dérivée locale retrouve
  `1`. Le même test échoue avec l'ancienne extension.

Les lois de couple ont des cassures. À leurs bornes, la tangente utilise
la moyenne des deux pentes, conformément à la convention centrée antérieure
de l'analyse. Cela ne rend pas ces lois classiquement différentiables.
Un gouverneur de couple maximal nul est constant : sa dérivée reste nulle,
y compris à l'arrêt. Un contrôle analytique couvre cet intervalle réduit
à un point. Les autres lois par morceaux suivent leur branche active.

## Temps et mémoire

La statique Princeton est exécutée en mode strict, avec `tol=1e-8`,
`iters=100`, `paliers_max=1`, à charge pleine. Les deux versions atteignent
la même tolérance ; les trois coordonnées finales du bout sont comparées
séparément. Les résidus portent sur tous les degrés de liberté libres.

Le tableau retient le contrôle croisé : trois processus par version et par
taille, ordre des versions alterné, même CPU Linux, un fil demandé. Le temps
exclut les imports et la construction du modèle. Le pic mémoire inclut le
processus Python et ses imports ; il est lu dans `VmHWM`. `getrusage` gardait
ici le pic du parent après lancement et n'est pas utilisé pour ce classement.

| Intervalles | Avant, médiane | Après, médiane | Gain | Pic avant | Pic après |
|---:|---:|---:|---:|---:|---:|
| 240 | 2,221 s | 0,474 s | 4,69× | 92,2 Mio | 61,9 Mio |
| 480 | 11,308 s | 0,943 s | 11,99× | 191,4 Mio | 66,5 Mio |

Les coordonnées du bout sont identiques dans toutes les répétitions.
Le travail de Newton est également identique : 48 évaluations en deux
tentatives, dont 36 itérations dans la dernière tentative.

La campagne initiale, à 60, 120, 240 et 480 intervalles, présentait davantage
de dispersion sur l'ancienne version. Ses résultats restent archivés ; le
contrôle alterné évite de retenir son gain le plus favorable à 240 intervalles.
Les étendues de chaque série sont disponibles dans les données brutes.

Pour `k_c_m_z()`, trois appels après échauffement, un fil demandé, avec
construction des listes Python incluse :

| Famille | Corps | Avant | Après | Gain |
|---|---:|---:|---:|---:|
| Couples | 8 | 0,371 ms | 0,145 ms | 2,55× |
| Couples | 32 | 6,062 ms | 2,697 ms | 2,25× |
| Couples | 128 | 145,825 ms | 90,864 ms | 1,60× |
| Contacts | 8 | 0,582 ms | 0,212 ms | 2,75× |
| Contacts | 32 | 9,410 ms | 1,986 ms | 4,74× |
| Contacts | 128 | 202,269 ms | 73,771 ms | 2,74× |
| Superéléments | 8 | 1,272 ms | 0,294 ms | 4,33× |
| Superéléments | 32 | 24,477 ms | 2,536 ms | 9,65× |
| Superéléments | 128 | 459,532 ms | 71,464 ms | 6,43× |

Les quatre projections déterministes de K et C sont comparées sur ces
tailles : l'écart relatif maximal reste inférieur à `1,13e-9`. Un contrôle distinct compare
toutes les entrées aux différences finies du champ complet sur les petits
modèles déformés.

## Coût dynamique des corrections

Ces petits témoins utilisent le même pas de `0,001 s` sur `0,25 s`, un même
CPU et cinq répétitions après échauffement. Ils servent à mesurer le coût
des corrections, sans établir un classement général de performance ou un
ordre de précision temporelle. Les états finaux des deux versions sont
comparés en position, orientation et vitesse.

| Témoin | Avant | Après | Variation du temps | Itérations Newton avant → après |
|---|---:|---:|---:|---:|
| Inflow uniforme | 2,208 ms | 2,260 ms | +2,3 % | 250 → 250 |
| Harmoniques | 2,418 ms | 2,686 ms | +11,1 % | 285 → 250 |
| Profil radial | 2,723 ms | 3,059 ms | +12,4 % | 250 → 250 |
| Carte | 3,747 ms | 4,093 ms | +9,2 % | 250 → 250 |
| Roulement chargé | 2,936 ms | 2,649 ms | −9,8 % | 500 → 250 |

L'état final uniforme est identique au bit près. Pour les autres cas,
l'écart maximal par composante est inférieur à `2,4e-15 m` en position,
`1,3e-15` sur les matrices d'orientation, `1,9e-14 m/s` en vitesse et
`1,4e-14 rad/s` en vitesse angulaire. Ces durées de quelques millisecondes
ne suffisent pas à interpréter le petit écart uniforme comme une régression.

L'extension intermédiaire propagait douze directions même en inflow
uniforme, avec environ 12 % de surcoût sur ce témoin. Le choix de neuf
directions dans ce cas conserve K et C au bit près. Les champs spatiaux
gardent leurs douze directions : leur coût supplémentaire reste publié.

## Validation et portée

`ci/local.sh --bancs` passe : format, Clippy, **31 tests Rust**, **59 tests
Python**, **41 groupes de vérification**, **46 bancs rapides**, **9 bancs de
contact** et référence d'API à **185 entrées**. Les temps de campagne sous
concurrence servent à la validation, pas aux tableaux de performance.

Les **400 paliers Princeton** à 10, 20, 40 et 60 intervalles, pour les deux
formulations, atteignent encore la tolérance stricte sans stagnation.
Les contrôles analytiques d'allongement `2^-60` passent également.

Le nouveau contrôle Rust compare le champ précontraint complet à ses
différences finies, avec liaisons, bielle, vis, cardan, engrenage à porteur
mobile, roulement et indices partagés. Il vérifie une tangente non symétrique
hors équilibre. Les tests Python vérifient aussi les interactions entre
superéléments, contacts, gyroscopique et amortissement.

Le câble condensé extrême reste en échec dans les deux initialisations du
diagnostic borné. Sa validité physique n'est pas rétablie par l'accélération
de l'assemblage.

`G`, `GGᵀ`, les matrices d'analyse exposées et certains replis sont toujours
denses. Les gains mémoire de Princeton ne prouvent donc pas une complexité
linéaire sur les mécanismes fortement contraints. Les changements ne résolvent
ni les limites de la poutre intégrée, ni le couplage complet des états
aérodynamiques, ni la couverture industrielle des logiciels concurrents.

## Reproduction

Le témoin utilise le paquet Python courant, avec la seule extension compilée
au commit `49bdaaa2f249dc008a8e93d6747b3752d8d1644e`. Dans la campagne,
les fichiers Python de ce paquet étaient des liens vers le dépôt courant.
Les deux binaires exécutent ainsi les mêmes modèles et contrôles ; leurs
empreintes distinctes sont archivées.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_assemblage_local.py \
  --sortie /tmp/assemblage-actuel.json --statique 60 120 240 480

# L'ancien-pythonpath doit contenir l'extension du commit témoin.
.venv314/bin/python ci/mesure_statique_alternee.py \
  --ancien-pythonpath /chemin/vers/le/temoin \
  --sortie /tmp/assemblage-alterne.json

RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_dynamique_locale.py \
  --sortie /tmp/assemblage-dynamique.json

.venv314/bin/python ci/bilan_assemblage_local.py \
  docs/bancs/assemblage-local-avant.json \
  docs/bancs/assemblage-local-apres.json /tmp/assemblage-comparaison.json

# Le validateur lit les archives assemblage-local-* dans docs/bancs.
# Après leur collecte et une campagne ci/local.sh --bancs journalisée :
.venv314/bin/python ci/valide_assemblage_local.py \
  --journal-ci /tmp/assemblage-ci.log
```

[Campagne ancienne](bancs/assemblage-local-avant.json),
[campagne actuelle](bancs/assemblage-local-apres.json),
[contrôle alterné](bancs/assemblage-local-alterne.json),
[comparaison calculée](bancs/assemblage-local-comparaison.json),
[400 paliers stricts](bancs/assemblage-local-princeton-strict.json),
[dynamique ancienne](bancs/assemblage-local-dynamique-avant.json),
[dynamique actuelle](bancs/assemblage-local-dynamique-apres.json),
[coût intermédiaire en douze directions](bancs/assemblage-local-dynamique-dual12.json),
[tangentes Princeton](bancs/assemblage-local-tangentes-princeton.json),
[statique Princeton à chaud](bancs/assemblage-local-statique-princeton.json),
[câble](bancs/assemblage-local-cable.json),
[journal CI](bancs/assemblage-local-ci.log),
[tests Python](bancs/assemblage-local-tests-python.log),
[contrôles négatifs sur l'ancienne extension](bancs/assemblage-local-regressions-avant.log),
[comparaison uniforme neuf/douze directions](bancs/assemblage-local-uniforme.json),
[empreintes et validation](bancs/assemblage-local-validation.json).

L'[objectif généraliste face à MBDyn, Simpack et les autres références](OBJECTIF_MBDYN.md)
reste actif.
