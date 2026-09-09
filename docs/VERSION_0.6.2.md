# Vinkulum 0.6.2 — multi-rythme restreint aux partitions

Les micro-pas factorisaient encore le système global, avec des lignes identité
pour les corps gelés. Newton et l'accélération initiale factorisent maintenant
uniquement les corps mobiles et les contraintes actives de leur partition.
L'API, les pas, les tolérances et le schéma de couplage sont conservés.

## Résolution et corrections

La correspondance entre indices globaux et locaux conserve six composantes
par corps mobile, puis les multiplicateurs actifs. Les lignes éliminées sont
résolues explicitement : leurs valeurs, éventuellement non nulles, contribuent
au second membre local. Les blocs de masse gelés de l'accélération initiale
sont résolus directement. Une partition vide ne lance aucune factorisation.
Le seuil dense/creux s'applique à la dimension réduite.

Chaque partition conserve son historique d'intégration et son cache de
jacobien. La rotation des vecteurs utilise les indices globaux, même pour des
corps mobiles non contigus. Les caches tiennent compte de la correspondance,
du masque, de la dimension et du motif. Pour les micro-pas seulement, une
variation relative du pas inférieure à `1e-10` n'invalide pas le cache : cela
absorbe les arrondis des dates interpolées, sans modifier le pas intégré.
Le contrôle de descente de Newton reste le juge de réutilisation.

La restriction seule ne suffisait pas : sur 16 maillons, une première mesure
passait de 3,79 à 3,37 s, dont 2,33 s de jacobien. Les caches par partition
suppriment ce recalcul systématique ; les calculs de résidu entièrement gelés
et les logarithmes de rotation répétés ont également été retirés. Les vecteurs
d'état restent globaux, donc une partie du coût croît encore avec le modèle.

Deux corrections accompagnent cette restriction :

- Les réactions lentes étaient remplacées par les zéros du dernier
  micro-pas rapide dans les trajectoires et `reactions()`. Elles proviennent
  désormais du dernier pas de leur propre partition. Un corps de 1 kg
  suspendu à une liaison lente conserve ainsi sa réaction de pesanteur.
- Une nouvelle détection de redondance conserve les lignes gelées masquées,
  et les masques recalculés de chaque partition survivent au macro-pas.

Les caches locaux sont abandonnés après un échec. La sauvegarde externe
restaure le dernier macro-pas accepté, ses réactions, ses masques et son
historique ; les sauvegardes intermédiaires redondantes ont été retirées.

## Mesures isolées

Même modèle de chaîne et ensemble vibrant à 300 Hz, `T = 0,3 s`, `h = 1 ms`,
`k = 100`. Python 3.14.7, NumPy 2.5.2, SciPy 1.18.1 ; Rayon, BLAS et OpenMP
à un thread, même affinité accessible. Chaque configuration fait un
échauffement puis cinq répétitions, sans compilation ni autre campagne en
parallèle. La référence est la roue v0.6.1 reconstruite depuis `181adad`.

| Maillons lents | Dimension globale → lente / rapide | Multi 0.6.1 | Multi 0.6.2, médiane [min–max] | Gain sur 0.6.1 | Monolithique fin 0.6.2 | Gain sur le fin |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 100 → 88 / 12 | 4,310 s | 1,281 s [1,278–1,281] | 3,36× | 1,313 s | 1,03× |
| 16 | 188 → 176 / 12 | 3,791 s | 1,637 s [1,625–1,642] | 2,32× | 1,689 s | 1,03× |
| 32 | 364 → 352 / 12 | 5,012 s | 2,316 s [2,303–2,317] | 2,16× | 3,434 s | 1,48× |

Le gain sur le fin est faible aux deux petites tailles ; ce tableau ne
promet pas une accélération universelle. Les dimensions réellement
factorisées sont contrôlées par des tests Rust, y compris de part et d'autre
du seuil dense/creux. Les jacobiens passent de 30 300 à 43 / 43 / 36 sur les
trois tailles. Les positions finales diffèrent de la v0.6.1 de moins de
`4e-10 m`, les vitesses de moins de `1e-8 m/s`. Les réactions lentes changent
intentionnellement, puisque leurs zéros étaient erronés.

La réutilisation fait davantage d'itérations de Newton sur ce modèle :
environ 177 000 contre 30 300. Elles coûtent moins cher que les jacobiens
évités, mais leurs résidus et vecteurs globaux limitent le gain, notamment
sur les petites chaînes. Les tolérances n'ont pas été relâchées.

Le script de mesure échantillonne toutes les trajectoires à 1 ms pour
comparer des coûts de sortie identiques. Ses amplitudes à cet échantillonnage
ne remplacent pas le banc physique : celui-ci échantillonne la référence fine
à 0,1 ms. Sur ce banc, l'amplitude multi-rythme vaut **95,9 %** de la référence,
contre **80,6 %** au pas grossier ; l'écart du bout de chaîne est **1,91 mm**
sur 2 m. Les critères de précision existants sont conservés.

Les [mesures archivées](bancs/multirythme-0.6.2.json) contiennent les temps de
chaque répétition, les chronomètres par phase, les états finaux déterministes,
les erreurs et les empreintes des sources et extensions.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_multirythme.py --sortie target/multi.json
# PYTHONPATH peut sélectionner une autre roue extraite ; le modèle reste commun.
```

## Validation et domaine

Les régressions couvrent l'élimination avec second membre non nul, les
rotations sur des corps non contigus, le repli SVD, les contraintes
redondantes, les partitions vides, `k=1`, les fins raccourcies, les changements
de partition et les échecs dans chacune des deux phases, y compris après
des macro-pas acceptés. Les réactions sont comparées au monolithique sur des
pendules indépendants.

La CI locale étendue a réussi : formatage, Clippy sans avertissement,
**25 tests Rust**, **41 tests Python**, API à jour (**183 entrées**),
vérification **41/41**, bancs rapides **46/46** et contacts **9/9**.
Les trois campagnes ont pris respectivement 22,32 s, 179,79 s et 167,88 s ;
ces durées sous concurrence ne sont pas des mesures de vitesse isolée.
`ci/installe_neuf.sh` a également réussi : roue reconstruite, installation
dans un venv vierge, vérification 41/41 hors dépôt, chargement des données
embarquées et bancs rapides exécutés depuis le répertoire de l'utilisateur.

Le couplage reste l'extrapolation des rapides pendant le pas lent puis
l'interpolation des lents pendant les micro-pas. Aucun correcteur ni ordre
supplémentaire n'est annoncé. Les liaisons traversantes restent refusées ;
l'aérodynamique, le contact non lisse, GGL et l'adaptatif restent hors domaine.
Il n'y a aucune migration d'API.
