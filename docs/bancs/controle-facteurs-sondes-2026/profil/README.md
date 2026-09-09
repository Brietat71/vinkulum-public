# Profil exploratoire des réponses contrôlées à 512 poutres

8 septembre 2026, dépôt `724236e96ef31f05b7bcf1203fb1f0b9d844166e`.
**Ce profil n'est pas une nouvelle campagne comparative.** Aucun solveur
concurrent n'est exécuté et les sources du dépôt restent inchangées.

Les entrées D/M natives, les six charges et les 257 fréquences 0–40 Hz
sont celles des campagnes précédentes. B/Phi et le certificat à 80 Hz
sont réutilisés depuis `n512-f40-inertie_controle-passage1`. QR, Krylov et
contrôleur sont reconstruits, puis ce même objet sert à quatre balayages :
chauffe, mesure interne sans instrumentation, cProfile, attribution par
lignes et par appels à `norme`. Il possède 40 directions Krylov, sept
coordonnées conservées et 3 072 coordonnées physiques.

CPU 8, BLAS/OpenMP/MKL/Rayon et bibliothèques complémentaires à un fil ;
venv figé Vinkulum 0.11.0. Versions, empreintes des entrées, B, scripts,
extension native et sources Python installées figurent dans
`environnement.json`. Le répertoire `sources/` est une copie exacte des
six modules locaux utilisés. `condensation_energie_distribuee.py` capture
l'implémentation installée derrière son petit module de compatibilité.

## Mesure et surcharge d'instrumentation

Le balayage sans instrumentation dure 0,3935 s ; celui sous cProfile
0,4507 s et celui sous trace 0,4477 s. Ces nombres servent à constater la
surcharge de mesure, pas à produire un nouveau ratio de solveurs.
Le temps affiché par cProfile pour son périmètre global, environ 0,485 s,
inclut aussi 0,032 s de calcul des empreintes après les réponses ; ce
travail est exclu des durées de balayage `wall_s`.

Les sorties des quatre balayages sont identiques bit à bit : champs,
coordonnées, normes et majorants par charge. Leur empreinte commune est
`9afa1c8b6434af3afbc8e914d41f813a945b65dbcc466ec87850700ce7682868`.
Ce contrôle détecte un changement de calcul dû à l'instrumentation ; il
ne remplace pas le juge physique ni une certification des réponses.

## Où passent les opérations

Attribution par lignes sous trace, sur les 257 fréquences. Les temps
des appels enfants sont inclus. Les numéros renvoient au snapshot
`sources/controle_complement.py`.

| Poste | Lignes | Durée attribuée |
|---|---:|---:|
| Deux évaluations spectrales identiques de DX, matrice 3072×7 | 158 et 187 | 0,0995 s |
| Deux corrections H appliquées à Yq, puis normes locales | 190 | 0,0802 s |
| Normes physiques des six champs, avec D et M originaux | 188–189 | 0,0660 s |
| Reconstruction X = T − W Y | 148 | 0,0488 s |
| Deux normes spectrales identiques du résidu Y, matrice 40×7 | 146 et 200 | 0,0264 s |
| Formation des Gram massique et matériel nécessaires au Schur | 150 et 152 | 0,0207 s |
| Applications D X et M X | 149 | 0,0177 s |
| Résolution du petit Y, 40×40 avec sept seconds membres | 144 | 0,0085 s |

Le chronométrage à l'intérieur de `norme` précise les doublons :
DX coûte 0,0467 puis 0,0503 s ; le résidu Y coûte 0,0115 puis 0,0130 s.
Ces durées recouvrent les lignes du tableau et ne s'ajoutent pas à elles.

cProfile confirme le diagnostic : 1 028 appels à la norme spectrale du
contrôleur, soit quatre par fréquence, et 1 285 SVD avec celle du petit
Schur. Le temps cumulé de `norme` est 0,130 s. Les 1 028 multiplications
creuses prennent environ 0,038 s avec leurs enveloppes Python. La
résolution de Y et la LU/résolution du Schur prennent ensemble environ
0,019 s. Les petits solveurs ne constituent donc pas le poste dominant.

## Gains envisageables et risques

1. **Réutiliser les deux normes répétées.** Une variable pour `norme(dx)`
   et une autre pour `norme(residu_y)` suppriment un appel de chaque type
   par fréquence sans changer l'identité mathématique. Le poste évitable
   représente environ 0,063 s sous instrumentation. Une comparaison
   bit à bit est adaptée pour vérifier cette simplification. Son gain
   réel reste à mesurer après modification.
2. **Conserver H sous forme de facteurs de petite taille.** Dans ce cas,
   H provient de C fois une matrice 1×40, mais ses images massique et
   matérielle sont conservées comme tableaux 3066/3072×40. Appliquer ces
   tableaux à Yq consomme une grande part des 0,080 s de correction locale.
   Le produit associatif `facteur_physique @ (petit_facteur @ Yq)` réduit
   fortement le nombre d'opérations lorsque le rang retenu reste faible.
   Il faut préserver les compensations dans le petit produit, conserver
   le modèle exact visé par H et contrôler le défaut des facteurs stockés.
   Le rang théorique n'autorise pas à tronquer des résidus flottants sans
   comptabiliser leur erreur. Les majorants de préparation doivent rester
   cohérents avec la représentation employée par fréquence.
3. **Réutiliser le Gram matériel pour majorer la norme de DX.** Le Schur
   demande déjà G = DXᵀ DX. En arithmétique exacte,
   `||DX||₂ = sqrt(lambda_max(G))` permettrait d'éviter les grandes SVD,
   avec un petit calcul spectral 7×7. La version machine doit tenir compte
   de l'erreur de formation du Gram et du calcul spectral : une estimation
   éventuellement inférieure n'est pas une majoration. Cette piste porte
   sur la norme d'opérateur globale ; elle ne justifie pas de remplacer
   les normes physiques des petites combinaisons de charges par un Gram
   projeté susceptible de perdre leurs compensations.
4. **Différer les transformations plus profondes de X et du Schur.**
   Les produits denses/creux physiques par fréquence restent visibles,
   mais réassocier `D(T−WY)` ou développer les formes quadratiques change
   les annulations et leurs erreurs. Ces transformations demandent des
   contre-épreuves indépendantes. Remplacer seulement le solveur 40×40
   laisse l'essentiel du coût en place.

Pour situer la priorité, la campagne précédente mesurait 0,8159 s au
total pour le contrôle contre 0,7114 s pour la LU, soit un écart de
0,1045 s. Ce contexte provient de la campagne archivée ; les temps
instrumentés ci-dessus ne sont pas directement soustractibles de ses
médianes. Les deux doublons seuls ne démontrent pas que l'écart est comblé.
Les facteurs de H et une majoration peu coûteuse de la norme de DX sont
des pistes concrètes supplémentaires, sans gain comparatif prouvé ici.

## Fichiers et reproduction

- `profiler_controle.py` : script exact du profil ; aucune modification
  des sources originales. La surcharge de `norme` existe seulement en
  mémoire pendant le dernier balayage, puis la fonction est restaurée.
- `profil.json` : statistiques complètes cProfile, attribution de toutes
  les lignes, appels de normes, durées et contrôles d'identité.
- `reponses_originales.prof`, `cprofile-cumulatif.txt`, `execution.log` :
  profil brut, affichage des fonctions et journal du lancement/arrêt.
- `bases_originales.npz`, `certificat_offert.json`, `environnement.json`
  et `sources/` : données de reconstruction et provenance.
- `manifest.json` : tailles et SHA256 de tous les autres fichiers.

Pour rejouer, copier le script dans un **nouveau répertoire**, puis utiliser
le même Python figé et les chemins d'entrée indiqués. Le script refuse
un dossier `sources/` préexistant ; ne pas relancer dans cette archive.
L'exécution originale était :

```sh
PYTHONDONTWRITEBYTECODE=1 /tmp/vinkulum-release-0.11.0-final/venv/bin/python /tmp/vinkulum-profil-controle-facteurs-2026/profiler_controle.py > /tmp/vinkulum-profil-controle-facteurs-2026/execution.log 2>&1
```
