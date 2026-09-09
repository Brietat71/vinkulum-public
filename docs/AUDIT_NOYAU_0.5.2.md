# Audit du noyau et de l'API — v0.5.2, 6 septembre 2026

Suite de l'[audit v0.5.1](AUDIT_NOYAU_2026-09-06.md), sur la base Git
`45b4a7f`. Le périmètre comprend mécanique, contraintes, intégrateurs,
contacts, analyses et interface Python. Les modèles aérodynamiques ne sont
pas réévalués contre des essais ; leurs bancs existants servent de contrôle
de non-régression aux changements de gestion d'état.

Les tests historiques ne suffisaient pas : les 17 tests Rust et les 12
régressions Python passaient avant les corrections ci-dessous. Les nouveaux
tests confrontent aussi le noyau à des références indépendantes : cinématique
fermée, équilibre d'un encastrement, dérivées d'un servo, SVD SciPy et
différences de modèles au paramètre perturbé.

## Défauts reproduits et corrections

| Priorité | Déclencheur et résultat initial | Correction et contrôle |
|---|---|---|
| P1 | Encastrement déplacé de **0,1 m**, masse `1e10`, gravité `−9,81` : `statique()` rend `(0.1, 0)` et laisse le déplacement | Le mérite sépare `‖r_force‖/échelle_force` et `‖Φ‖`. La grande charge n'autorise plus une grande erreur géométrique ; contrôle aussi des réactions |
| P1 | Servo à consigne `0,5t`, modèle à `t=0`, analyse à `t=1` : équilibre à **0 rad** au lieu de **0,5 rad** | Forces, accélération consistante et tangentes utilisent la date d'analyse. Le test dérive également la saturation `tanh` et la limitation de vitesse |
| P1 | Rotule dont le bras vaut `1e5` : base admissible de dimension **4 au lieu de 3**, `‖GZ‖≈1` | QR pivoté sur `Gᵀ` après normalisation des lignes ; complément orthogonal complet, sans former `GᵀG`. Comparaison des projecteurs à SciPy sur plusieurs échelles et avec contraintes répétées |
| P1 | Engrenages de rapports **1,5 et −2,3**, vitesses constantes : refus vers **1,59 s** pour un faux saut d'angle | Le recentrage par tours n'est appliqué qu'aux rapports entiers. Contrôle des orientations et vitesses analytiques sur 20 s, rapports rationnels et irrationnel, en série et sur quatre instances concurrentes |
| P1 | Sensibilité d'une poutre avec un superélément déformé et un contact : une composante attendue à **0,1** vaut **10,1**, et des efforts apparaissent sur un nœud étranger à la poutre | Évaluation de la seule poutre paramétrée. `d_raideur_poutre` assemble sa tangente locale par duaux emboîtés. Vérification par différences de modèles complets |
| P1 | Masse `1e-300`, superélément de raideur `1e100`, données finies : débordement modal suivi d'une **PanicException** | Contrôle des matrices avant décomposition, décompositions propres/Schur bornées, validation des valeurs propres ; `RuntimeError` interceptable. Reproduction isolée avec délai maximal |
| P2 | Simulation jusqu'à **5e-13 s** : succès sans progression, temps restant à zéro | Les trois intégrateurs avancent jusqu'au temps final représentable. Pas d'arrêt par soustraction d'un seuil absolu ; pas non représentable refusé avec restauration |
| P2 | `modes(0)` rend un mode ; dates NaN acceptées sur certains chemins modaux | Liste vide pour zéro mode, contrôles de date et d'état avant retour ; modèle vide couvert |
| P2 | `bloque_r=[0,1,2,0]` accepté, puis chemin pouvant paniquer dans le calcul du bassin | Axes répétés refusés à la déclaration, avant insertion de l'élément |
| P2 | Un contact sans `b`, avec un maillage porté par le corps en contact, contourne le refus de contact avec soi-même | Vérification répétée après résolution du porteur implicite du maillage |
| P2 | Axes finis dont la norme déborde, limites PD nulles, paramètres de contact non finis et inertie minuscule relativement asymétrique acceptés | Validation des domaines avant mutation ; symétrie relative de l'inertie et Cholesky sur tenseur normalisé ; table Rust refusant NaN/Inf |

Deux durcissements complètent ces reproductions : les erreurs de contrainte
sur les poses perturbées sont propagées par `raideur`, et la sauvegarde de
l'état inclut les références d'énergie et de moment des projections.

## Architecture, API et optimisation

Les signatures et formats de sortie sont conservés. Ce lot est une version
corrective **0.5.2**, avec numéros Rust/Python synchronisés. Les arguments
invalides nouvellement détectés produisent `ValueError` ; les échecs d'analyse
produisent `RuntimeError`. Les dates d'analyse ne déplacent pas l'horloge du
modèle. Le résidu retourné par `statique` garde son format historique.

Les opérations numériques communes sont isolées dans `src/numerique.rs`.
Une extraction générale du monolithe n'était pas nécessaire pour ces
corrections et aurait compliqué leur comparaison physique.

Optimisations du lot :

- Deux copies complètes des corps et une copie des inflows supprimées dans
  `actualise` ; les données sont empruntées jusqu'à leur dernière lecture.
- Les sensibilités de poutre n'évaluent plus tous les éléments pour toutes
  les colonnes et ne copient plus le modèle complet.
- La base admissible utilise les réflecteurs QR, évitant le produit normal
  et sa décomposition spectrale ; elle reste une analyse dense.

Le banc `ci/audit_noyau.py` conserve les empreintes des extensions, les
répétitions, médiane, dispersion MAD, pic RSS et résultats physiques. Il
couvre les seuils dense/creux (17/18 corps), le seuil de parallélisme creux
(221/223 corps), 900 corps, poutres, superéléments, contacts et analyses.
Le [relevé brut](bancs/audit-noyau-0.5.2.json) contient **96 processus et
672 mesures**, sans compilation ni autre suite de tests concurrente. La
référence a été reconstruite depuis `git archive 45b4a7f`, puis sa roue
extraite dans un répertoire indépendant ; les deux extensions sont mesurées
sous le même Python 3.14.7. BLAS/OpenMP sont bornés à un thread ; Rayon est
réglé à un puis quatre threads.

Médianes des 14 répétitions mesurées par version/configuration, **temps total
de l'appel chronométré**, en millisecondes. Les chaînes sont intégrées sur
20 pas ; `modal30` mesure `k_m_z`, donc K, M et la base, et non les seules
valeurs propres. Les autres durées sont définies dans le script.

| Cas | 1 thread, avant → après (ms) | Variation | 4 threads, avant → après (ms) | Variation |
|---|---:|---:|---:|---:|
| Pendule | 45,361 → 42,495 | −6,3 % | 42,446 → 41,322 | −2,6 % |
| Chaîne 17 | 2,927 → 2,922 | −0,2 % | 3,011 → 2,960 | −1,7 % |
| Chaîne 18 | 2,676 → 2,618 | −2,2 % | 2,749 → 2,658 | −3,3 % |
| Chaîne 30 | 5,106 → 4,991 | −2,3 % | 5,148 → 5,049 | −1,9 % |
| Chaîne 221 | 247,767 → 254,251 | **+2,6 %** | 240,947 → 248,011 | **+2,9 %** |
| Chaîne 223 | 254,339 → 256,062 | +0,7 % | 244,689 → 251,216 | **+2,7 %** |
| Chaîne 900 | 175,373 → 172,789 | −1,5 % | 149,770 → 150,357 | +0,4 % |
| Poutre 8 | 142,769 → 142,340 | −0,3 % | 144,765 → 141,796 | −2,1 % |
| Superélément 8 | 1,919 → 1,878 | −2,1 % | 1,921 → 1,906 | −0,8 % |
| Contact | 0,334 → 0,338 | +1,2 % | 0,335 → 0,332 | −0,7 % |
| K/M/base, 30 corps | 9,201 → 7,280 | **−20,9 %** | 9,200 → 7,296 | **−20,7 %** |
| Sensibilité K, 30 corps | 1,500 → 1,080 | **−28,0 %** | 1,503 → 1,069 | **−28,8 %** |

Les gains K/M/base et sensibilité sont reproduits dans les deux ordres et
dépassent leur dispersion MAD (0,4 à 1,9 %). Les gains apparents du pendule
ne sont pas robustes à l'ordre : la première comparaison à un thread donne
−46,1 %, la seconde +2,3 %. Ils ne sont pas retenus comme preuve de gain.
Les autres petits écarts doivent être lus avec leurs dispersions et les
ralentissements publiés. **Aucune accélération générale n'est démontrée.**

Le profil interne confirme néanmoins une réduction du coût de fin de pas :
sur 900 corps à un thread, le cumul passe de **3,812 à 1,903 ms** pour les
20 pas ; le nombre d'itérations et de jacobiens est identique. La résolution
linéaire et l'assemblage dominent toujours. Pour 221 corps, environ 200 ms
sont dépensées hors de la boucle de pas, principalement dans l'initialisation
avec QR dense de détection des redondances. Cette détection est sautée à
900 corps, ce qui explique pourquoi le banc court à 900 corps peut être plus
rapide que celui à 221 corps. Les chronos ne permettent pas d'attribuer
précisément le petit ralentissement d'initialisation au changement de code.

Le pic RSS à 900 corps reste proche de **107 000 Kio à un thread** dans les
deux versions ; aucune baisse générale de mémoire résidente n'est démontrée.
Les comparaisons physiques de toutes les répétitions passent (`rtol=1e-8`,
`atol=1e-9` pour états, énergie, contraintes et sorties selon le cas).

## Validation et reproduction

La vérification complète a réussi en **71 s**, avec les **12/12 cas statiques
convergents**. **`ci/local.sh --bancs` est réussi**, y compris les bancs
étendus et la suite dédiée aux contacts (Hertz, restitution, non lisse,
friction et CCD sur primitives/maillages). Les nouveaux contrôles sont
intégrés à `test_noyau`, chargé par la vérification complète.

Après les derniers ajouts de tests et de validation des tables Rust :
**20/20 tests Rust**, **24/24 tests Python**, Clippy sans avertissement,
formatage conforme et documentation d'API synchronisée (178 entrées).
Sur la v0.5.1 reconstruite, cette même suite Python donne **10 échecs
d'assertion et une erreur d'intégration**, contre zéro sur la v0.5.2.
La CI étendue a précédé les deux derniers tests Python et la garde des tables
Rust ; ces ajouts ont ensuite été contrôlés par les suites ciblées et la
reconstruction du paquet. Aucun gain de performance n'est déduit des temps
des suites de validation.

```bash
cargo fmt --check
cargo clippy --release --offline --all-targets -- -D warnings
cargo test --release --offline
.venv/bin/maturin develop --uv --release --offline
.venv/bin/python -m unittest vinkulum.test_noyau -v
CARGO_NET_OFFLINE=true ci/local.sh --bancs
RAYON_NUM_THREADS=1 .venv/bin/python ci/audit_noyau.py chaine900
RAYON_NUM_THREADS=4 .venv/bin/python ci/audit_noyau.py chaine900
```

## Limites explicites

- Le contrôle géométrique statique est indépendant des charges, mais garde
  les unités de chaque ligne de Φ (longueurs et radians). Ce lot ne construit
  pas un système universel de tolérances adimensionnelles. Le repli historique
  de stagnation peut accepter un résidu de force relatif jusqu'à `1e-6` ; il
  ne relâche désormais plus le contrôle géométrique.
- Pour les engrenages non entiers, les angles restent déroulés : l'erreur
  d'arrondi peut croître sur de très longues simulations. Le test de 20 s
  ne démontre pas une précision indépendante du nombre de tours.
- `modes_complexes` publie les paires oscillantes ; les racines purement
  réelles sont exclues. Une liste sans σ positif ne démontre donc pas la
  stabilité complète. Cette limite est désormais indiquée dans l'API.
- La détection préalable des redondances reste limitée à 900 contraintes.
  Les analyses modales restent denses. Ce lot ne fournit ni QR creuse de rang
  ni solveur spectral partiel pour des millions de degrés de liberté.
- Les structures Rust publiques permettent toujours de contourner les
  constructeurs validants en modifiant directement leurs champs. Il ne
  s'agit pas d'une garantie d'absence de panique pour toute structure Rust
  arbitrairement mal formée.
- La concurrence est contrôlée sur des instances indépendantes avec le GIL
  libéré dans les intégrateurs. Aucun interpréteur Python 3.14t n'est utilisé.
  Aucun profileur d'allocations natives n'est disponible : les suppressions
  d'allocations sont établies par le code, le RSS est mesuré séparément.
- Les garanties de CCD et les domaines des intégrateurs restent ceux déjà
  documentés. Les bancs de contact ne constituent pas une preuve exhaustive
  pour toutes les géométries et trajectoires possibles.
