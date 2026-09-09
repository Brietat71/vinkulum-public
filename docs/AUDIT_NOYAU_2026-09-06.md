# Audit approfondi du noyau — 6 septembre 2026

Base examinée : `440b58c` (v0.5.0), Rust et interface Python. Les modèles
aérodynamiques et leur confrontation aux essais ne font pas partie de cet
audit. Les contrôles historiques passaient avant correction : 13 tests Rust,
Clippy, formatage, vérification Python en 76 s et campagne statique 12/12.

Les corrections sont livrées dans **v0.5.1**. Cette version corrective conserve
les signatures publiques et les formats de sortie, avec les refus explicites
de configurations auparavant ignorées ou numériquement invalides décrits
ci-dessous. Les références à v0.5.0 dans les mesures désignent la base examinée.

## Défauts et corrections

| Déclencheur | Défaut reproduit | Correction |
|---|---|---|
| Deux masses de 1 kg, K = 100 N/m, allongement initial 0,1 m | Allongement 0,099005 m à 0,1 s au lieu de 0,015594 m | Les forces des superéléments entrent désormais dans le résidu dynamique et sa tangente |
| Superélément en rotation rigide avec beta non nul | Forces d'amortissement sans déformation | Vitesse de déformation objective, transport et réactions sur le premier nœud |
| Effort NaN ou matrice de superélément NaN | Statique annoncée convergée, parfois avec toutes les poses NaN | Refus à la déclaration ; contrôles de finitude avant normes, factorisations et critères de succès |
| Corps libre, J = diag(1,2,3), omega = (1,2,3), h = 3 s | Énergie multipliée par 16,47, sans erreur | Newton énergie-moment conserve le meilleur itéré et refuse une stagnation hors du plancher d'arrondi |
| Pendule avec newton_max = 1 et pas_contact actif | Rejeu sans fin : le pas réduit est rétabli au début de chaque tentative | Le pilotage du contact propose un nouveau pas après progression temporelle ; les rejeux sont bornés |
| Échec de la partition lente du multi-rythme | Corps toujours gelé lors d'une simulation classique ultérieure | Restauration du début du macro-pas, y compris les masques et états du schéma |
| Superélément déformé passé à simule_em | Effort ignoré, masses immobiles | Refus de domaine pour les superéléments et l'appariement automatique |
| Rotation relative exacte de pi | Couple de torsion nul ; panique dans les modes d'une poutre | Logarithme commun, axe extrait de la partie symétrique ; refus des tangentes sur la coupure |
| Modification de la loi d'appariement sur des paires existantes | Ancienne raideur conservée | Mise à jour immédiate des coefficients et invalidation de la symbolique |
| Résolution sous jacobien tourné avec repli SVD | Solution exprimée dans le repère de factorisation | Retour au repère courant aussi après le repli SVD, contrôlé sur un système singulier |

Les lois de couple tronquées, directions invalides, pas/tolérances non finis,
bornes adaptatives invalides et échantillonnage nul sont également refusés.
Les dernières SVD sans plafond de l'assemblage et de la statique ont été
remplacées ; un échec de résolution des réactions statiques est propagé.

## Formulation des superéléments

Le déplacement relatif `u(q)` est exprimé dans le repère du premier nœud.
Avec `B = du/dq`, les forces sont `-Bᵀ K u`, gradient de `uᵀ K u / 2`.
Les réactions du premier nœud incluent le bras de levier des forces nodales.
L'amortissement est `-beta Bᵀ K B v` : une vitesse rigide donne `Bv = 0`.
La dérivée spatiale de ce champ, y compris le transport et l'amortissement,
est calculée par nombres duaux ; les colonnes GGL reçoivent la même
contribution de raideur.

La matrice K conserve l'hypothèse constitutive du superélément initial :
déformations relatives petites, repère suivant le premier nœud. Cette
correction ne transforme pas une réduction linéaire en modèle de matériau
non linéaire. K doit être finie et symétrique ; les très petites asymétries
acceptées dans la tolérance de déclaration sont symétrisées. Aucune positivité
de K n'est imposée : une raideur tangente instable reste représentable.

`alpha != 0` est désormais refusé, car ce paramètre n'avait aucun effet.
`beta` doit être fini et positif ou nul. Les nœuds d'un superélément doivent
être distincts. L'énergie retournée par `energie()` garde sa définition
cinétique plus gravitationnelle ; elle ne comprend pas les potentiels
élastiques. Les projections d'énergie et `simule_em` refusent donc ces modèles.

Un test différentiel du jacobien ne suffit pas : auparavant, résidu et
jacobien omettaient tous deux le superélément et leur comparaison était exacte.
Le nouvel oscillateur est confronté à une solution indépendante, amortissement
compris, et son erreur doit décroître à l'ordre 2 lorsque le pas est divisé.

## Échecs, compatibilité et performances

Les signatures Python et les formats des trajectoires sont conservés.
Les arguments invalides produisent `ValueError`, les indices invalides
`IndexError`, et les échecs numériques des intégrateurs/statique `RuntimeError`.
Un appel refusé à la validation ne modifie pas l'état physique.

En cas d'échec d'intégration, le noyau restaure le dernier pas accepté ; en
multi-rythme, le dernier macro-pas accepté. Les positions, orientations,
vitesses, multiplicateurs, gels, angles déroulés, états de fluide, contacts,
moyennes et longueurs d'historique sont restaurés ; les caches symboliques
sont invalidés. Les compteurs de coût conservent le travail effectué, y
compris les tentatives rejetées. `audit_jacobien` restaure aussi l'état.

Le critère énergie-moment descend jusqu'au plancher d'arrondi et conserve
le meilleur itéré. Une remontée importante est un échec, sans accepter le
pas. Pour le Newton dynamique de contact non lisse, le repli de recherche
non monotone préexistant est conservé : interdire toute remontée échouait sur
le cas de contact du corpus des jacobiens. Un essai non fini n'est jamais
accepté. Les seuils d'arrondi historiques du Newton dynamique restent présents.

Optimisations : suppression de l'allocation faer inutile en résolution dense,
résolution creuse directement dans le vecteur résultat, espace de travail
creux conservé dans la factorisation, options de diagnostic lues hors des
descentes, et GIL libéré pour les trois intégrateurs. La sauvegarde des corps
ne copie que les poses et vitesses ; celle des liaisons ne copie que les angles
déroulés. Les grandes géométries et matrices de superéléments ne sont pas
recopiées à chaque pas.

### Mesures comparatives

Deux passes dans des ordres inversés (avant/après puis après/avant), chacune
avec un échauffement et sept mesures, `RAYON_NUM_THREADS=1`, processus séparés.
Les données proviennent de l'extension initiale conservée avant reconstruction
et de l'extension corrigée, sous le même Python 3.14.7. Les mesures sont
effectuées après la fin des autres bancs et compilations lancés pour cet audit.
Le [relevé brut](bancs/audit-noyau-2026-09-06.json) conserve les temps et les
contrôles physiques de chaque répétition.

| Cas | Médiane avant, ms/pas | Médiane après, ms/pas | Variation |
|---|---:|---:|---:|
| Pendule | 0,014173 | 0,014200 | +0,2 % |
| Chaîne, 30 corps | 0,198354 | 0,202668 | +2,2 % |
| Poutre, 8 éléments | 0,114713 | 0,117453 | +2,4 % |
| Chaîne creuse, 900 corps | 9,958087 | 10,584551 | +6,3 % |

Ces mesures ne démontrent **aucune accélération globale**. Les petites
durées varient fortement selon l'ordre : le pendule initial a donné de
0,009238 à 0,022618 ms/pas. Les nouvelles garanties de restauration et de
finitude ont un coût ; les réductions d'allocations ne le compensent pas
sur tous ces cas. Le pic RSS cumulé maximum de ces processus est passé de
119 072 à 119 704 Kio, sans baisse démontrée de la mémoire résidente.
L'amélioration validée de ce lot est d'abord la fiabilité physique et la
gestion des erreurs. Les optimisations de stockage conservent les résultats
et réduisent les allocations temporaires visibles dans le code ; leur gain
de temps isolé n'a pas été établi par ces mesures du lot complet.

## Validation réalisée

Sur Python 3.14.7 et le profil Rust release :

- 17 tests Rust réussis, dont les contrôles différentiels et mécaniques des
  superéléments ; formatage conforme et Clippy sans avertissement.
- 12 tests Python de régression réussis, intégrés à la vérification complète.
- Vérification complète réussie en 73 s, avec 12/12 cas statiques convergents.
- Bancs étendus `bancs rapide` réussis : contact, CCD, adaptation du pas, GGL,
  énergie-moment, multi-rythme et cas aérodynamiques existants.
- Référence d'API vérifiée à jour : 178 entrées.

Cette validation n'inclut pas une exécution avec Python sans GIL (`3.14t`).

## Reproduction

```bash
cargo fmt --check
cargo clippy --release --offline --all-targets -- -D warnings
cargo test --release --offline
.venv/bin/maturin develop --uv --release --offline
.venv/bin/python -m unittest vinkulum.test_noyau -v
.venv/bin/python -m vinkulum.verification
.venv/bin/python -m vinkulum.bancs rapide
.venv/bin/python -m vinkulum.doc
RAYON_NUM_THREADS=1 .venv/bin/python ci/mesure_noyau.py
```

La vérification complète inclut les tests de l'audit. Les cas susceptibles de
boucler sont exécutés dans un sous-processus avec délai maximal et contrôle
du code de sortie. Les tests Rust contrôlent aussi le travail virtuel, les
résultantes internes, les tangentes en pose/vitesse sur deux perturbations,
le signe de la puissance d'amortissement, le logarithme à pi et le repli SVD.
