# Poutre intégrée : précision améliorée, option expérimentale

Étape archivée. Le [diagnostic suivant](TANGENTE_MIXTE.md) précise le cas
câble et identifie aussi une dépendance du critère statique à son échelle
de départ ; la limite décrite ici n'est pas considérée comme résolue.

Une formulation supplémentaire est accessible par
`n.poutre(..., formulation="integree")`. La formulation `"milieu"` reste le
défaut. L'objectif est de réduire le maillage nécessaire à une précision
donnée ; le remplacement général de l'élément historique n'est pas validé.

À dix intervalles sur Princeton, le prototype réduit l'écart à la référence
MBDyn raffinée de **254,87 µm à 5,33 µm**, soit environ **48 fois moins
d'erreur**. Une console d'un mètre sous moment pur retrouve aussi l'arc
circulaire avec un seul élément. Ces résultats ne classent pas encore les
deux moteurs en vitesse à précision commune.

## Pourquoi changer la formulation

La mesure de raideur au repos retrouve le défaut connu de l'élément à
intégration au milieu : en flexion élancée, un élément donne le facteur
1/4 dans la flèche sous force terminale, au lieu de 1/3. Ce comportement
est décrit dans le [cours de Timoshenko de TU Delft](https://teachbooks.tudelft.nl/computational-modelling/structural_linear/timoshenko.html).
Notre diagnostic mesure une sous-estimation de 24,93 % pour L = 1 m,
EI = 100 N·m² et GA = 100000 N, incluant la déformation de cisaillement.

Le second défaut est géométrique : prendre la corde entre nœuds pour une
longueur d'arc produit une erreur même sous moment constant, lorsque la
courbure exacte est constante.

## Construction de l'option

Les deux corrections ci-dessous sont combinées dans une énergie objective.
Cette combinaison est une formulation développée ici ; elle ne doit pas
être présentée comme un élément non linéaire exact pour tout chargement.

Pour la rotation relative matérielle θ et la rotation moyenne R_m, on
intègre R(s) = R_m exp((s/L − 1/2)[θ]×). L'opérateur J_sym de cette intégrale
agit comme l'identité sur l'axe θ et comme sinc(|θ|/2) dans le plan normal.
On calcule donc :

```
γ = J_sym(θ)⁻¹ R_mᵀ (r_B − r_A)/L − e₁
κ = θ/L
```

L'inverse utilise une série régulière près de zéro, y compris pour les
variations calculées par duaux. La coupure du logarithme à π garde les
limites et refus existants.

Les flexibilités transverses condensées sont :

```
1/C_y = 1/GA_y + L²/(12 EI_z)
1/C_z = 1/GA_z + L²/(12 EI_y)
U = L/2 (EA γ_x² + C_y γ_y² + C_z γ_z² + κᵀ C_M κ)
```

La limite linéaire retrouve la matrice nodale de Timoshenko, avec ses deux
rigidités anisotropes. Les paramètres GA et EI fournis restent les
paramètres physiques ; C_y et C_z sont des coefficients discrets. Pour
cette option, les γ retournés par `poutres()` sont des déformations
**généralisées** : les multiplier directement par les GA d'entrée ne donne
pas les efforts transverses de l'élément condensé.

Les sensibilités appliquent la règle de chaîne à ces coefficients avant
la dérivation par duaux. L'ancienne hypothèse de linéarité de l'énergie en
chaque paramètre ne suffit plus.

## Contrôles

- Matrice 12×12 au repos comparée à la matrice analytique anisotrope.
- Arc sous moment pur, rotations finales 0,01, 1 et 2,5 rad, puis même
  modèle après rotation et translation dans le monde.
- Dérivées des forces et des raideurs pour EA, GA, GJ, EI et chaque EI
  transverse, comparées à des différences finies indépendantes à pose fixée.
- Princeton à 10, 20 et 40 intervalles avec les deux formulations : équilibre,
  contraintes et convergence spatiale vérifiés.

Formatage et Clippy passent, ainsi que **25 tests Rust**, **47 tests Python**,
**41/41 groupes de vérification**, **46/46 bancs rapides** et **9/9 bancs de
contact**. La référence d'API est à jour (183 entrées). La campagne générale
utilise le défaut historique ; les nouveaux tests exercent explicitement
l'option intégrée. L'adjoint dynamique est aussi rejoué avec cette option
sur un pendule flexible de 300 pas et confronté à trois différences finies.
L'extension release est chargée dans le venv Python 3.14 ; l'installation
par roue n'est pas incluse dans ces contrôles.

## Mesures définitives

Même binaire, un fil demandé, un échauffement puis trois répétitions par
maillage. Le temps inclut le chargement et la résolution statique, hors
construction du modèle et démarrage Python. L'écart est celui de la position
du bout à la référence MBDyn à 160 intervalles ; il ne borne pas l'erreur
sur tous les déplacements, contraintes ou efforts intérieurs.

| Intervalles | Erreur milieu | Erreur intégrée | Temps milieu | Temps intégrée |
|---|---|---|---|---|
| 10 | 254.87 µm | 5.330 µm | 0.0262 s | 0.0266 s |
| 20 | 63.72 µm | 1.346 µm | 0.0630 s | 0.0622 s |
| 40 | 15.93 µm | 0.337 µm | 0.1394 s | 0.1436 s |
| 60 | 7.08 µm | 0.150 µm | 0.2322 s | 0.2419 s |

Sous le seuil de 7,1 µm au bout, les configurations mesurées à 60 intervalles
pour le milieu et 10 pour l'intégrée donnent **0,2322 s contre 0,0266 s**,
soit **8,73×**. À dix intervalles, l'erreur est réduite de **47,82×**.
À maillage égal, l'option n'accélère pas systématiquement la résolution.
Les deux chemins en cinquante paliers convergent aussi à 10/20/40
intervalles ; leurs temps uniques sont des contrôles, pas des médianes.

- [Mesures milieu](bancs/princeton-poutre-milieu.json) et
  [intégrées](bancs/princeton-poutre-integree.json).
- [Diagnostic analytique et câble, milieu](bancs/poutre-diagnostic-milieu.json)
  et [intégré](bancs/poutre-diagnostic-integree.json).
- [Gains recalculés, empreintes et journaux](bancs/poutre-integree-bilan.json).

## Limite qui interdit le remplacement par défaut

Le corpus interne contient une chaîne dont les rotations nodales sont
bloquées et EI vaut 1e-6 N·m². La flexibilité condensée change fortement ce
modèle : à cette souplesse, la déformation interne représentée par la
condensation linéaire peut devenir grande. La campagne avec remplacement
global a été interrompue pendant ce cas ; les autres groupes terminés ne
valent pas validation complète de cette substitution.

Un diagnostic séparé reprend exactement la chaîne de vingt éléments et
borne chaque processus à vingt secondes. La formulation historique termine ; l'option intégrée dépasse ce délai.
Les deux résultats sont conservés, sans exclusion du cas. Il s'agit d'une limite de robustesse
observée ; sa cause complète et la validité de la condensation dans ce
régime restent à établir. Il faut traiter la flexion interne non linéaire
avant de prétendre que l'option remplace l'élément sur tous les domaines.

## Reproduction

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_poutres.py \
  --formulation integree --sortie /tmp/poutre-integree.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_statique.py \
  --formulation integree --sortie /tmp/princeton-integre.json
```

Pour rejouer le contrôle de l'adjoint dynamique :

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_poutres.py \
  --formulation integree --adjoint-worker
```

Passer `--formulation milieu` pour le contrôle historique. Les résultats,
paramètres, empreintes du binaire et sources sont enregistrés dans les JSON.
L'[objectif généraliste face à MBDyn](OBJECTIF_MBDYN.md) reste ouvert.
