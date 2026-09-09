# Tangente de poutre : blocs de translation et couplages exacts

Étape archivée. Le [correctif suivant](ECHELLE_STATIQUE_LIBRE.md) exclut les
réactions de la normalisation et supprime le faux succès décrit ci-dessous.
Il ne déclare pas le câble extrême résolu.

La tangente utilise maintenant la force de translation analytique et sa
dérivée par duaux. Les blocs qui couplent translation et rotation sont
obtenus par transposition ; le bloc rotation/rotation garde ses différences
finies centrées. Le changement s'applique aux deux formulations de poutre,
sans modifier leur énergie ni les forces physiques.

## Calcul

Avec B = J_sym⁻¹ R_mᵀ pour l'option intégrée (B = R_mᵀ pour l'élément
historique), γ = B(r_B − r_A)/L − e₁. La force sur A est Bᵀ C γ et celle
sur B est son opposée. Différencier cette expression fournit les six
lignes de translation en un passage de duaux.

Les perturbations de translation dans le monde et de rotation gauche
commutent : les blocs mixtes sont transposés. Cela ne suppose pas que le
bloc rotation/rotation soit symétrique hors équilibre. Six colonnes
rotatives restent calculées par différences finies au pas existant ; les
six perturbations translatoires et leurs soustractions de forces sont
supprimées.

Le nouveau test compare tous ces blocs à la dérivation de l'énergie par
duaux emboîtés, pour les deux formulations, avec rotations relatives et
translations ordinaires puis amplifiées d'un million. Il contrôle une
erreur par entrée sous 1e-10*(1+|référence|).

## Mesures et validation

Formatage et Clippy passent ; **26 tests Rust**, **47 tests Python**,
**41/41 groupes de vérification**, **46/46 bancs rapides** et **9/9 bancs de
contact** passent. L'API est à jour (183 entrées). L'extension release est
chargée dans le venv Python 3.14 ; l'installation par roue n'est pas rejouée.

Les mesures de Princeton utilisent un fil,
un échauffement puis trois répétitions, sans autre campagne concurrente,
avec le même protocole que l'étape précédente. Le gain porte sur le calcul
statique ; la précision spatiale ne change pas. Les positions finales
varient de moins de 2e-13 m sur ces essais.

| Formulation | Intervalles | Avant | Après | Temps économisé |
|---|---|---|---|---|
| milieu | 10 | 0.0262 s | 0.0198 s | 24.4 % |
| milieu | 20 | 0.0630 s | 0.0498 s | 20.9 % |
| milieu | 40 | 0.1394 s | 0.1126 s | 19.2 % |
| milieu | 60 | 0.2322 s | 0.1910 s | 17.8 % |
| integree | 10 | 0.0266 s | 0.0201 s | 24.4 % |
| integree | 20 | 0.0622 s | 0.0487 s | 21.7 % |
| integree | 40 | 0.1436 s | 0.1138 s | 20.8 % |
| integree | 60 | 0.2419 s | 0.1978 s | 18.3 % |

[Bilan calculé, empreintes et contrôles](bancs/tangente-mixte-bilan.json).
Mesures statiques : [milieu](bancs/princeton-mixte-milieu.json) et
[intégrée](bancs/princeton-mixte-integree.json). Matrices d'analyse :
[avant](bancs/tangente-mixte-avant.json) et [après](bancs/tangente-mixte-apres.json).

Les matrices d'analyse `k_c_m_z()` n'accélèrent pas autant : à 120
intervalles, leurs médianes sont 89,23 et 88,99 ms. L'assemblage global et
la matérialisation de matrices denses limitent le gain sur cette API.

## Le diagnostic du câble sépare deux problèmes

Le modèle du corpus bloque les rotations de tous les nœuds. Ses efforts
internes deviennent donc linéaires en translation : on peut assembler
indépendamment R C Rᵀ/L et résoudre uniquement les translations libres.
Le [diagnostic reproductible](../ci/diagnostic_cable.py) conserve le modèle
original et limite les tentatives Newton à huit itérations et un palier.

| Formulation | Conditionnement du bloc libre | Flèche de l'équilibre linéaire |
|---|---|---|
| Milieu | 162,5 | 0,50371 m |
| Intégrée | 8,72e10 | 320350,4 m |

Le bloc du noyau diffère de l'assemblage analytique de moins de 1,5e-16 en
relatif. La flèche immense est une propriété de ce modèle condensé avec
rotations bloquées, pas une preuve de comportement physique correct d'un
câble. Le modèle historique ne reproduisait pas non plus la caténaire.

Depuis l'état original, l'option intégrée échoue encore dans le budget
borné. Depuis l'équilibre linéaire, `statique()` annonce un succès avec
un résidu libre d'environ 4,12e-6 N. **Cela ne valide pas la tolérance
initiale** : l'échelle des forces et moments internes utilisée par le
solveur passe de 14,715 à 1,24e8. Rapporté à la première échelle, le résidu
vaut 2,80e-7, au-dessus de la tolérance 1e-10.

Le diagnostic archive donc aussi l'échelle de départ de chaque essai et
la vérification à échelle fixe. Il faut traiter cette dépendance du critère
statique à l'état initial, puis la validité de la condensation dans les
régimes extrêmes. L'option intégrée reste expérimentale et le défaut
historique est conservé.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_cable.py --sortie /tmp/cable-lineaire.json
```

[Résultats, positions de référence et résidus](bancs/cable-tangente-mixte.json).
L'[objectif généraliste](OBJECTIF_MBDYN.md) reste ouvert.
