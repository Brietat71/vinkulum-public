# Statique : les réactions ne doivent pas masquer le déséquilibre libre

Étape archivée. Le [rapport et le mode strict](STATUT_STATIQUE.md) traitent
ensuite la distinction entre tolérance atteinte et arrêt sur stagnation.

Le critère statique est maintenant normalisé par le résidu initial hors
contraintes, et non par les forces totales avant calcul des réactions.
Ce correctif concerne les deux formulations de poutre et les autres
éléments du noyau ; il ne modifie ni leur énergie ni leurs lois physiques.

## Défaut reproduit

Un corps guidé selon z est relié à un support par un ressort de 10 N/m.
Une force libre de 1 N impose une flèche exacte de 0,1 m. Ajouter 1e12 N
sur l'axe x bloqué ne doit pas changer cet équilibre.

| Charge sur l'axe bloqué | Ancien solveur | Correctif |
|---|---|---|
| 0 N | 0,1 m en une itération | 0,1 m en une itération |
| 1e12 N | succès sans itération, déplacement nul, résidu libre 1 N | 0,1 m en une itération, résidu libre 2,22e-16 N |

Le test couvre aussi un départ à 0,09 m. Il échoue avec le binaire
précédent et passe avec le correctif. La référence est F/k ; elle ne
provient pas d'un autre solveur.

## Critère

Le calcul des réactions par moindres carrés fournit déjà λ, puis le
résidu libre r = f − Gᵀλ. Au premier passage de chaque tentative Newton :

```
avant : échelle = max(1, ||f_initial||∞)
après : échelle = max(1, ||r_initial||∞)
```

Le test de convergence reste `max(||r||∞/échelle, ||Φ||∞) <= tol`.
Le même dénominateur sert au mérite en forces et au contrôle de stagnation.
La trace `VINKULUM_TRACE=1` publie désormais l'échelle libre.

**La portée du correctif est précise :** une force éliminée par les
contraintes ne peut plus, par sa seule amplitude, relâcher la tolérance
sur les directions libres. Le critère reste relatif au déséquilibre libre
initial, avec un plancher de 1 ; il n'est donc pas indépendant de toute
initialisation, ni une nouvelle tolérance absolue. La projection elle-même
reste soumise aux limites de l'arithmétique flottante.

Le mécanisme de stagnation existant peut encore accepter après dix
itérations un résidu sous 1e-6 fois cette échelle, avec contraintes sous
`tol`. Il n'a pas été supprimé ni présenté comme une garantie de tenir
`tol` dans tous les cas. Clarifier ce statut et le contrôle d'erreur absolue
reste nécessaire pour renforcer le contrat de convergence.

## Retour au câble

Le [diagnostic précédent](TANGENTE_MIXTE.md) avait observé un faux succès
depuis l'équilibre linéaire calculé séparément : les forces et moments
totaux atteignaient 1,24e8, tandis que le résidu libre restait autour de
4,12e-6 N.

Avec le correctif, l'échelle de cette tentative vaut **1**. Le solveur
refuse la convergence dans le budget de huit itérations et un palier.
L'échec restaure l'état de départ, ce qui explique que le résidu archivé
après l'appel diffère légèrement de celui de la dernière itération dans
le message d'erreur. Depuis l'état original, le cas échoue également.

Ce refus corrige le verdict ; il ne résout pas le modèle extrême. Sa flèche
linéaire de 320 km et son conditionnement de 8,72e10 restent incompatibles
avec l'utilisation de ce cas comme validation d'une caténaire physique.
L'option de poutre intégrée reste expérimentale.

## Validation et reproduction

Formatage et Clippy passent ; **26 tests Rust**, **48 tests Python**,
**41/41 groupes de vérification**, **46/46 bancs rapides** et **9/9 bancs de
contact** passent. La référence d'API reste à jour (183 entrées). L'extension
release est chargée dans le venv Python 3.14 ; l'installation par roue
n'est pas incluse. Le nouveau test vérifie aussi la restauration exacte
après un refus dû à un budget insuffisant.

[Bilan, empreintes et journaux](bancs/echelle-statique-libre-bilan.json).

## Incidence sur Princeton

À charge pleine, un échauffement et trois répétitions, un fil demandé,
hors construction et démarrage Python :

| Formulation | Intervalles | Avant | Après |
|---|---|---|---|
| milieu | 10 | 0.0198 s | 0.0200 s |
| milieu | 20 | 0.0498 s | 0.0501 s |
| milieu | 40 | 0.1126 s | 0.1147 s |
| milieu | 60 | 0.1910 s | 0.1947 s |
| integree | 10 | 0.0201 s | 0.0203 s |
| integree | 20 | 0.0487 s | 0.0491 s |
| integree | 40 | 0.1138 s | 0.1156 s |
| integree | 60 | 0.1978 s | 0.2012 s |

Les positions finales de ces essais restent inchangées. Les deux chemins
en cinquante paliers convergent également. Pour la formulation intégrée à
quarante intervalles, l'essai unique en paliers passe de 2,57 à 3,50 s :
exclure les réactions resserre le critère et peut demander davantage de
travail. Ce chemin n'est pas une médiane de performance, et ce coût n'est
pas masqué dans le bilan.

Mesures : [milieu](bancs/princeton-echelle-libre-milieu.json) et
[intégrée](bancs/princeton-echelle-libre-integree.json).

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv314/bin/python -m unittest \
  vinkulum.test_noyau.AuditNoyau.test_statique_reaction_ne_masque_pas_equilibre
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_cable.py --sortie /tmp/cable-echelle-libre.json
```

[Diagnostic du câble et empreintes](bancs/cable-echelle-libre.json).
L'[objectif généraliste](OBJECTIF_MBDYN.md) reste ouvert.
