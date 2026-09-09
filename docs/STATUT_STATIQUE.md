# Statique : tolérance, stagnation et mode strict

`statique()` expose maintenant la différence entre une tolérance atteinte
et une acceptation de secours. Le nouveau `statique_info()` donne le bilan
du dernier calcul, et `strict=True` refuse l'acceptation sur stagnation.
Le couple de retour historique `(résidu, compteur)` reste disponible.

## Comportement

| Résultat | Mode habituel (`strict=False`) | Mode strict |
|---|---|---|
| Tolérance atteinte à tous les paliers retenus | Retour normal, statut `tolerance` | Identique |
| Arrêt accepté sur stagnation | `RuntimeWarning`, statut `stagnation` | Reprises et continuation possibles, puis refus si la tolérance reste hors d'atteinte |
| Échec du calcul | Exception et restauration de l'état | Identique |

L'avertissement est émis une fois par appel concerné, selon les filtres
Python. Un filtre qui le transforme en exception provoque aussi une
restauration. Le test dédié vérifie ce cas après un déplacement effectif
produit par la détente d'un ressort.

```python
residu, compteur = n.statique(tol=1e-8, strict=True)
bilan = n.statique_info()
assert bilan["statut"] == "tolerance"
```

Le défaut reste `strict=False` : les modèles qui utilisent le secours
numérique continuent à fonctionner, mais leur arrêt n'est plus silencieux.
Cela ne démontre pas qu'ils atteignent la tolérance demandée.

## Contenu du bilan

Le bilan est un instantané historique indépendant des modifications
ultérieures du modèle. Le dictionnaire retourné est une copie. Il vaut
`None` si aucun rapport n'est disponible. Les validations internes du
nouvel appel le réinitialisent ; une erreur de conversion d'argument dans
la liaison Python, avant l'entrée dans le noyau, ne constitue pas un calcul.

- `statut` : `tolerance`, `stagnation` ou `echec`.
- `tol`, `strict` : réglages du calcul.
- `residu_libre`, `contraintes`, `echelle_force`, `residu_relatif` : dernière
  évaluation Newton valide. Le résidu mécanique inclut les composantes de
  force et de moment dans les coordonnées du noyau.
- `tolerance_finale_atteinte` : test de la tolérance à cette évaluation.
  Il peut être vrai alors qu'un palier antérieur a utilisé le secours.
- `evaluations`, `tentatives` : nombres d'évaluations principales du résidu
  et de tentatives Newton, y compris celles abandonnées. Les évaluations de
  recherche amortie ne sont pas incluses dans le premier compteur.
- `paliers`, `palier`, `paliers_stagnation` : dernière subdivision tentée,
  dernier palier tenté et nombre de paliers acceptés sur stagnation dans
  cette subdivision. Une subdivision abandonnée ne pollue pas ce dernier
  compteur.
- `etat_restaure`, `message` : restauration après échec et motif éventuel.

Après un échec, les résidus décrivent la tentative avant restauration,
**pas l'état courant restauré**. Une tentative interrompue avant sa première
évaluation valide donne des champs de résidu à `None`.

Le mode strict utilise le critère existant : résidu libre relatif à son
échelle initiale, avec plancher 1, et contraintes sous `tol`. Il ne promet
ni une erreur absolue sur les déplacements ni une invariance à toutes les
unités. Il désactive précisément le secours de stagnation à 1e-6 de
l'échelle ; il n'assouplit aucun autre contrôle.

## Ce que les essais révèlent

Un corps libre soumis à une force constante de 5e-7 N n'a aucun équilibre.
L'ancien binaire rendait pourtant `(5e-7, 10)` sans avertissement, avec
`tol=1e-12`. Le mode habituel rend désormais le statut `stagnation` et un
avertissement ; le mode strict refuse et restaure l'état.

Sur Princeton, à `tol=1e-8`, avec 100 itérations et un seul palier interne
autorisés par appel :

Ce tableau archive le noyau avant les
[translations compensées](POSITIONS_COMPENSEES.md), qui font désormais
passer ces parcours stricts.

| Maillage | Chargement | Habituel | Strict |
|---|---|---|---|
| 10 intervalles | Charge pleine | Tolérance atteinte | Tolérance atteinte |
| 10 intervalles | 50 paliers externes | Tolérance atteinte aux 50 paliers | Identique |
| 40 intervalles | Charge pleine | Tolérance atteinte | Tolérance atteinte |
| 40 intervalles | 50 paliers externes | Calcul terminé, 49 arrêts sur stagnation | Refus au premier palier |

Au premier palier strict à quarante intervalles (0,0987 % de la charge),
le dernier résidu relatif vaut 1,90e-8, au-dessus de 1e-8, après 111
évaluations principales sur deux tentatives. L'état de départ est restauré.

Les deux parcours ne prouvent donc pas la même chose, même si leurs
positions finales paraissent proches. Le diagnostic conserve chaque
avertissement, les rapports, les paliers refusés et la vérification de
restauration. Un calcul strict refusé n'est pas comparé à la référence de
la charge pleine comme s'il l'avait atteinte.

## Validation et reproduction

Formatage et Clippy passent ; **26 tests Rust**, **50 tests Python**,
**41/41 groupes de vérification**, **46/46 bancs rapides** et **9/9 bancs de
contact** passent. La référence d'API est à jour (184 entrées). L'extension
release est chargée dans le venv Python 3.14 ; l'installation par roue
n'est pas rejouée.

Les nouveaux tests couvrent le rapport et son cycle de vie, une convergence
stricte sur l'équilibre analytique du ressort, le cas sans équilibre,
l'avertissement transformé en exception et la restauration après un
mouvement effectif. Les campagnes générales conservent le mode habituel ;
elles ne prouvent pas que tous leurs modèles passent en mode strict.

[Bilan, témoin ancien, empreintes et journaux](bancs/statut-statique-bilan.json).
Le [source antérieur du module principal](bancs/statut-statique-ancien-lib.rs.gz)
est conservé pour ce témoin ; son empreinte est dans le bilan.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_statut_statique.py \
  --sortie /tmp/statut-statique.json
```

Le programme de mesure `ci/mesure_statique.py` enregistre également les
rapports de chaque appel et accepte `--strict`. Il conserve les mesures
déjà écrites si un calcul strict ultérieur échoue.

[Diagnostic détaillé](bancs/statut-statique-diagnostic.json).
Le [diagnostic ultérieur de précision](PRECISION_STATIQUE.md) isole
l'effet de l'arrondi des positions sur le premier palier. Sa reconstruction
porte uniquement sur les translations à rotations fixées ; elle ne fait
pas encore converger le solveur strict complet.
L'[objectif généraliste face à MBDyn](OBJECTIF_MBDYN.md) reste ouvert : ce
changement rend les verdicts contrôlables, mais ne résout pas les cas
stricts qui échouent ni la limite physique de la poutre intégrée.
