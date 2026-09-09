# Petits déplacements conservés dans le noyau

Le noyau conserve désormais la partie d'un incrément de translation qui
disparaissait lors de son addition à une position absolue en f64. Les
poutres utilisent cette information dans leurs déformations, leurs forces
et leurs tangentes. Les superéléments et les contraintes de liaison,
distance et vis la prennent également en compte. La géométrie de contact
est évaluée après recentrage des positions des corps.

**Princeton passe maintenant ses cinquante paliers en mode strict** pour
les maillages de 10, 20, 40 et 60 intervalles, dans les deux formulations
de poutre. Aucun secours de stagnation n'est utilisé. Les tolérances,
les budgets et les règles d'acceptation ne sont pas assouplis.

## Preuves d'équilibre

Le [diagnostic précédent](PRECISION_STATIQUE.md) avait seulement reconstruit
les translations à rotations fixées. Le présent résultat vient de Newton
dans le noyau complet, avec les rotations et les contraintes résolues.

Réglages : premier à cinquantième palier de la rampe cosinus,
`tol=1e-8`, `iters=100`, `paliers_max=1`, `strict=True`.

| Formulation | Intervalles | Paliers à la tolérance | Plus grand résidu relatif accepté |
|---|---:|---:|---:|
| milieu | 10 | 50/50 | 5,85e-9 |
| milieu | 20 | 50/50 | 5,13e-9 |
| milieu | 40 | 50/50 | 6,88e-9 |
| milieu | 60 | 50/50 | 7,51e-9 |
| intégrée | 10 | 50/50 | 7,86e-9 |
| intégrée | 20 | 50/50 | 6,39e-9 |
| intégrée | 40 | 50/50 | 8,04e-9 |
| intégrée | 60 | 50/50 | 6,45e-9 |

À quarante intervalles, le noyau antérieur refusait dès le premier palier
strict. Le contrôle indépendant en précision étendue aux **positions
réellement conservées** donne maintenant 1,32e-9 N sur les translations
de ce premier palier, pour la formulation `milieu`.

Une barre analytique de longueur 1, de raideur axiale `2^30` et soumise à
`2^-30` doit s'allonger de `2^-60`. Sa position arrondie reste exactement
1. Le nouveau noyau conserve l'allongement, équilibre la force et passe
`tol=1e-14` en mode strict. Le témoin est vérifié avec les deux poutres
et avec un superélément. **Le binaire antérieur échoue dans les trois cas**,
avec le même programme et le même budget de quarante itérations.

## Représentation et restauration

Une position contient deux vecteurs f64 : `r` et `r_bas`. La somme réelle
des deux représente la translation. Une addition compensée retient le
reste des petits incréments ; les équations sensibles soustraient leur
référence avant d'ajouter la petite contribution. Il ne s'agit pas d'un
passage de tout le solveur en précision étendue.

Cette représentation traverse Newton, les reprises, l'assemblage initial,
les intégrateurs et les interpolations de sous-pas. Les sauvegardes
internes conservent les deux parties. Les dérivées des poutres et des
superéléments portent la même déformation que les forces.

```python
instantane = n.etat_precis()
# ... calculs ou modifications de l'état ...
n.pose_etat(*instantane)
```

`etat_precis()` fournit les six champs historiques suivis de la liste des
parties basses. Avec ce septième argument, `pose_etat()` valide les
rotations sans les réorthonormaliser et normalise les positions en deux
parties. Les entrées invalides sont refusées avant toute mutation.

`etat()`, `pose()` et les trajectoires historiques restent des vues
arrondies. Restaurer les six champs de `etat()` remet la partie basse à
zéro et conserve la réorthonormalisation historique des rotations.
`etat_precis()` porte l'état cinématique et les `v_i` d'inflow ; il ne
remplace pas une sauvegarde des multiplicateurs, des historiques de
contact ou de tous les états aérodynamiques.

## Contrôles hors Princeton

Les nouveaux tests vérifient l'allongement sous le pas représentable,
ses forces et ses tangentes, la conservation des petits incréments puis
leur annulation, et la restitution exacte après un refus statique.
Ils contrôlent aussi les validations atomiques de l'état précis.

Un mouvement libre de `2^-60` est conservé par α-généralisé,
énergie–moment et multi-rythme, avec restauration des orientations non
triviales. Liaison, bielle et vis détectent puis corrigent un déplacement
de `2^-60`. Un contact donne les mêmes efforts près de l'origine et
après une translation de `2^40`, avec support fixe ou mobile.

La vis–écrou vérifie encore l'inertie ramenée et sa contrainte de vitesse.
Son évaluation retire désormais la position de référence avant le petit
terme hélicoïdal et conserve l'erreur de soustraction des positions.
Les gradients des contraintes sont confrontés aux duaux emboîtés et aux
différences finies, avec des parties basses amplifiées pour exposer une
éventuelle omission.

L'adjoint dynamique de la poutre intégrée est également rejoué, y compris
le pendule flexible à 300 pas et les trois différences finies en EI.
Le câble condensé extrême reste en échec dans les deux initialisations du
diagnostic borné ; sa flèche linéaire de 320 km reste une limite du modèle.

## Coût mesuré et validation finale

Temps du calcul statique direct et de la lecture de son bilan, hors imports
et construction du modèle : médiane de trois répétitions après échauffement,
un fil demandé, sans autre campagne Vinkulum en cours. Les deux binaires
utilisent le même programme et le mode strict.

| Formulation | Intervalles | Avant (s) | Après (s) | Variation |
|---|---:|---:|---:|---:|
| milieu | 10 | 0,01997 | 0,02095 | +4,9 % |
| milieu | 20 | 0,05041 | 0,05437 | +7,8 % |
| milieu | 40 | 0,11356 | 0,12121 | +6,7 % |
| milieu | 60 | 0,19419 | 0,20471 | +5,4 % |
| intégrée | 10 | 0,02018 | 0,02195 | +8,8 % |
| intégrée | 20 | 0,04939 | 0,05298 | +7,3 % |
| intégrée | 40 | 0,11513 | 0,12472 | +8,3 % |
| intégrée | 60 | 0,20489 | 0,21223 | +3,6 % |

Le calcul direct coûte donc davantage dans cette série. Les positions
finales et la précision spatiale sont conservées à l'arrondi près. Les
programmes exécutés avec l'ancien binaire terminent leurs mesures directes,
puis échouent sur la rampe stricte à quarante intervalles ; leurs mesures
déjà écrites et leurs sorties d'erreur restent archivées. Un échec n'est
pas classé comme un calcul rapide réussi. Ces temps n'ont pas le périmètre
du processus MBDyn complet de la confrontation initiale.

Formatage et Clippy passent, ainsi que **28 tests Rust**, **55 tests
Python**, **41/41 groupes de vérification**, **46/46 bancs rapides** et
**9/9 bancs de contact**. L'API générée est à jour avec 185 entrées.
Les quatre tests analytiques du diagnostic numérique passent également.
L'extension release est chargée dans le venv Python 3.14 ; l'installation
par roue n'est pas rejouée.

## Reproduction et portée

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_positions_compensees.py \
  --sortie /tmp/positions-compensees.json
```

[Rapports, états précis et contrôle indépendant](bancs/positions-compensees-diagnostic.json).
[Témoin exécuté avec l'ancien binaire](bancs/positions-compensees-ancien.json).
[Sources antérieures archivées](bancs/positions-compensees-sources-avant.tar.gz).
[Diagnostic du câble](bancs/cable-positions-compensees.json).
[Mesures, empreintes et validation](bancs/positions-compensees-bilan.json).

Ces cas établissent un progrès de robustesse à tolérance inchangée. Ils
ne prouvent pas une précision universelle pour toutes les géométries,
toutes les unités ou tous les types d'interaction. La normalisation du
critère statique reste celle décrite dans [son contrat](STATUT_STATIQUE.md).
La précision spatiale des formulations de poutre n'est pas améliorée par
ce changement, et aucune nouvelle supériorité chronométrique face à MBDyn
n'est déduite de cette étape.
