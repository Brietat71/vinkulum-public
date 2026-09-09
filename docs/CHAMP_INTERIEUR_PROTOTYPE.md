# Contrôler le champ physique, puis enrichir les intérieurs

7 septembre 2026. Prototype de recherche sur modèles linéaires conservatifs,
séparé de l'API de simulation. La version publique reste 0.9.0.

La précision du Schur ne suffisait pas : l'erreur du fonctionnel énergétique
est quadratique en résidu, celle du champ linéaire. Le prototype contrôle
désormais déplacements en norme massique, déformations en norme énergétique
et observables linéaires. Il enrichit les bases en réutilisant les facteurs
statiques. Les preuves et leur portée se trouvent dans
[le carnet mathématique](CHAMP_INTERIEUR_PREUVES.md), les antériorités dans
[les cinq sources primaires](CHAMP_INTERIEUR_SOURCES.md).

## Ce qui est démontré

Pour l'intérieur coercif A = D_IᵀD_I − ω²M_II, avec
D_IᵀD_I ≥ λ_* M_II et α = λ_* − ω² > 0, un candidat de résidu R dans
le facteur D d'entrée vérifie, à port imposé :

\[
\|e\|_M\le\frac{\|R\|_{M^{-1}}}{\alpha},\qquad
\|De\|\le\frac{\sqrt{\lambda_*}\,\|R\|_{M^{-1}}}{\alpha}.
\]

Les bornes locales sont uniformes sur la bande de construction grâce à
l'enveloppe résiduelle déjà établie. Après assemblage, on ajoute l'erreur
des ports et la norme du relèvement. La marge globale utilise la plus petite
valeur singulière du Schur réduit, sa perturbation bornée et le résidu de la
résolution réduite. Elle reste nécessaire même si chaque intérieur est
coercif. Les masses M_IS, les demi-masses d'interface et les masses externes
sont incluses dans les normes.

La borne relative vise la solution exacte : si un majorant absolu B est
strictement inférieur à la norme du candidat, le rapport est
B/(norme_candidate−B). Une solution proche de zéro ne reçoit pas un faux
succès relatif. Les observables disposent en plus d'une correction primal-dual
dont le reste est borné par le produit des deux résidus, divisé par α.

Ces énoncés sont exacts sous leurs hypothèses. Le code les évalue en doubles ;
il ne certifie ni les arrondis ni les constantes spectrales fournies.

## Expériences reproductibles

Les cinq expériences sont archivées avec les sources exécutées dans
[champ-interieur-2026](bancs/champ-interieur-2026/manifest.json).
Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, environnement de la roue figée
0.9.0, un fil BLAS. Il s'agit d'expériences de précision, sans classement
chronométré ni nouvelle comparaison avec un solveur externe.

Chaîne discrète : n masses m=1, n ressorts k=1, racine fixée, force physique
terminale 1/√n. Chaque segment est construit indépendamment et les masses
d'interface sont recomposées. Bande locale commune :
Ω = 0,8·2 sin(π/(2n)). La tolérance initiale du Schur assemblé est 1e−10.
La tolérance demandée sur **les deux normes physiques relatives** vaut
1e−10. L'oracle analytique du système discret complet ne participe pas à
l'arrêt. Les essais ont été définis avant la campagne archivée.

Trois fréquences demandées : 0, 0,37Ω et Ω. Les directions indiquées sont
les intérieurs réduits ; ajouter le nombre de ports globaux (4, 8 ou 32).

| n / segments | Directions avant → après | Erreur massique avant → après | Erreur de déformation avant → après | Borne relative finale de déformation |
|---|---:|---:|---:|---:|
| 64 / 4 | 24 → 32 | 1,42e−11 → 2,89e−15 | 5,23e−10 → 1,26e−13 | 1,68e−12 |
| 2048 / 8 | 32 → 64 | 6,27e−10 → 1,47e−14 | 3,70e−8 → 2,89e−13 | 2,81e−11 |
| 2048 / 32 | 128 → 192 | 1,75e−13 → 8,63e−14 | 3,58e−11 → 2,81e−13 | 7,42e−11 |

Ce sont des maxima sur les fréquences demandées. Les trois critères sont
acceptés par leurs bornes. Sur le dernier cas, l'oracle indiquait déjà une
erreur assez petite ; le majorant initial ne permettait pas de l'affirmer.
L'enrichissement conserve les mêmes objets facteurs QR et toutes les
directions précédentes. Il procède par rondes sur tous les intérieurs,
sans marquage sélectif ni affirmation de dimension minimale.

## Les refus près de la résonance restent visibles

Les deux autres expériences demandent 17 fréquences uniformément espacées
de 0 à Ω, donc incluent un point proche du premier mode global
ω₁ = 2 sin(π/(4n+2)). Après trois rondes autorisées :

| n / segments | Directions avant → après | Erreur massique finale | Erreur de déformation finale | Borne relative finale de déformation | Statut |
|---|---:|---:|---:|---:|---|
| 2048 / 8 | 32 → 80 | 2,53e−12 | 2,55e−12 | 3,42e−9 | budget atteint, refus |
| 2048 / 32 | 128 → 320 | 3,55e−10 | 3,55e−10 | 9,25e−8 | budget atteint, refus |

Le premier refus est conservateur : l'oracle passe, la borne échoue. Le
second accompagne une erreur réellement supérieure à la tolérance demandée.
Ajouter des directions intérieures ne supprime pas le défaut de Schur
amplifié par la résonance. Aucun théorème d'impossibilité de mieux calculer
n'est déduit de ces observations ; améliorer les bornes et la précision
de l'assemblage constitue une question distincte.

**Portée de l'acceptation :** les enveloppes locales couvrent la bande ;
le contrôle relatif global est ponctuel, aux fréquences explicitement
demandées. La grille ne certifie pas les fréquences intermédiaires.

## Contrôles et transfert natif

Onze tests ajoutés confrontent les champs à l'assemblage dense complet de
petits systèmes indépendants : masses couplées, énergie, identité primal-dual,
dual exact, refus de marge, critères périmés après enrichissement, budgets,
force nulle et réutilisation QR/LU. Les 14 tests Krylov et 15 tests des ports
relevés précédents passent également. Le vérificateur d'archive conserve les
refus et distingue sa marge de comparaison flottante de toute certification.

L'[extraction matérielle native des poutres](FACTEURS_ENERGIE_NOYAU.md)
fournit en Rust les six déformations pondérées z et leur jacobien D, avec
U = ½‖z‖² et f = −Dᵀz. Elle n'obtient pas D par factorisation de K.
Au repos sans contraintes internes, DᵀD restitue la raideur élastique.
Sous précontrainte, les contributions géométriques doivent être conservées.
Cet accès Rust expérimental ne fournit encore ni API Python ni raccordement
des ports aux trajectoires multicorps.

Reproduction des expériences et contrôles :

```sh
OPENBLAS_NUM_THREADS=1 python ci/test_champ_interieur.py
OPENBLAS_NUM_THREADS=1 python ci/experience_champ_interieur.py --destination /tmp/champ-interieur-reproduction
python ci/experience_champ_interieur.py --verifier docs/bancs/champ-interieur-2026
```

La CI locale exécute les nouveaux tests et vérifie l'archive. Aucune
supériorité sur Exudyn, MBDyn ou Simpack n'est déduite de ces expériences.
