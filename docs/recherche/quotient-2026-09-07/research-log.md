# Du programme fondamental à une expérience réfutable

Le complément porte sur le premier axe du rapport du 7 septembre : quotient
des contraintes générales, complexité structurelle et stabilité du rang.
La première dérivation, les lectures et les contrôles ont été réalisés
directement après la limite de fils rencontrée lors du rapport fondamental.
À la nouvelle demande de raisonnement mathématique pur, la délégation a
pu reprendre : trois analyses indépendantes ont porté sur la condensation
massique, le certificat d'erreur et les singularités géométriques. La
limite dynamique par comptage des pôles a ensuite été soumise à une
contre-lecture. Les résultats figurent dans `QUOTIENT_LOCAL_PREUVES.md`.

Trois lectures primaires précisent l'antériorité et les limites. Le PDF de
Foster–Davis via CiteSeer ne s'ouvre pas ; le manuscrit disponible sur la
page de l'auteur a été consulté. Il possède des métadonnées provisoires ;
la référence publiée est confirmée par la liste de publications de l'auteur.
La version 1 de spaQR a été étudiée, sans extrapoler sa borne conditionnelle
au problème multicorps. Aucun code de solveur concurrent n'a été consulté.

Le prototype a été dérivé indépendamment comme élimination à droite par
Householder. La lecture de Foster–Davis confirme que le mécanisme de
perturbation des suppressions est antérieur : aucune nouveauté revendiquée.
Le programme expose le remplissage et un budget de travail, puis calcule
projection, injection/restriction implicites et réactions sélectionnées.
Ces réactions ne sont pas de norme minimale.

Onze contre-épreuves couvrent SVD de petits problèmes, permutations,
rotations locales, normalisation des équations, métrique de masse,
cinématique indépendante des cascades natives, grande chaîne sans
conversion dense, cas vides, seuil ambigu, budget et mobilité fausse
malgré un résidu petit. Les chronométrages exploratoires ont été écartés
de l'archive finale après ajout des contrôles de résultats non finis.

La campagne finale emploie la roue figée 0.9.0 pour les matrices de cascade.
Le prototype modal encore présent dans l'arbre de travail ne participe pas.
Chaque variante utilise un processus frais sur le CPU 8, un fil demandé,
avec un échauffement conservé et trois répétitions. L'ordre des variantes
est alterné. Les durées incluent normalisation, factorisation, projection
et réactions sur trois vecteurs ; elles excluent imports, assemblage du
modèle et oracle indépendant. Le pic RSS avant oracle inclut le modèle.
Les sources exécutées, empreintes et données brutes sont archivées.

Le sondage compare une brique au QR dense des normales, pas Vinkulum aux
meilleurs solveurs. Les cas défavorables sont conservés. Le constat de
remplissage du facteur R malgré des réflexions de support deux précise
la prochaine étape : ordonnancement par séparateurs, compression contrôlée
des facteurs et validation du rang. La borne physique complète et le
remplacement du calcul natif de réactions restent ouverts.

La contre-lecture croisée du carnet a corrigé les hypothèses de projecteur
orthogonal, de changement massique carré inversible et de compression
symétrique hors noyau. Elle a précisé le caractère ponctuel des injections
et l'absence de preuve globale de stabilité après report des pivots.
Elle a également distingué la statisation du résolvant de la congruence
conservant sa masse : cette dernière donne une queue d'ordre Ω⁴ au lieu
d'Ω² sous les hypothèses indiquées. L'équivalence avec Craig–Bampton
est explicitée, avec vérification de la documentation officielle Exudyn.
Ces dérivations ne sont pas présentées comme des méthodes nouvelles.
