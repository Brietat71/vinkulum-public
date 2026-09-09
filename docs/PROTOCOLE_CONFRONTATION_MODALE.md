# Protocole de confrontation modale à précision commune

Ce protocole compare les interfaces publiques de Vinkulum 0.17.0, Exudyn
1.11.0 et une compilation identifiée de MBDyn develop. Il porte sur la
recherche des premières fréquences d'une chaîne articulée à l'équilibre.
Il ne permet pas de classer globalement des solveurs multicorps généralistes.

## Modèle et critère fixés avant les mesures finales

La chaîne possède `n` barres, longueur totale 1 m et masse totale 1 kg.
Chaque barre mesure `1/n`, pèse `1/n` et a une inertie autour du pivot Y
`m*l²/12`. Les pivots sont sans amortissement ni gravité ; leur ressort de
rotation vaut `10*n` N·m/rad. La base est fixée. Les énergies linéarisées sont
identiques dans les trois descriptions publiques ; aucune équivalence des
lois non linéaires à grandes rotations n'est revendiquée.

Les tailles sont **4, 16, 64, 256 et 1 024**. La demande commune porte sur
`min(6,n-1)` fréquences. La référence indépendante est l'assemblage des
énergies en coordonnées angulaires de `ci/modeles_modes.py`, déjà employé
pour la qualification modale. Sa construction ne fait pas partie du temps
d'un solveur. Cette référence est flottante, sans borne spectrale formelle.

Le seuil d'erreur relative de fréquence est fixé à **1e-8**. Pour le spectre
complexe MBDyn, la partie réelle divisée par la pulsation doit aussi rester
sous **1e-8** : le modèle analytique est conservatif et strictement stable.
Les fréquences sont ordonnées, sans appariement au plus proche qui pourrait
masquer une fréquence manquante ou répétée.

## Configurations conservées

- Vinkulum : `Noyau.modes_creux`, réglages publics par défaut.
- Exudyn : coordonnées cartésiennes avec contraintes et solveur dense ;
  arbre cinématique avec solveur dense ; arbre cinématique avec solveur creux.
- MBDyn : UMFPACK pour le système linéaire ; LAPACK avec permutation ;
  ARPACK avec recherche courte ; ARPACK avec recherche presque complète.

L'interface cartésienne creuse d'Exudyn refuse les contraintes algébriques
dans les sondes préparatoires ; la variante **arbre creux est bien incluse**.
Les options viennent de l'[interface publique Exudyn 1.11.0](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/cInterface/MainSystem.html).

MBDyn fournit un spectre complexe transformé. Les résultats `alpha` et
`dCoef` sont reconvertis suivant le [manuel public, analyse propre](https://www.mbdyn.org/userfiles/documents/mbdyn-input-1.7.3.pdf).
Le paramètre de transformation est fixé à **0.02** sur toute la grille.
La fenêtre de sortie est définie par la référence : moitié de la première
fréquence, jusqu'au milieu entre la dernière demandée et la suivante.
Cette information a priori réduit les sorties de vecteurs imposées par
l'interface pour obtenir des valeurs avec 17 chiffres. Elle est déclarée,
et non déduite après lecture d'un résultat concurrent.

Avec `k=min(6,n-1)`, la recherche ARPACK courte utilise `nev=2*k+4`,
`ncv=max(32,4*k+8)`, `tol=1e-12`. La recherche élargie utilise `nev=N-4`,
`ncv=N`, où `N=17*n+18` est la dimension du système MBDyn de ce modèle.
Ces deux réglages sont conservés même lorsqu'ils échouent. Les sorties hors
fenêtre sont comptées. Les triplets nuls et valeurs infinies des contraintes
ne sont pas assimilés à des fréquences physiques.

## Mesures et refus

Chaque configuration et taille utilise **un échauffement et trois répétitions**,
chacun dans un processus neuf. L'ordre des moteurs tourne entre répétitions.
Tous les processus utilisent le même CPU autorisé et un seul fil demandé.
Le temps couvre le processus complet : imports ou lecture du modèle,
construction, assemblage, analyse et sortie. Les différences d'interface font
partie de cette mesure ; elle n'isole pas le coût de la factorisation.
Les temps internes Python sont conservés séparément, sans leur attribuer
une équivalence avec un temps interne MBDyn.

Le plafond commun est **60 s** et **4 Gio d'espace virtuel** par exécution.
Un dépassement est conservé comme refus sous ce budget, sans préjuger de la
capacité du solveur avec davantage de ressources. Le pic RSS est mesuré par
GNU time ; un processus tué peut ne pas le fournir. La charge du système est
conservée. Les journaux natifs MBDyn, susceptibles de contenir l'environnement,
ne sont pas publiés ; les résultats numériques et les modèles le sont.

Une configuration n'entre dans le classement d'une taille que si
**l'échauffement et les trois répétitions sont qualifiés**. La médiane et
l'étendue des trois temps mesurés sont alors rapportées avec le pic RSS.
Aucun temps d'une sortie imprécise, incomplète, refusée ou interrompue ne
sert à annoncer un gain. Trois répétitions décrivent cette campagne ; elles
ne prouvent pas une garantie de répétabilité universelle.

## Qualification de la collecte

Une campagne préparatoire a été invalidée : sans extension explicite, le
binaire MBDyn tronque le chemin spectral au dernier point rencontré, même
lorsque ce point appartient à un répertoire. Les essais écrivaient alors
sur le même fichier hors dossier. Le producteur emploie maintenant un nom
avec extension explicite, vérifié dans un répertoire contenant des points.
Ce défaut de collecte n'est pas compté comme un défaut de précision de MBDyn.
La campagne finale repart intégralement ; les seuils numériques sont conservés.

Le vérificateur d'archive recalcule les références, les fréquences et les
décisions, puis le classement. Il vérifie l'exhaustivité de la grille et les
empreintes des fichiers. Il ne peut pas prouver a posteriori un temps ou
une consommation mémoire : ce sont des observations de la machine identifiée.
