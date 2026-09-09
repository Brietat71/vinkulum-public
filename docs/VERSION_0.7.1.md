# Vinkulum 0.7.1 — raffinement creux des équilibres redondants

Cette corrective accélère le Newton statique des mécanismes dont les
contraintes présentent des dépendances générales. Elle corrige aussi
l'amplification numérique causée par une diagonale presque nulle dans
l'ancien équilibrage. Les signatures Rust et Python, modèles physiques,
unités et critères publics de convergence sont conservés.

Le numéro **0.7.1** suit les [règles du projet](VERSIONNEMENT.md) : il s'agit
d'une correction et d'une optimisation internes, sans extension de l'API.

## Changements

L'équilibrage utilise les normes des blocs de trois coordonnées spatiales,
avec des facteurs en puissances de deux. Il réutilise les échelles entre
itérations quand elles restent représentables, et évite ainsi le surcoût
des premières tentatives d'implémentation.

Après un échec de résolution directe sur un grand système contraint,
une factorisation auxiliaire creuse permet de corriger la solution par
raffinement. Chaque résidu est recalculé avec produits et sommes compensés
sur la matrice originale. La tentative est rejetée si elle stagne ou
diverge ; le repli orthogonal reste disponible. Le LU dense intermédiaire
devenu inutile sur ce chemin est supprimé.

Le [rapport technique](RAFFINEMENT_STATIQUE.md) présente les équations,
références scientifiques, matrices capturées, essais négatifs et contrôles
physiques. Le [bilan des mesures](bancs/raffinement-bilan.json) couvre
26 configurations, dont les cascades tournées, avec trois répétitions
alternées entre la roue 0.7.0 et la nouvelle extension.

La cascade connexe de **96 corps passe de 2,038 s à 37,65 ms**, soit
**54,1×**, avec **86,0 → 40,9 Mio** de pic mémoire. Le même modèle tourné
gagne **3,34×** et les 32 parallélogrammes indépendants **79,6×**. Les
156 équilibres atteignent la tolérance stricte et passent les contrôles
physiques indépendants. Des surcoûts restent mesurés, jusqu'à **6,8 %**
sur les autres configurations ; aucun gain universel n'est annoncé.

## Compatibilité

Aucune adaptation des appels Python ou Rust n'est nécessaire depuis
0.7.0. Python **3.14 ou plus** reste requis et l'API compte **185 entrées**.
Les dépendances ne changent pas. Le paquet doit être reconstruit ou la
nouvelle roue installée pour employer le nouveau noyau.

Une transformation d'échelle qui perdrait un coefficient, une composante
du second membre ou de la correction est désormais explicitement refusée.
Ce refus suit la restauration statique existante. Les résultats numériques
ne sont pas promis identiques bit à bit ; le verdict final reste fondé
sur les résidus mécaniques complets.

## Validation et livraison

La CI étendue passe sur le code 0.7.1 : formatage, Clippy sans avertissement,
**57 tests Rust**, **67 tests Python**, **41/41 vérifications**,
**46/46 bancs mécaniques**, **9/9 contacts** et référence d'API à jour.
Les 400 paliers Princeton restent stricts et les trois essais de translation
de `2⁻⁶⁰ m` conservent exactement leur petite composante.

Les [preuves du noyau](bancs/raffinement-validation.json) identifient
les sources, extensions et artefacts contrôlés. Le
[relevé de livraison](bancs/version-0.7.1.json) identifie la roue installée
dans un environnement neuf hors du dépôt, ses métadonnées, les empreintes
`RECORD` et les contrôles d'usage. Le tag annoté `v0.7.1` désigne le commit
de cette livraison. Le paquet reste privé, sans publication sur un index.

## Limites

Les projecteurs des composantes réellement déficientes restent denses,
ce qui limite notamment le gain du cas tourné. Le raffinement n'est pas
une garantie de norme minimale pour toute correction intermédiaire ;
les réactions finales conservent leur projecteur QR/COD. Les modes libres
et les contraintes incompatibles peuvent imposer le repli.

Les matrices d'analyse publiques gardent leurs limites, dont l'échec de
`k_c_m_z()` sur la cascade initiale de 96 corps. La poutre intégrée reste
expérimentale ; le câble condensé extrême garde son échec connu. Les
surcoûts mesurés sont publiés avec les gains.

La confrontation MBDyn exécutée reste celle de 0.6.2 ; aucune exécution
Simpack n'est archivée. L'[objectif généraliste](OBJECTIF_MBDYN.md) reste
ouvert et les gains internes ne prouvent pas une supériorité globale.
