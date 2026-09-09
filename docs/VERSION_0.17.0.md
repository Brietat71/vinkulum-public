# Vinkulum 0.17.0 — analyse modale creuse qualifiée

`Noyau.modes_creux` rend accessible le chantier modal préservé dans le
checkout de travail, avec une qualification et des contrôles supplémentaires.
Il calcule les valeurs propres signées proches d'un décalage, les formes
normalisées en masse et les réactions de contraintes.

Le raffinement de Rayleigh–Ritz accumule le petit problème dans les matrices
physiques originales avec précision étendue lorsqu'elle est disponible.
Sur la console à 1 024 éléments, l'erreur relative de fréquence face à la
référence indépendante passe d'environ `2,46e-8` dans le brouillon à `4,9e-11`.
Les matrices et les résultats avant/après sont archivés. Le brouillon n'était
pas une API publiée dans une version antérieure.

La fonction vérifie aussi l'orthogonalité dans la masse originale et refuse
une masse couplée entre corps, incompatible avec sa normalisation. Le rang
et sa réduction sont déclarés numériques ; aucun certificat spectral n'est
annoncé. Le [contrat modal](MODES_CREUX.md) décrit les replis denses, la
précision dépendant de la plateforme et les limites de l'interprétation physique.

## Qualification isolée

Roue : `vinkulum-0.17.0-cp314-cp314-manylinux_2_39_x86_64.whl`.
SHA-256 : `c2835cd2b900f2ef2e1610e28f3960d4b4ed51f5246d3770d4ee6cce532a58a9`.
Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1.

- 209 tests Python réussis depuis le paquet installé, dont 14 tests modaux.
- 46/46 bancs rapides et 9/9 cas de contact réussis sur cette roue isolée.
- 18 cas à référence indépendante, trois familles et six tailles jusqu'à
  1 024 éléments, passent le seuil relatif commun de fréquence `1e-8`.
- Les résidus physiques et l'orthogonalité sont recalculés depuis les archives,
  sans appeler le solveur modal ; trois tests de relecture et d'altération.

Le [manifeste](bancs/modes-creux-0.17.0/qualification.json) identifie la roue,
son extension native, les sources de qualification et les 18 archives.
Il s'agit de qualification numérique, pas d'une preuve d'inclusion spectrale.
Les contrôles historiques et les preuves restent obligatoires avant push.

## Travaux ouverts

Un complément documentaire compare désormais cette roue à MBDyn et Exudyn
sur [140 essais modaux](CONFRONTATION_MODALE.md), à précision commune et
configurations identifiées. Le diagnostic des repères à temps long reste ouvert. Cette livraison ne clôt pas le plan
de fiabilité et ne certifie pas le noyau entier.
