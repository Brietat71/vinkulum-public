# Vinkulum 0.8.1 — globalisation de Newton statique

7 septembre 2026 · Version corrective, API conservée.

La statique peut accepter pendant sa première tentative un pas qui réduit
la correction de pose prédite, même si les forces augmentent, à condition
de contrôler séparément la fermeture des liaisons. Cette règle évite les
reprises systématiques sur les paliers Princeton. Le filtre géométrique
préserve les chaînettes que la seule contraction des corrections faisait
diverger dans une variante exploratoire.

Sur les six rampes Princeton de la campagne interne, le temps médian est
divisé par 3.17 à 3.29 par rapport à la roue 0.8.0. Les cascades de
32 boucles prennent environ 13 à 14 % de temps en moins. Les chaînettes
restent proches de leur coût précédent. Ces rapports proviennent de
306 résolutions chronométrées sur 51 configurations, avec trois
répétitions et un échauffement par processus.

Les signatures et les critères physiques de convergence sont conservés.
`strict=True` exige la tolérance ; les refus restaurent l'état. La
recherche de pas peut suivre un autre chemin et le nombre d'itérations
change. Le mode sans amortissement et l'option `sigma_lie` restent dans
leurs domaines documentés en 0.8.0.

Les [équations, sources et limites](GLOBALISATION_STATIQUE.md) distinguent
les expériences internes des confrontations externes. Deux configurations
de petite échelle restent non résolues au réglage resserré du diagnostic,
comme en 0.8.0, et coûtent environ 61 % de temps supplémentaire avant le refus.
Aucun gain général face à MBDyn ou à Simpack n'est déduit de ce seul lot.

Deux régressions Python ajoutent des contrôles indépendants de résultante,
de moment et de polygone funiculaire discret. La CI locale et son modèle
GitHub exécutent désormais les quatre modules de tests Python, afin que
ces contrôles accompagnent les futures modifications du noyau.

Validation : 70 tests Rust, 8 du prototype, 78 tests Python, 41 groupes
de vérification, 46 bancs mécaniques et 9 bancs de contact. La roue finale
est également vérifiée dans un environnement neuf, hors du dépôt.
[Journaux et empreintes](bancs/version-0.8.1.json).
