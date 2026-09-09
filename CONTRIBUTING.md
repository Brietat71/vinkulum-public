# Contribuer à Vinkulum

**New here?** Start with the [contributor projects](docs/CONTRIBUTOR_PROJECTS.md):
concrete first deliverables for CAD, dynamics, Pinocchio, numerical verification,
Rust performance, Qt and documentation. Contributions can be discussed in English
or French. Describe the problem, your reference and how the result can be checked.

Vinkulum est un noyau multicorps généraliste en phase alpha. Les contributions
peuvent porter sur un défaut reproductible, une référence physique indépendante,
une preuve, une documentation ou une amélioration mesurée des performances.

Pour Studio, suivre [l'installation locale](apps/studio/README.md) et exécuter
`PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh`. Les modifications de la CAD
doivent utiliser [l'environnement OCCT 8 adapté](docs/STUDIO_CAD.md). Les retouches
Python/Qt se vérifient depuis les sources ; réserver le binaire autonome aux
livraisons. Pour un nouveau solveur, commencer par un contrat de conversion et
un cas de référence, comme dans le [plan Pinocchio](docs/PINOCCHIO_INTEGRATION.md).

Ouvrir une issue avec la version, la plateforme, un modèle minimal, le résultat
observé et le résultat attendu avec sa source. Pour un problème numérique,
préciser les unités, les tolérances et les pas de temps. Ne pas joindre de
secrets ni de données confidentielles.

Installer le projet selon le README et préparer Lean 4.19.0/Mathlib selon
[`preuves/README.md`](preuves/README.md), puis exécuter `ci/local.sh` avant une
pull request. Pour une modification mécanique, ajouter `ci/local.sh --bancs`
et les contre-épreuves pertinentes. Le workflow GitHub Actions appelle la
même CI avec `--bancs` et teste une roue installée hors du dépôt. Les campagnes externes
peuvent exiger des solveurs installés séparément. Décrire les contrôles exécutés
et ceux qui n'ont pas pu l'être.

Les implémentations originales s'appuient sur la littérature, les interfaces
publiques, les modèles d'entrée et les observations. Ne pas copier de code
d'un solveur tiers pour implémenter Vinkulum. Identifier la provenance et les
licences de tout fichier tiers proposé.

Les contributions intentionnellement soumises sont proposées sous Apache-2.0,
conformément à sa section 5, sauf indication explicite contraire. Ne soumettre
que du contenu que vous êtes autorisé à distribuer.
