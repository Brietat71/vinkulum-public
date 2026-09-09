# Studio 3D — audit des critères de livraison

Audit du 9 septembre 2026, Studio 0.2.0 / noyau 0.19.0. La qualification
complète reste ouverte : les essais sur bureaux physiques ne sont pas disponibles
dans cet environnement Linux avec Xvfb. Le dernier tour de développement a
produit du code, des corrections de régressions et des preuves ; ce n'était
pas une simple attente.

| Critère du plan | Preuve disponible | Conclusion |
|---|---|---|
| Document rigide immuable, UUID, historique, JSON, import G0 | `test_document.py`, installation indépendante | Vérifié dans le domaine borné du document |
| Boîte, cylindre, sphère, sélection et pose | `test_viewport3d.py`, `test_editor3d.py`, captures des deux mécanismes | Fonctionnement automatisé vérifié ; confort et gestes physiques à qualifier |
| Repères indépendants, quatre types de liaison | Adaptateur `mechanism.py`, validations document, contre-épreuves natives | Implémenté ; conservation des repères et refus des ancrages incohérents testés |
| Lois constante, linéaire, table ; commandes et charges | `test_document.py`, `test_mechanism.py`, `test_charges_temporelles.py` | Loi et physique testées ; parcours humains de chaque dialogue restant ouverts |
| Tangente de force déportée et domaines d'analyse | Contre-épreuve Rust, audit jacobien, référence DOP853, refus statique/conservation | Vérifié sur les cas définis ; pas de théorème général supplémentaire |
| Calcul dans un worker, annulation et panne | Tests G0 du contrôleur et recettes 3D de calcul/annulation | Vérifié automatiquement |
| Résultat associé à son projet, scène/arbre/propriétés cohérents | Régression de remplacement de résultat ; saisie et déplacement pendant calcul | Vérifié automatiquement, y compris le défaut détecté dans les captures |
| Positions/rotations/vitesses, coordonnées de liaison, animation et export | Tests adaptateur, tableaux immuables, recettes graphiques, CSV | Vérifié ; angle principal, sans comptage de tours |
| Échanges bornés, sans pickle | Validation NPZ et test d'en-tête déclarant 10¹² éléments avec petite charge utile | Rejet avant allocation vérifié ; pas de limite mémoire système revendiquée |
| Références physiques indépendantes | Lagrange double pendule, fermeture bielle-manivelle, translation et rotation analytiques | Tolérances et convergence satisfaites pour ces références |
| CI noyau complète et bancs | `bancs/studio-3d/ci-local.log` | Réussie localement |
| Roues installées hors checkout | Manifestes, concordance des sources, 29 tests Studio, 41 vérifications noyau et 8 nouveaux tests | Réussie sous Linux x86_64 / CPython 3.14 |
| Versions, API et licences | Manifestes 0.19.0/0.2.0, génération API, vérificateur des notices | Vérifié |
| Rendu 1080p et objectif 60 images/s | `desktop.json`, 120 mesures après chauffe, Mesa llvmpipe | Objectif observé sur scènes de 3/32 corps sans liaison ni charge ; aucune garantie générale |
| Linux interactif sur bureau physique | Aucune recette humaine disponible | **Non vérifié** |
| macOS ARM64, rendu et interaction réels | Archive d'essai et commandes fournies ; aucune exécution Mac reçue | **Non vérifié** |
| CI GitHub de cette évolution | Workflow préparé et actionlint valide ; pas de push de cette évolution | Exécution distante non effectuée |

La recette manquante se trouve dans [STUDIO_3D.md](STUDIO_3D.md). Le résultat
automatisé sous Xvfb ne constitue pas une preuve d'exécution sur macOS ni une
preuve de confort d'utilisation. Une réponse concernant la disponibilité d'un
Mac a été demandée ; l'absence de réponse n'est pas une validation.

L'implémentation n'est pas déclarée achevée au sens du plan multiplateforme.
Les artefacts sont prêts à être essayés ; les limites de qualification restent
explicites. Aucun ancien résultat scientifique ni aucune donnée de recette
humaine n'est remplacé par une hypothèse de succès.

## Vérification de l'archive d'essai

L'archive `vinkulum-studio-0.2.0-essai.tar.gz` (4 294 718 octets) a été extraite
dans un dossier distinct. Ses 50 fichiers sources recensés par le manifeste
correspondent aux empreintes testées. La construction hors ligne des deux roues
a réussi depuis cette extraction, avec les dépendances de compilation déjà
en cache : [noyau](bancs/studio-3d/archive-build.log),
[Studio](bancs/studio-3d/archive-studio-build.log).

SHA-256 de l'archive :
`684a3da66bf2471848537abe0bee81723f106c4dd7ddd6ef4775ca02dd639aed`.
Les anciennes archives scientifiques volumineuses sont omises, comme annoncé
par `LIRE_AVANT_ESSAI.txt`. Cette archive permet la compilation et la recette
Studio ; elle ne permet pas la CI historique complète sans les archives du dépôt.
Ce contrôle de construction Linux ne démontre pas la compilation macOS.
