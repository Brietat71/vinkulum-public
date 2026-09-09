# Numérotation des versions

Vinkulum utilise trois nombres : **MAJEURE.MINEURE.CORRECTIVE**, par exemple
`0.5.1`. Le paquet Rust et le paquet Python portent toujours le même numéro.
Le tag Git correspondant ajoute le préfixe `v` : `v0.5.1`.

## Choisir le numéro

La règle dépend de l'effet pour les utilisateurs, pas du nombre de lignes
modifiées ni du temps passé. Si un lot contient plusieurs types de changement,
le changement qui exige la plus forte augmentation détermine le numéro.

| Changement | Pendant la série `0.x` | À partir de `1.0.0` |
|---|---|---|
| Correction d'un défaut, optimisation ou réorganisation interne conservant le contrat public | Corrective : `0.5.1` → `0.5.2` | Corrective : `1.2.3` → `1.2.4` |
| Nouvelle fonctionnalité ou extension compatible de l'API | Mineure : `0.5.1` → `0.6.0` | Mineure : `1.2.3` → `1.3.0` |
| Rupture de compatibilité d'un usage documenté et valide | Mineure : `0.5.1` → `0.6.0`, avec indications de migration | Majeure : `1.2.3` → `2.0.0`, avec indications de migration |
| Documentation, commentaires ou tests seuls, sans modification du paquet distribué | Commit ; pas de nouvelle version nécessaire | Commit ; pas de nouvelle version nécessaire |

Une augmentation de la mineure remet la corrective à zéro. Une augmentation
de la majeure remet les deux autres nombres à zéro.

La série `0.x` correspond au développement initial : les ruptures restent
possibles, mais elles sont annoncées et regroupées dans une nouvelle mineure.
Le passage à **1.0.0** est une décision de stabilisation du contrat public,
appuyée par sa documentation et ses validations. Il ne découle ni du nombre
de fonctionnalités ni d'un nombre donné de versions `0.x`.

## Compatibilité et résultats numériques

Le contrat public couvre les API Rust et Python destinées aux utilisateurs,
les formats d'entrée et de sortie, les unités, les conventions de repère,
les comportements documentés et les environnements pris en charge.
Supprimer une méthode, modifier une unité ou relever la version minimale de
Python constitue une rupture. Une modification d'une valeur par défaut qui
change un usage valide constitue également une rupture.

Une correction peut modifier des résultats numériques auparavant erronés :
elle reste corrective si elle rétablit le comportement attendu et documenté.
Elle doit être accompagnée d'un cas de reproduction et d'une validation
indépendante lorsque la physique est concernée. Une optimisation ne promet
pas des résultats identiques bit à bit ; elle doit conserver la précision et
les propriétés physiques vérifiées. Un changement volontaire de modèle
physique ou de convention pour des entrées valides relève d'une rupture.

Refuser explicitement une configuration invalide ou non prise en charge,
auparavant acceptée en silence, peut faire partie d'une corrective. Ce refus
doit être signalé dans les notes de version, avec son déclencheur et sa raison.
Cette exception ne permet pas de retirer une fonctionnalité prise en charge.

**Exemple : `0.5.0` → `0.5.1`.** L'[audit du noyau](AUDIT_NOYAU_2026-09-06.md)
corrige la dynamique des superéléments et la gestion des échecs, en conservant
les signatures et les formats de sortie. Il explicite aussi les refus ajoutés,
dont `alpha != 0` pour les superéléments, auparavant sans effet. Ce lot justifie
une corrective. Une nouvelle capacité du solveur justifierait une mineure.

## Préparer une version

1. Choisir le numéro selon les règles ci-dessus. Consigner les changements,
   leurs effets sur les utilisateurs, les éventuelles migrations et les
   validations dans une note de version ou une note d'audit liée au README.
2. Modifier ensemble `Cargo.toml` (`package.version`) et `pyproject.toml`
   (`project.version`). Reconstruire et réinstaller le paquet, par exemple
   avec `.venv/bin/maturin develop --uv --release` ; Cargo actualise l'entrée
   `vinkulum` de `Cargo.lock`, à inclure dans le commit.
3. Vérifier que `vinkulum.__version__` et les métadonnées du paquet installé
   affichent le nouveau numéro. `__version__` vient de ces métadonnées ;
   `bancs.VERSION` le reprend automatiquement. Aucun de ces deux champs ne
   doit recevoir une copie manuelle du numéro.
4. Actualiser les mentions de version courante du README et régénérer
   `docs/API.md` avec `.venv/bin/python -m vinkulum.doc --ecrire`, après
   réinstallation. Conserver les numéros des références historiques et
   l'identité des versions réellement mesurées dans les rapports de bancs.
5. Exécuter `ci/local.sh` pour la livraison. Pour un changement du solveur,
   des contacts ou de la physique, ajouter les régressions pertinentes et
   les bancs étendus (`ci/local.sh --bancs`). Consigner les résultats et les
   limites observées, y compris les régressions de performance.
6. Créer le commit de version après validation, puis un tag annoté sur ce
   commit : `git tag -a vX.Y.Z -m "Vinkulum X.Y.Z"`, en remplaçant `X.Y.Z`
   par le numéro choisi. Le tag doit désigner le code et les métadonnées
   effectivement validés.

## Commits, tags et diffusion

Un **commit** enregistre une modification ; une **version** identifie un état
livrable du paquet ; un **tag** désigne son commit. Plusieurs commits peuvent
préparer une même version. Une modification documentaire après `v0.5.1` peut
donc être commitée sans créer `0.5.2` ni déplacer `v0.5.1`.

Un tag de version est immuable : ne pas le déplacer pour y ajouter une
correction. Toute livraison modifiée du paquet reçoit un nouveau numéro ;
ne pas remplacer une roue ou une archive déjà diffusée sous le même numéro
par une révision différente. Des constructions du même état pour plusieurs
plateformes ou interpréteurs partagent naturellement le numéro de version.

Créer un commit et un tag locaux ne pousse rien vers un dépôt distant et ne
publie aucun paquet. La diffusion est une opération distincte. Le paquet est
public à partir de la 0.18.1, sous Apache-2.0 pour ses éléments originaux,
et n'est publié sur aucun index de paquets.
