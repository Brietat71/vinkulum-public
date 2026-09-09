# Éléments tiers

Apache-2.0 couvre les éléments originaux de Vinkulum. Elle ne remplace pas les
licences des fichiers tiers ci-dessous. Le dossier `python/vinkulum/licences/`
est distribué avec le paquet Python et contient les textes de licence.

## Modèles témoins

Le [manifeste des 36 fichiers XML/URDF](python/vinkulum/licences/modeles.json)
donne pour chaque fichier le dépôt public, le commit, le chemin d'origine,
l'empreinte SHA-256 et les notices associées. Les fichiers ont été comparés
octet par octet, ou par leur identifiant de blob Git, aux sources indiquées.
Ils sont redistribués sans modification ; les maillages externes ne sont pas
embarqués.

- 23 XML du [MuJoCo Menagerie](https://github.com/google-deepmind/mujoco_menagerie),
  sous les licences individuelles de chaque modèle, conservées dans
  `python/vinkulum/licences/modeles/`. Le projet Menagerie précise que les
  modèles ont des licences distinctes ; sa licence générale ne suffit pas.
- 11 URDF du dépôt [Bullet](https://github.com/bulletphysics/bullet3), avec
  sa notice générale et les notices spécifiques KUKA, Husky, humanoïde,
  quadrupède, racecar, Panda et Laikago dans `python/vinkulum/licences/bullet/`.
  Les notices déjà présentes dans les fichiers sont conservées.
- `acrobot.xml` provient de [dm_control](https://github.com/google-deepmind/dm_control)
  et `humanoid.xml` de [MuJoCo](https://github.com/google-deepmind/mujoco), sous
  Apache-2.0. Les deux textes officiels sont conservés séparément.

Attribution demandée par le témoin Laikago : Erwin Coumans et Yunfei Bai,
*PyBullet, a Python module for physics simulation for games, robotics and
machine learning*, 2016–2018, <http://pybullet.org>.

Les fichiers `fourbar.json` et `multibarmech.json` contiennent des résultats
numériques de référence issus des bancs MBDyn identifiés par leur champ
`source`. Les fichiers `.mbd` du dossier de confrontation modale sont générés
par les pilotes Vinkulum. Aucun code d'implémentation de MBDyn, Exudyn ou
MuJoCo n'est embarqué pour constituer le noyau.

## Dépendances de construction et d'exécution

Le [registre Rust](python/vinkulum/licences/rust.json) recense les 150 paquets
tiers du graphe Cargo verrouillé, y compris les dépendances de construction
et celles d'autres plateformes. Les notices sont conservées dans
`python/vinkulum/licences/rust/`. Les expressions SPDX à choix multiples
conservent les possibilités offertes par les auteurs ; le composant
`r-efi` offre notamment MIT ou Apache-2.0, sans imposer le choix LGPL.

NumPy, SciPy et les outils optionnels sont installés séparément par le
gestionnaire de paquets. Ils ne sont pas recopiés dans les sources Vinkulum
et conservent les licences fournies par leurs distributions. Les solveurs
externes utilisés dans les comparaisons s'installent séparément.

Les rapports citent des publications scientifiques et des documentations
publiques. Ces citations ne transfèrent aucun droit sur les œuvres citées.
