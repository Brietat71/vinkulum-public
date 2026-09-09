# Vinkulum Studio 0.1.0 — prototype G0

Interface de bureau expérimentale PySide6 : paramètres d'un pendule, calcul réel
Vinkulum dans un processus séparé, animation 2D, angle en fonction du temps,
arrêt, paramètres JSON versionnés et export CSV avec provenance.

Le noyau reste indépendant. Ce paquet ne modifie ni sa physique ni sa version.
La vocation généraliste est celle de la suite ; le pendule limite seulement G0.

## Installation et lancement

Python **3.14** et le noyau **0.18.1** sont nécessaires. Dans un environnement
virtuel, depuis la racine du dépôt :

```sh
python3.14 -m venv .venv-studio
. .venv-studio/bin/activate
python -m pip install .
python -m pip install ./apps/studio
vinkulum-studio
```

La première installation compile le noyau et nécessite Rust/Cargo. Une roue
Vinkulum 0.18.1 **compatible avec la plateforme et Python 3.14** peut remplacer
`pip install .`. La roue Linux publiée ne fonctionne pas sur macOS. PySide6
est installé séparément, sous les licences de son éditeur ; aucun installateur
embarquant Qt n'est fourni ici. Le code original de Studio est sous Apache-2.0,
comme le dépôt. L'installation nécessite un accès réseau aux dépendances ;
l'application installée fonctionne localement sans réseau.

## Premier parcours

1. Ouvrir Studio puis cliquer **Lancer le calcul**.
2. Observer la courbe, utiliser **Lecture / Pause** et le curseur temporel.
3. Modifier la longueur et relancer. Pendant le calcul les champs restent éditables.
4. Cliquer **Arrêter** : le dernier résultat réussi reste affiché et identifié.
5. Enregistrer les paramètres, les modifier, puis recharger le JSON.
6. Exporter le CSV : sa première ligne commentée contient les entrées du résultat
   affiché, les unités et la provenance. La sauvegarde JSON concerne le brouillon.

Un seul calcul est admis à la fois. `Ctrl+Entrée` lance le calcul. Les champs,
boutons et le curseur sont accessibles au clavier. Aucun calcul n'est lancé
à l'ouverture ni au chargement d'un fichier.

## Domaine et limites

Le montage reprend le pendule de `python/vinkulum/verification.py` : masse avec
inertie isotrope de 10⁻⁸ kg·m², rotule au point fixe, gravité 9,80665 m/s²,
vitesse initiale nulle, intégrateur α-généralisé avec ρ∞ = 0,9.

Bornes de saisie propres à G0 : longueur 0,05–20 m, masse 0,01–100 kg,
angle initial −170° à 170°, durée 0,001–600 s, pas 10⁻⁶–1 s, pas ≤ durée,
au plus 20 000 pas demandés. **Ces bornes ne garantissent pas la convergence
ni la précision** ; le noyau peut refuser un calcul. L'interface vérifie la
finitude, la chronologie, la géométrie et l'association aux entrées, sans
en déduire une borne d'erreur physique. Le statut scientifique est `NotAssessed`.

Les samples exportés sont ceux du solveur avec l'état initial ajouté explicitement.
La lecture sélectionne les échantillons à l'horloge réelle sans interpolation ;
la courbe relie visuellement les points. Le dessin est remis à l'échelle selon
la longueur capturée dans le résultat. Il ne représente pas une géométrie CAO.

Arrêter demande la terminaison du processus, puis force son arrêt après une
seconde si nécessaire. Un arrêt ne crée ni résultat réussi ni checkpoint.
La fermeture termine le worker. Le processus n'est pas un bac à sable de sécurité.
Deux threads Rayon et un thread BLAS/OpenMP sont demandés. La taille des fichiers
d'entrée et de résultat est bornée ; ce mécanisme ne constitue pas une limite
stricte de mémoire CPU imposée par le système d'exploitation.

Le JSON de paramètres utilise `format: vinkulum-studio-parameters`,
`schema_version: 1` et les cinq champs de `Parameters`. Champs inconnus, clés
dupliquées, valeurs non finies et versions inconnues sont refusés. Les sauvegardes
et exports remplacent le fichier seulement après écriture complète et `fsync`
du fichier temporaire ; aucune garantie universelle après perte d'alimentation
du disque n'est revendiquée. Le dernier résultat reste en mémoire pendant la
session ; seuls les paramètres et l'export sont persistés.

## Vérification

```sh
python -m pip install './apps/studio[test]'
QT_QPA_PLATFORM=offscreen python -m unittest discover -s apps/studio/tests -v
```

Voir `docs/STUDIO_G0.md` à la racine pour la recette, les preuves produites et
les plateformes effectivement testées. **La recette macOS ARM64 reste à faire
sur une machine réelle.** G0 n'apporte aucune certification générale du noyau.
