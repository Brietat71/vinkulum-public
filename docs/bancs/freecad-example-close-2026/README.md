# Fermer l’exemple pendant son cadrage dans FreeCAD

Le 10 septembre 2026, l’extension installée issue de `0493d34` reproduit un
**SIGSEGV natif** lorsqu’une fermeture de document est traitée pendant le cadrage
animé de l’exemple. La même recette passe après le correctif de `open_example()`.
Ce changement complète l’[extension FreeCAD existante](../../../apps/freecad/README.md).

Le cadrage animé de FreeCAD entre dans une boucle d’événements imbriquée. Une
fermeture et une nouvelle ouverture peuvent alors être traitées avant le retour
du premier cadrage. La trace conservée atteint `View3DInventorViewer::viewAll`
après destruction de la vue. L’extension effectue désormais l’orientation et
le cadrage initial immédiatement, puis restaure le réglage d’animation de la vue.

## Comparaison conservée

[Données et paquets testés](record.zip) · [Empreinte SHA-256](SHA256SUMS)

| Contrôle | Extension de référence | Extension corrigée |
|---|---|---|
| Cinq ouvertures avec fermeture programmée pendant le cadrage | SIGSEGV, code de sortie 1 ; la recette reste inachevée | Cinq ouvertures revenues, cinq documents fermés, code de sortie 0 |
| Douze contrôles existants de l’extension installée | Voir la [qualification initiale](../freecad-extension-010/README.md) | Réussis, code de sortie 0 |

La recette programme une fermeture 50 ms après chaque demande d’ouverture,
via les API réelles de FreeCAD et sa boucle d’événements Qt. Elle attend le
retour de toutes les ouvertures avant d’autoriser un verdict positif. Le lanceur
vérifie aussi le code de sortie du processus FreeCAD : un rapport incomplet ou
un crash ne peut pas être admis comme un succès.

Les douze contrôles conservés couvrent le menu dans un atelier natif, le calcul
réel, l’affichage, la réouverture d’un calcul et du `.FCStd`, le nettoyage de la
copie d’animation lors de la sauvegarde, la géométrie modifiée, un lancement
impossible, le refus d’un calcul concurrent, l’annulation et la fermeture du
document. Les deux calculs équivalents par changement de repère conservent des
écarts maximaux de 1,6461 × 10⁻¹³ m en position et 3,0182 × 10⁻¹³ sur les
composantes de rotation. Ce correctif ne modifie pas le calcul mécanique.

Les deux paquets sont installés dans des profils XDG séparés. L’archive de
référence provient du commit propre `0493d34e939d76e0e29ecb8a2e8e1ac535e04d26`.
Le paquet corrigé a été construit avant commit : son manifeste indique donc
`source_dirty: true` et identifie chaque fichier par SHA-256. Seul `host.py`
diffère parmi les sources distribuées ; son empreinte est
`1bddf14f1efddc5a9bec6c1aee96326d786c2989641cc2d5a501365bc9e017f2`.
Il s’agit d’un paquet de qualification, pas d’une nouvelle publication de version.

## Reproduire

Avec un FreeCAD Linux dédié et l’environnement moteur configuré :

```sh
python apps/freecad/package.py /tmp/freecad-candidate.zip
python ci/freecad_extension.py \
  --freecad /absolute/path/to/FreeCAD \
  --engine-python /absolute/path/to/engine/bin/python \
  --output /tmp/freecad-close-fresh \
  --archive /tmp/freecad-candidate.zip --recipe example-close
```

Utiliser un dossier de sortie nouveau. Le lanceur crée le profil isolé, installe
le paquet et utilise Xvfb si aucun écran n’est disponible. Retirer
`--recipe example-close` pour rejouer les douze contrôles de l’extension dans
un autre dossier neuf. Pour le contre-exemple, extraire `record.zip` et fournir
`before/Vinkulum-FreeCAD-0.1.0a1.zip` : le lanceur doit rendre un code non nul.

La qualification utilise l’AppImage officiel FreeCAD 1.1.3 Linux x86-64,
Python 3.11.14 et OCCT 7.8.1, sans modification de FreeCAD. Le moteur séparé est
Vinkulum 0.20.0 / services Studio 0.6.0a2.dev6 / OCCT 8.0.1.0. Les empreintes,
rapports, journaux de processus et recettes figurent dans `record.zip`.
Les autres versions et plateformes ne sont pas qualifiées par cette comparaison.

Recette et correctif sont des contributions originales sous
[Apache-2.0](../../../LICENSE). Les paquets conservent la licence et les notices
de l’extension ; le runtime FreeCAD n’est pas inclus.
