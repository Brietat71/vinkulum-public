# Studio 3D — suivi d'implémentation

Objectif accepté : Studio 0.2.0, éditeur de mécanismes rigides 3D, avec noyau
0.19.0. Primitives boîte/cylindre/sphère, sélection et manipulateurs, liaisons
pivot/rotule/glissière/encastrement, lois de mouvement et charges temporelles,
document persistant, annulation/rétablissement, simulation isolée et provenance.

Les déplacements de conception ne résolvent pas les contraintes : les
incohérences sont diagnostiquées. CAO, URDF, collaboration, contact et flexibles
ne font pas partie de ce lot. Les documents historiques G0 restent lisibles.

## Contrôles de livraison à compléter

- [x] Essai initial Qt 6.11.2 / VTK 9.7.0 sous Xvfb/OpenGL, Python 3.14.7 Linux.
- [x] Document mécanique, lois, historique, sauvegarde et migration G0.
- [x] Vue 3D, sélection, navigation et manipulateurs : contrôles automatisés réussis.
- [x] Liaisons natives à deux repères explicites et contre-épreuves.
- [x] Efforts temporels natifs, tangentes, diagnostics et références physiques.
- [x] Éditeur complet : arbre, propriétés, création et diagnostics.
- [x] Worker multicorps, résultats, animation, courbes et exports.
- [x] Recettes double pendule, mécanisme motorisé/glissière et charge temporelle.
- [x] Roues, versions, API, notices, CI complète et documentation.
- [ ] Recette bureau Linux interactif et macOS ARM64.

Ce suivi n'est pas une déclaration de livraison. La qualification sans écran
ne remplace pas celle des bureaux réels. Le checkout privé reste préservé.

## Environnement

Checkout : `/tmp/vinkulum-publication`, base `d2191b4`.
Environnement de développement Studio : `/tmp/vinkulum-studio-env`.
Qt/X11 nécessite libxcb-cursor0, libxcb-icccm4, libxcb-keysyms1,
libxcb-image0, libxcb-render-util0 et libxcb-util1. Dans le sandbox courant,
l'accès au socket X11 de Xvfb requiert une exécution escaladée.
Premier essai : rendu `vtkXOpenGLRenderWindow`, image 640 × 480 vérifiée.

## Qualification locale terminée

Voir [la recette détaillée](STUDIO_3D.md) et ses journaux : CI noyau complète,
29 tests Studio sur roues installées hors checkout, références physiques,
captures et mesure de rendu à 1080p. Les bureaux humains restent à qualifier ;
le workflow modifié est préparé mais pas encore exécuté sur GitHub.

L’[audit des critères](STUDIO_3D_AUDIT.md) distingue les fonctionnalités testées
des validations de bureau encore absentes.
