# FreeCAD comme interface de Vinkulum

La priorité de développement est de s’appuyer sur FreeCAD. Les nouveaux parcours
de modélisation et d’analyse sont construits dans un atelier FreeCAD ; le
développement de l’interface maison de Studio est mis de côté.

FreeCAD porte les esquisses, la CAO paramétrique, le document, la sélection et la
navigation. Vinkulum apporte son moteur mécanique, les adaptateurs de calcul,
les références physiques et la traçabilité des résultats. Les fonctions utiles
des paquets existants sont réutilisées sans ouvrir l’interface Studio. Le
processus de calcul garde son propre Python et ses bibliothèques natives.

## Première livraison

L’[atelier installable](../apps/freecad/README.md) propose une pièce rigide avec
un pivot explicitement défini dans le repère global. Un exemple PartDesign
permet de changer une dimension dans FreeCAD, capturer la géométrie, recalculer
les propriétés physiques et lire les échantillons natifs sur une copie de résultat.
Les fichiers de calcul sont conservés et référencés depuis le document FreeCAD.

Le travail de cette étape porte également sur l’annulation, la fermeture des
documents, le refus des résultats devenus obsolètes et la sauvegarde/réouverture.
La qualification doit porter sur l’archive réellement installée, un moteur
externe réel et des références physiques indépendantes.

## Étapes suivantes

- Représenter les corps et les liaisons d’un mécanisme dans le document FreeCAD,
  avec leurs identités, unités et repères explicites. Qualifier d’abord un
  assemblage à deux corps et deux pivots contre une référence indépendante.
- Définir la politique de conservation des références aux faces et arêtes après
  modification de la CAO. Une référence perdue doit être signalée et réparée.
- Porter les parcours de maillage, conditions aux limites, CalculiX et Pinocchio
  dans les mécanismes d’extension de FreeCAD, en conservant les contrats
  d’admission et les archives de calcul existants.
- Simplifier l’installation du moteur et la distribution sur les plateformes
  ciblées ; rendre le document et ses calculs transférables ensemble.
- Extraire progressivement les services partagés encore situés dans le paquet
  `vinkulum_studio`, à mesure que les parcours FreeCAD les utilisent.

Le dépôt FreeCAD amont reste la base de l’interface. Un fork du cœur devra être
motivé par une limitation concrète des points d’extension, avec un exemple
reproductible et un coût de maintenance identifié. La première intégration
utilise FreeCAD standard.
