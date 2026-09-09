# Refonte de Vinkulum Studio — septembre 2026

Objectif : une mise à niveau majeure de l'interface de travail, depuis l'édition
d'un mécanisme jusqu'à l'examen de résultats traçables. La qualification porte
sur des parcours réels dans l'application et dans son paquet macOS, pas sur une
maquette ni sur une revendication de parité universelle avec les outils CAO.

## Exigences et preuves attendues

| ID | Livraison attendue | Preuve de recette |
|---|---|---|
| UX-01 | Espace 3D dégagé, panneaux réorganisables et restaurables, menus natifs, commandes regroupées | Captures conception/résultats, fenêtres 1280×800 et 1440×950, restauration après fermeture |
| UX-02 | Palette de commandes filtrable, raccourcis natifs et accès clavier aux fonctions essentielles | Tests clavier réels ; commandes indisponibles non exécutables |
| UX-03 | Explorateur filtrable, sélection cohérente arbre/scène/propriétés, visibilité et isolation réversibles | Recherche par nom/type ; objets masqués conservés dans le modèle et le calcul |
| UX-04 | Inspecteur structuré, unités et repères explicites, modifications en attente visibles, validation sans perte des champs | Recette saisie invalide, changement de sélection, annulation et résultat reçu pendant la saisie |
| UX-05 | Pilotage 3D précis : orientation, cadrage, projection, modes de manipulation et informations de navigation | Essais VTK/Qt, manipulation numérique et souris, rendu macOS et Linux |
| UX-06 | Espace résultats : courbes interactives, tableau des échantillons, lecture contrôlée, comparaison identifiée | Sélection temporelle synchronisée scène/courbe/table ; modèles et unités des calculs visibles |
| UX-07 | États intelligibles du document et du calcul, diagnostics navigables, provenance et erreurs accessibles | Calcul réussi, modèle invalide, erreur, annulation et résultat antérieur conservé |
| UX-08 | Thèmes contrastés, mise à l'échelle, focus visible, noms accessibles, disposition persistante | Captures à 100 % et 200 %, parcours clavier ; limites des essais d'accessibilité documentées |
| UX-09 | Exécution réactive, budgets mémoire conservés, fichiers et résultats intègres | Tests existants et nouveaux, mesures de rendu et de grandes courbes |
| UX-10 | Version publique, documentation à jour, DMG Apple Silicon réellement exécuté | CI macOS, rapport du bundle, empreinte du DMG et release publique |

Les exigences GUI-01 à GUI-08 et UI-01 à UI-06 du cahier des charges v1.1 restent
applicables. La refonte ne change pas les garanties scientifiques du noyau.
Un modèle dont les contrôles d'entrée passent n'est pas une trajectoire certifiée.
Les fonctions futures de CAO, de collaboration et de multiphysique du cahier des
charges ne deviennent pas disponibles par un changement de présentation.

## Références consultées le 9 septembre 2026

- [FreeCAD 1.1, publié le 24 mars 2026](https://freecad.github.io/Website/download/releases/1-1/) : transformation précise, navigation, recherche, thèmes et retours visuels. Ces directions motivent UX-03 à UX-05 ; elles ne prouvent pas leur implémentation dans Studio.
- [ParaView 6, personnalisation](https://docs.paraview.org/en/v6.0.0/ReferenceManual/customizingParaView.html) : recherche des propriétés, séparation simple/avancé, préférences persistantes et restauration. Application à UX-01 et UX-04.
- [Qt 6.11, accessibilité](https://doc.qt.io/qt-6/accessible.html) : navigation clavier, adaptation de la taille, contrastes et sémantique des composants. Application à UX-02 et UX-08.

## État de réalisation

Travail en cours. La version 0.2.0 publiée constitue la référence avant refonte.
Chaque exigence devra être confrontée aux preuves effectives avant de déclarer
la mise à niveau terminée.

## Direction précisée par l’utilisateur

La première présentation a été jugée « cheap et pas sérieuse ». Une simple
recoloration de formulaires ne satisfait donc pas la cible. La scène doit dominer,
l’inspecteur doit présenter des valeurs lisibles, les commandes doivent avoir une
iconographie cohérente et la densité doit s’adapter à l’atelier actif. Les captures
réelles servent à éliminer les débordements et les incohérences de hiérarchie.

Les références demandées orientent des choix précis. Elles ne sont ni des
bibliothèques intégrées, ni une promesse de reproduire toutes leurs capacités :

| Référence primaire | Principe retenu pour Studio | Traduction dans la refonte |
|---|---|---|
| [Plasticity : interface](https://doc.plasticity.xyz/plasticity-essentials/plasticity-interface/user-interface-overview), [palette](https://doc.plasticity.xyz/plasticity-essentials/plasticity-interface/command-palette) | Scène dominante, outils compacts, accès direct aux commandes | Chrome graphite, icônes vectorielles, palette recherchable |
| [Blender : outils et ateliers](https://docs.blender.org/manual/en/4.5/interface/tool_system.html) | Un outil actif et un contexte de travail explicite | Ateliers Modéliser / Simuler / Examiner, sélection / déplacement / rotation |
| [Shapr3D : interface adaptative](https://support.shapr3d.com/hc/en-us/articles/7873882619548-Adaptive-user-interface) | Actions pertinentes pour la sélection | Outils contextuels et commandes indisponibles explicitement désactivées |
| [Abaqus/CAE : gestion et visualisation](https://www.3ds.com/fileadmin/Products/Simulia/PDF/datasheets/Abaqus_CAE_Datasheet.pdf) | Séparer modèle, exécution et examen des résultats | Instantanés en lecture seule, historique des calculs, provenance |
| [Onshape : recherche d’outils](https://cad.onshape.com/help/Content/Home/search_tools.htm) | Une commande reste découvrable sans connaître sa position | Recherche par mots, raccourcis visibles et activation clavier |
| [NX : accès aux commandes](https://blogs.sw.siemens.com/designcenter/designcenter-x-nx-tips-and-tricks-copilot/) | Aider à trouver l’opération pertinente dans un outil riche | Registre commun de commandes et contexte ; aucun assistant IA fictif |
| [STAR-CCM+ 2606](https://blogs.sw.siemens.com/simcenter/simcenter-star-ccm-2606-released/) | Faciliter l’examen des différences entre simulations | Superposition de séries, identité des objets et différences de réglages |
| [SolidWorks : raccourcis et menus contextuels](https://blogs.solidworks.com/products/solidworks/useful-keyboard-shortcuts-workflow-customizations-solidworks/) | Réduire le trajet jusqu’aux outils courants | Palette contextuelle S, menus natifs et outils près de la scène |
| [3DEXPERIENCE : arbre, zone 3D et barre d’action](https://3dswym.3dexperience.3ds.com/wiki/solidworks-news-info/getting-started-with-3dexperience-simulation-solidpractices_rFZtKhrBSdO2cIlYITOksg) | Lier sélection et contexte d’action | Arbre, scène et inspecteur partagent les identités stables |
| [Autodesk Fusion : interface](https://help.autodesk.com/view/fusion360/ENU/?contextId=LP-STEPS-P13N-SNP-GS-OTH-CRD-1) | Ateliers et navigation spatiale directement accessibles | Ateliers persistants, orientation interactive de caméra, isolation |
| [Rhino : Gumball](https://www.rhino3d.com/en/docs/guides/user-guide/gumball-basics/), [dispositions](https://www.rhino3d.com/features/user-interface/window-layouts/) | Concilier manipulation et précision numérique | Modes du manipulateur, composantes X/Y/Z, panneaux restaurables |
| [Creo : recherche](https://support.ptc.com/help/creo/creo_pma/r12/usascii/fundamentals/fundamentals/to_search_a_command.html) | Retrouver une commande par son nom et son aide | Recherche sans distinction d’accents, actions issues du même registre |

Les principes plus anciens toujours utiles ne sont pas présentés comme des
inventions de 2026. Les pages consultées sont un état de documentation, pas une
évaluation exhaustive de licences commerciales exécutées localement.

## Première itération implémentée

- Espace principal Qt/VTK réorganisé, menus natifs, ateliers et panneaux.
- Palette globale Ctrl/Cmd+K et palette contextuelle S ; filtres de l’explorateur.
- Thèmes graphite et clair, icônes originales vectorielles et focus visible.
- Champs vectoriels par composante. L’affichage compact conserve chaque valeur
  binaire originale tant que sa composante n’est pas modifiée.
- Orientation de caméra interactive, sélection / translation / rotation,
  masquage et isolation sans changement du document ou du calcul.
- Historique de session : huit calculs et 128 Mio de tableaux au maximum. Les plus
  anciens sont évincés lorsque le budget est atteint. Les exports restent explicites.
- Courbes avec sélection temporelle, zoom, déplacement et lecture à vitesse réglable ;
  tableau paresseux des échantillons, sans copie de chaque cellule.
- Comparaison par UUID et grandeur physique, avec chaque série sur ses temps natifs.
  Les noms peuvent changer ; les identités différentes ne sont pas appariées par
  position. Aucun rééchantillonnage ni calcul d’écart interpolé n’est implicite.
- Sauvegarde proposée avant abandon, document modifié signalé, champs invalides
  conservés et fin de calcul incapable d’écraser une saisie en attente.

Tests Linux réalisés : 38 tests, comprenant les 29 tests antérieurs et neuf tests
supplémentaires. Les tests de GUI surveillent aussi les exceptions Qt différées,
car un résumé unittest « OK » seul ne suffit pas à les détecter. Un défaut de
nettoyage du wrapper VTK a été corrigé après sa détection dans les journaux.

Restent à qualifier avant clôture : ergonomie finale des outils contextuels,
accessibilité et écrans à forte densité, performances mesurées, robustesse de la
restauration des panneaux sur plusieurs écrans, exécution du nouveau bundle macOS
et livraison publique de la nouvelle version. Les fonctions et preuves ci-dessus
ne déclarent pas le grand objectif achevé.

Recette à 200 % : [rapport Linux](bancs/studio-gui-2026/linux-hidpi.json),
[modélisation](bancs/studio-gui-2026/modeling-dark.png),
[comparaison](bancs/studio-gui-2026/comparison-dark.png),
[écran portable](bancs/studio-gui-2026/comparison-light-1280.png).
La capture mesure 2880×1900 pixels pour une fenêtre logique 1440×950.
Avec trois corps et deux courbes, le rendu logiciel de cet environnement mesure
36,2 ms par image en médiane et 56,1 ms au 95e centile. Ce n’est ni une mesure GPU
sur le Mac utilisateur ni une garantie de fréquence. La recette vérifie aussi la
présence des commandes de lecture dans une fenêtre logique de 1280×800.
