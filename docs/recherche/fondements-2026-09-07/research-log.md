# Recherche fondamentale — 7 septembre 2026

## Périmètre

Approfondissement de la recherche Deep Research demandée : identifier les résultats
fondamentaux de 2001 à septembre 2026 qui peuvent changer le coût à erreur physique
imposée d'un solveur multicorps généraliste. Référence livrée : Vinkulum 0.9.0,
commit 7850f6d. Le prototype modal présent dans l'arbre de travail reste expérimental.
La précédente étude de 58 références et la pièce jointe « Ideal Theoretical Stack
for Multibody Dynamics_ Brick-by-Brick Mathematical State of the Art and
Composability.md » sont des points de départ, pas des preuves d'exclusivité.

Plan : découverte ciblée ; vérification des théorèmes et contre-indications ;
synthèse d'un programme mathématique falsifiable ; production et contrôle PDF.
`update_plan` est absent du catalogue et son appel échoue ; plan maintenu ici.
La délégation a échoué sur la limite de fils agents ; recherche directe.
Aucun code d'implémentation de solveur concurrent consulté.

## Matrice de lacunes après la première fusion

| Question | Évidence acquise | Lacune décisive / prochaine lecture |
|---|---|---|
| Peut-on obtenir un coût quasi linéaire ? | Owhadi 2017 et Owhadi–Zhang 2017 : localisation et bornes pour opérateurs coercifs | Relire hypothèses ; le KKT non linéaire général n'est pas coercif |
| Que faut-il réduire aux interfaces ? | Smetana–Patera 2016 : opérateur compact de transfert, largeur optimale | Noyaux rigides, charges internes, constantes, domaine paramétrique ; ne pas confondre statique et dynamique |
| Peut-on supprimer la dépendance aux raideurs ? | FEEC 2010 et préconditionnement opérateur 2011 | Dériver le complexe et l'inf-sup propres aux poutres/liaisons ; aucun transfert automatique |
| Une variété non linéaire suffit-elle ? | Cohen et al. 2022 : largeurs stables et entropie ; Pagliantini 2021 : base symplectique évolutive | Stabilité encodeur/décodeur, rang évolutif et coût d'évaluation du champ |
| Qu'est-ce qui manque après réduction ? | GENERIC 2025 : compression du bain et fermeture markovienne sous hypothèses | Lire théorème et conditions ; distinguer bain infini de structure flexible finie |
| Contact et convergence rapide ? | SCD semismooth* : théorème 4.4, convergence locale sous régularité | La source était déjà recensée ; approfondir la condition et les gradients, aucune garantie globale |
| Coût des boucles ? | LCABA/proxBBO 2025 : élimination anticipée, boucles locales vs couplées | Relire coût par passe, itérations proximales et singularités |
| Différenciation concurrentielle ? | MBDyn CMS, Exudyn arbre/FFRF, Simpack réduction non linéaire déjà publics | Lire guides officiels précis ; absence de détail public ne prouve pas absence d'algorithme |

## Vérifications ciblées et arrêt

Les résultats structurants ont été relus : théorème 1 de Fürer et al. 2025
(opérations de corps, décomposition fournie) ; FEEC 3.8–3.9 ; gamblets 5.5–5.6
et extension temporelle 3.1 ; ports 3.5/3.8 ; réduction randomisée 3.2/3.7 ;
largeurs stables 3.3 ; lente–rapide A1–A3 et théorème 1 ; GENERIC hypothèse 1.2
et théorème 4.10 ; SCD 4.4. Les normes et les domaines d'application sont
explicités. La borne résidu/erreur et l'élimination avec mémoire sont
redérivées dans le rapport, sans revendication de nouveauté.

La vérification d'antériorité LCABA a révélé une version TRO 2026 et un
manuscrit Delassus 2026. Les textes auteur ont été lus ; le contraste entre
produit linéaire et inverse amorti à coût de pire cas supérieur est conservé.
Une ultime lecture de V-C révèle une divergence mémoire : O(n+m²) dans
l'introduction, O(n+dm²) dans la section détaillée. Elle est explicitement
signalée dans le rapport ; aucune borne mémoire générale n'est retenue.
La formule de LCABA est reliée à son paramètre de voisinage,
distinct de la largeur arborescente du résultat d'algèbre exacte.

La page KinematicTree échoue deux fois à l'ouverture ; son contenu officiel
est visible dans le résultat de recherche. Le DOI d'une version publiée du
manuscrit Delassus n'est pas accessible ; le rapport cite explicitement le
manuscrit auteur. Aucune absence chez Simpack n'est inférée de sa documentation.

Arrêt des recherches : chaque décision structurante dispose d'un résultat
primaire avec conditions ou d'une lacune explicite. Les suites susceptibles
de modifier les décisions sont désormais les dérivations spécifiques et les
expériences Vinkulum. Ce n'est pas une revue exhaustive mondiale.

## Production et contrôle terminés

Découverte, vérification, synthèse et production terminées. `update_plan` a
été réessayé, toujours indisponible. La source canonique produit un PDF de
14 pages avec ReportLab 5.0.1, DejaVu et le moteur de rendu du dépôt :

`python ci/rend_recherche.py --source docs/recherche/fondements-2026-09-07/report-source.md --sortie docs/FONDEMENTS_MATHEMATIQUES_VINKULUM_2026.pdf --titre 'Vinkulum : les fondements d’une rupture de performance'`

Les 183 blocs ou cellules de la source sont retrouvés dans le texte extrait.
Les 47 annotations de lien couvrent les 21 sources. La géométrie des mots
et des liens est contrôlée sur toutes les pages. Les 14 pages ont été
examinées visuellement à 100 dpi ; la définition d'une norme a été ajoutée
page 10, puis cette page a été rendue et relue. La divergence bibliographique
Delassus a ensuite été explicitée page 3, rendue et relue. Aucun débordement ou glyphe
manquant observé. Les empreintes et contrôles figurent dans `verification.json`.

La livraison concerne la recherche, sa traçabilité, les priorités des
documents du projet et une option de titre des métadonnées du moteur PDF.
La version logicielle reste 0.9.0. Le prototype modal n'est pas inclus.
