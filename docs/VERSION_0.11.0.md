# Vinkulum 0.11.0 — réponses groupées avec contrôles physiques

Six charges à une même fréquence répétaient la construction du contrôle,
la factorisation du Schur et les relèvements. `ReductionMaterielle.reponses`
partage maintenant ce travail et rend un dictionnaire par colonne, avec
champ physique, normes et bornes. La nouvelle méthode justifie une version
mineure ; les signatures et les sorties de `reponse` sont conservées.

Les appels successifs à `reponse` partagent eux aussi la dernière fréquence.
Le cache est borné à un état de fréquence et invalidé après enrichissement.
Deux appels simultanés sur un modèle inchangé conservent leur propre état,
même si l'autre appel remplace le cache entre-temps. Toute mutation du
modèle doit être synchronisée par l'appelant.

Les normes sont évaluées sur les champs reconstruits : un Gram réduit
peut annuler une petite combinaison non nulle par arrondi. Les tableaux
retournés restent indépendants. Un `ControleChamp` d'assemblage périmé
après enrichissement ou changement d'applications, masse ou raideur externe
est explicitement refusé et doit être reconstruit.

## Mesures et portée

La [campagne dédiée](REPONSES_GROUPEES_0.11.0.md) conserve 216 essais,
neuf cas et les mêmes références physiques que la confrontation 0.10.0.
À 128 poutres jusqu'à 20 Hz, préparation et 257 × 6 réponses avec bornes
passent de 2,029 à 0,358 s : **5,68×**. Le témoin HCB Exudyn retenu prend
1,619 s, la LU corrigée 0,220 s. À 512 poutres sur cette bande, le total
passe de 2,978 à 1,046 s, dont 0,821 s de préparation.

Les six configurations Vinkulum acceptées conservent une erreur relative
d'opérateur inférieure à 5,56 × 10⁻⁹ dans les normes testées. Les trois
bandes à 40 Hz restent refusées. Les bornes de champ restent calculées en
doubles et ne constituent pas une certification machine. La minoration
spectrale automatique conserve sa certification des données binary64.
Ni déflation du complément spectral, ni extension aux tangentes non
linéaires, contacts ou intégration temporelle n'est livrée dans ce lot.

## Validation et livraison

Huit régressions couvrent les masses couplées, charges mixtes et nulles,
l'oracle complet indépendant, les petites déformations annulantes,
l'invalidation, les entrées invalides, l'indépendance des résultats et
l'entrelacement concurrent déterministe. Onze contre-épreuves contrôlent
la provenance et la qualification des mesures archivées.

La CI étendue `ci/local.sh --bancs`, les tests de la roue installée hors
dépôt et les exemples documentés sont consignés dans le
[rapport de livraison](bancs/version-0.11.0.json). Le hook `pre-push`
rejoue la CI locale sur le commit livré. Les environnements chronométrés
restent figés ; une roue de livraison séparée incorpore le README final,
avec vérification de l'identité de tous les modules Python et de l'extension
native par rapport à la roue 0.11.0 mesurée.

Le [guide](REDUCTION_PORTS.md) donne les formes, unités, refus et règles
de concurrence de la nouvelle méthode. Aucun gain universel sur MBDyn,
Exudyn ou Simpack n'est déduit de ces cas matériels linéaires.
