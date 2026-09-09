# Vinkulum 0.12.0 — complément contraint dans le paquet installable

La réduction contrainte, le contrôle des masses couplées et l'assemblage
des compléments étaient accessibles par les scripts de recherche. Ils sont
maintenant disponibles dans une roue avec
`vinkulum.reduction_contrainte.ReductionContrainte` et `AssemblageContraint`.
Cette nouvelle API justifie une mineure ; les API 0.11 restent inchangées.

L'extension embarque les certificats C++ d'intervalles derrière une frontière
Rust validant les entrées. Aucun compilateur n'est nécessaire au premier
appel de la roue. La construction depuis les sources exige désormais C++17.
La masse complète est contrôlée, y compris ses couplages aux ports et les
coordonnées explicitement sans masse.

La sélection des directions utilise la masse intérieure complète ; le
complément réellement obtenu est vérifié par inertie au seuil demandé.
Cela permet de franchir certaines résonances intérieures en conservant
leurs directions dans le système réduit. Les marges d'assemblage portent
sur le bloc global, ce qui permet aussi un Schur local singulier lorsque
le système assemblé reste régulier.

Les certificats de masse et de complément sont distincts des **bornes de
réponse numériques, non certifiées par machine**. Ce lot ne fournit pas de
nouvel intégrateur temporel, de contact, d'amortissement, ni de preuve d'une
supériorité générale sur Exudyn, MBDyn ou Simpack.

[Contrat complet et exemples](REDUCTION_CONTRAINTE.md) ·
[Preuves du transport massique](MASSE_COMPAREE_PREUVES.md) ·
[Preuves d'assemblage](ASSEMBLAGE_COMPLEMENTS_PREUVES.md)

## Validation

Les tests du paquet confrontent les champs à une résolution rationnelle du
système physique complet, vérifient les pôles locaux/globaux, la masse
complète, les copies, les refus et le pont natif. Les six preuves archivées
(masse et complément, trois tailles) sont reproduites exactement par
l'extension. Un calcul à directions identiques reproduit aussi les champs,
normes et majorants du prototype bit à bit.

La qualification de livraison compare séparément les sélections automatiques
de l'API à des références Decimal 70/90 chiffres déjà convergées : trois
pièces seules et trois assemblages à trois branches, 257 pulsations de 0 à
40 Hz et six charges. C'est un contrôle de précision, pas une campagne de
classement de vitesse. Les performances historiques restent attribuées
aux prototypes et aux versions qui ont effectivement été mesurés.

Le bilan de livraison et les journaux sont conservés dans
[`bancs/version-0.12.0.json`](bancs/version-0.12.0.json).

Les six cas sont admis au seuil 10⁻⁶, pour les champs et leurs majorants.
Les premières tentatives restent archivées : environnement non installé
après un échec réseau, puis refus des métriques historiques légèrement
dissymétriques. La qualification explicite ensuite leur triangle inférieur,
comme le Cholesky historique ; D, M et les références restent inchangés.
La revue a également corrigé le tri interne après permutation des DDL,
avec une contre-épreuve à deux directions retenues.

Les grands champs NPY sont identifiés mais exclus de l'archive compacte.
La CI en vérifie les inventaires et requalifie les diagnostics ; elle ne
prétend pas rejuger ces champs absents. Les journaux et les empreintes de
la roue qualifiée sont conservés avec le bilan de version.
