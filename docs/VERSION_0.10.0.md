# Vinkulum 0.10.0 — réduction matérielle et minoration spectrale vérifiée

Cette version raccorde les facteurs d'énergie du noyau à une réduction
par interfaces. Un modèle de poutres construit avec `Noyau` peut fournir
directement son modèle matériel `K=DᵀD`, sa masse et ses contraintes ;
la réduction calcule une réponse harmonique réelle, reconstruit le champ
et enrichit sa base aux fréquences choisies. Le facteur statique et son
certificat spectral sont réutilisés pendant l'enrichissement.

`Noyau.facteurs_materiels_poutres()` et le module
`vinkulum.reduction_ports` étendent l'API sans modifier les appels existants.
Cette extension compatible justifie une version mineure. L'option
d'installation `reduction` ajoute SciPy. L'[exemple et le contrat public](REDUCTION_PORTS.md)
détaillent les indices physiques, la métrique des ports et les forces.

La constante spectrale intérieure peut être obtenue automatiquement sans
calcul de modes : la récurrence classique d'inverse sélectionnée évalue
une trace, et des intervalles à arrondis dirigés encadrent cette trace
ainsi que l'écart entre D d'entrée et sa factorisation QR calculée.
Le [carnet de preuves](INVERSE_SELECTIONNEE_PREUVES.md) dérive la fermeture
des dépendances, la borne d'ordre matriciel et la minoration finale.
La récurrence sélectionnée et les algorithmes de contrôle du champ sont
embarqués dans la roue ; les anciens imports `ci/` restent des adaptateurs.

Les contre-épreuves comprennent des calculs rationnels exacts, des
arrondis Decimal ambiants défavorables, des inerties couplées, des
permutations et des rotations de repère. Elles ont aussi conduit au
refus des coefficients creux dupliqués avant conversion et à la
reconstruction des dépendances après modification d'un objet d'inverse.
Les tests d'adaptation comparent les déplacements à un oracle dense
et vérifient la conservation du facteur et du certificat.

La portée de chaque garantie est explicite :

- `certification_machine=True` concerne uniquement le minorant spectral
  du modèle défini par D et M en binary64. La borne de trace peut réduire
  la bande utilisable ; son coût dépend du remplissage du graphe.
- Les contrôles de réponse et de champ restent évalués en doubles,
  avec `certification_machine=False`. Une tolérance relative acceptée
  porte sur les fréquences demandées. Près d'une résonance globale,
  l'algorithme peut refuser de conclure.
- L'extraction désigne la partie matérielle des poutres. Une précontrainte,
  les autres forces et les effets dynamiques ne sont pas absorbés dans
  DᵀD. Les trajectoires non linéaires réduites restent à développer.
- La voie native exige des contraintes fixes holonomes éliminables par
  suppression de coordonnées. La certification automatique de masse
  intérieure couvre des blocs connexes SPD de taille au plus six.

La CI étendue passe **83 tests Rust**, **8 tests du prototype de poutre
mixte**, **107 tests Python publics**, **41 groupes de vérification**,
**46 bancs mécaniques**, **9 bancs de contact** et **189 entrées d'API**.
Les suites spécialisées comprennent les **12 tests du minorant spectral**.
Les [contrôles de livraison et empreintes](bancs/version-0.10.0.json)
consignent aussi l'installation de la roue hors du dépôt, le rejeu des
107 tests publics et des 12 tests spectraux, ainsi que les exemples.
Le [témoin natif reproductible](bancs/reduction-native-0.10.0.json)
décrit ses tailles, paramètres, précision et coûts de préparation.
Sur une console d'aluminium de 1 m, section 40 × 6 mm, sollicitée par
0,1 N en flexion faible, l'adaptation passe de 6 à 12 directions
intérieures aux fréquences 0, 1 et 2 Hz ; les six ports sont conservés.
Les maxima ci-dessous portent sur ces trois points, pour la norme
de déformation :

| Poutres | DDL intérieurs avant réduction | Directions finales | Borne relative évaluée | Écart relatif à l'oracle |
|---|---:|---:|---:|---:|
| 32 | 186 | 12 | 1,23e−9 | 1,85e−11 |
| 128 | 762 | 12 | 7,65e−9 | 1,87e−11 |
| 512 | 3 066 | 12 | 1,03e−7 | 2,09e−11 |

L'oracle dynamique assemble ses blocs de Timoshenko à partir des paramètres
physiques et résout la chaîne en Decimal à 70 et 90 chiffres. Sa limite
statique est aussi confrontée à la formule analytique de console. Cet
oracle contrôle la précision observée ; il ne certifie pas les arrondis
du champ calculé. Le seuil demandé est 1e−6 pour les deux normes physiques.
Le facteur et le certificat sont conservés pendant l'enrichissement.

Pour reproduire avec les sources et l'environnement 0.10.0 :

```bash
python ci/experience_reduction_native.py --n 32 128 512 --output /tmp/reduction-native.json
python ci/experience_reduction_native.py --verifier --output /tmp/reduction-native.json
```

Le vérificateur exige les mêmes empreintes de code et de binaire que le
lot contrôlé, puis recalcule ses oracles. Les coûts archivés sont ceux
d'un seul passage, sans protocole de mesure isolée, et ne permettent
pas un classement de vitesse.
Aucune nouvelle comparaison de performances aux implémentations de
réduction d'Exudyn, MBDyn ou Simpack n'est revendiquée.
