# Vinkulum 0.13.0 — unicité et bornes vérifiables de l'erreur linéaire

Cette version ajoute une certification facultative des systèmes linéaires
et des accélérations initiales calculées par le noyau. Elle fournit un
certificat transportable qui établit **l'unicité de la solution et une borne
de son erreur**, sur les coefficients numériques exportés.

Les trois nouvelles fonctions publiques de `vinkulum.certification` sont
`certifier_systeme`, `certifier_initialisation` et `verifier_certificat`.
Cette nouvelle capacité justifie le passage à **0.13.0**. L'appel ordinaire
au simulateur conserve ses critères ; la certification dense est demandée
explicitement et n'est pas exécutée automatiquement à chaque pas.

## Garantie obtenue

Pour A, b, x et une inverse approchée R, les coefficients binary64 sont
interprétés exactement. Le certificat vérifie en arithmétique exacte :

    η = ||I − RA||∞ < 1
    ||x* − x||∞ ≤ ||R(b − Ax)||∞ / (1 − η)

La première inégalité établit l'inversibilité de A. Les majorants publiés
sont arrondis vers le haut et contrôlés en rationnels, avec une borne
par composante. Une exigence `erreur_max` est refusée lorsque la borne ne
la satisfait pas. Une inverse numérique incorrecte ne peut pas, à elle
seule, provoquer une acceptation.

Le [contrat et sa démonstration](CERTIFICATION_LINEAIRE.md) détaillent
les hypothèses, les arrondis et la preuve par contraction. Le vérificateur
autonome recalcule les inégalités par un second chemin rationnel, sans
Rust, NumPy ou Vinkulum :

```bash
python3 -I -S python/vinkulum/_verification_lineaire.py docs/bancs/certificat-pendule-0.13.0.json
```

La démonstration n'est pas une preuve formelle du code dans un assistant
de preuve. Les opérations entières et rationnelles, les conversions,
les compilateurs, le matériel et la fidélité de l'export constituent
encore une base de confiance explicitée.

## Deux exemples qui distinguent résidu et précision

Pour le [pendule documenté](bancs/certificat-pendule-0.13.0.json), la
borne globale vaut environ **1,387778780781446 × 10⁻¹⁶**, avec des bornes
par accélération et multiplicateur. Cette valeur concerne l'erreur de
résolution du système initial exporté ; elle ne borne pas la trajectoire.
Les coordonnées conservent leurs unités propres : leur maximum numérique
n'est pas une norme physique invariante aux changements d'unités.

Pour [A = diag(1, 2⁻⁶⁰), b = (1, 2⁻⁶⁰), x = (1, 0)](bancs/certificat-residu-trompeur-0.13.0.json),
le résidu est inférieur à 10⁻¹² alors que l'erreur exacte vaut **1**.
La certification donne cette borne et refuse une précision demandée de
10⁻¹². Un résidu minuscule ne suffit donc pas à établir une solution précise.

## Rattachement au noyau et périmètre

Un témoin facultatif capture les contributions, le second membre et la
solution réellement produits par le chemin natif `acc_init`. Le calcul
est exécuté sur une copie du modèle. NumPy propose l'inverse approchée ;
il ne remplace pas les accélérations natives dans ce certificat.

Les contributions de matrice sont additionnées exactement. Si leur somme
ne peut pas être représentée exactement en binary64, l'export est refusé.
Le vérificateur contrôle aussi cette somme. L'état du modèle reste intact,
y compris lors d'un échec.

Le budget par défaut est de 64 inconnues, avec une limite de 128.
Contraintes désactivées ou redondantes, contacts non lisses, partitions
gelées et schémas de redémarrage imposés sont actuellement exclus.
Un refus de contraction ne démontre pas à lui seul une singularité.

**Le noyau entier n'est pas certifié.** Les erreurs des coefficients
mécaniques, la géométrie, les quotients de contraintes, les événements
et l'erreur globale de temps restent à traiter. Le dossier exécutable
maintient neuf obligations ouvertes ; `--exiger-couverture-complete`
retourne le code **2** attendu.

Le garde exact de vitesse de la 0.12.2 reste en place. Sa documentation
de portée est précisée : il intervient aussi lors des réassemblages
déclenchés par un changement du masque des contraintes redondantes.
Les autres projections temporelles ne sont pas couvertes par ce garde.

## Vérification et traçabilité

- 15 nouveaux tests : référence par élimination rationnelle indépendante,
  systèmes indéfinis, permutations, échelles extrêmes, mauvais
  conditionnement, sous-normaux, refus et certificats altérés.
- Une borne abaissée d'un seul ulp sous 1/3 est refusée. La contraction
  `1 − 2⁻¹⁰⁷⁴ < 1` reste décidable exactement, malgré son arrondi à 1
  en binary64. La solution native est comparée aux accélérations et
  réactions enregistrées au départ du premier pas de simulation.
- 12 certificats de campagne relus indépendamment, dont six pendules
  de longueurs 10⁻³, 1 et 10³ ; critères fixés avant exécution.
- 13 314 contre-épreuves du garde entier, aucune divergence avec l'oracle,
  et sept entrées invalides refusées.
- 160 tests Python sur une roue installée dans un environnement neuf,
  87 tests Rust et huit tests du prototype de poutre.
- 41 cas de vérification, 46 bancs étendus et neuf cas de contact.
  Trois intégrations Exudyn restent ignorées, faute de binaire installé.
- 2 304 projections et six trajectoires de pendule : résultats identiques
  à ceux de la roue 0.12.2, avec les mêmes critères. Les deux explorations
  de conditionnement extrême restent des refus avec restitution de l'état.

Le [manifeste](bancs/version-0.13.0.json) et l'[archive des preuves](bancs/version-0.13.0-preuves.json.gz)
conservent les sources identifiées, la roue, les certificats complets,
les journaux et les limites de la livraison. Les deux exemples JSON ont
également été relus avec `python3 -I -S`, sans importer de paquet tiers.

Roue validée : `vinkulum-0.13.0-cp314-cp314-linux_x86_64.whl`, SHA-256
`58208fccdaa6df777f2b71958021f218f3a66321a6a1db5a37eb445e6651e8a3`.
Environnement : CPython 3.14.7, NumPy 2.5.3, SciPy 1.18.1, Rust 1.98.0,
Linux x86_64. Cette livraison ne qualifie pas les autres plateformes.

L'inverse et la vérification denses ont un coût cubique. Les temps
enregistrés par la campagne servent au diagnostic et ne constituent
pas une comparaison de performances entre solveurs. La méthode ne
revendique aucun gain de vitesse ni une certification des grands systèmes.
