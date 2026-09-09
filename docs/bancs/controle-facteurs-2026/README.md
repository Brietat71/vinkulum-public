# Contrôle par facteurs — campagne du 8 septembre 2026

Archive compacte de **60 processus terminés** : trois maillages natifs
(32, 128 et 512 éléments), cinq variantes, une chauffe et trois mesures
par configuration. Chaque réponse couvre 257 fréquences de 0 à 40 Hz et
six charges. Tous les résultats et refus de précision sont conservés.
Aucun résultat défavorable n'a été remplacé par une relance.

Les variantes sont `facteurs_controle`, `inertie_controle`, `lu_corrigee`,
`hcb_standard` et `hcb_energie`. Le nouveau contrôle remplace des évaluations
de normes d'opérateur par un majorant du Gram déjà calculé et des facteurs
minces de réparation audités. Les résolutions et normes physiques gardent
leur chemin précédent. Les SHA des champs des deux candidats sont égaux
pour les **12 paires maillage/passage**. Les erreurs, bornes et diagnostics
de facteurs restent conservés séparément pour chaque essai.

## Qualifications

Le juge impose `1e-6` aux erreurs des colonnes et opérateurs dans les
normes masse, déformation et port. **50/60 essais** passent ce juge :
37/45 mesures et 13/15 chauffes. **12/15 configurations** passent leurs
quatre essais. Les 24 essais candidats et les 12 essais LU corrigée passent.

HCB standard à 128 éléments échoue sur deux des trois mesures ; sa chauffe
et une mesure passent. Les deux HCB à 512 éléments échouent sur les quatre
passages. Leurs résultats figurent dans le [bilan complet](bilan.json), mais
ces trois configurations ne produisent aucun ratio de performance admis.

Les **24 contrôles** passent aussi les marges globales, les majorants
relatifs et la confrontation des majorants avec les oracles. Le contrôleur
par facteurs atteint les mêmes erreurs physiques que le contrôle précédent.
Le lot ne comprend pas de nouvelle mesure du certificat par trace.

Les médianes totales du nouveau contrôle sont 0,154800 s, 0,253325 s et
0,671176 s aux trois tailles ; les médianes du contrôle précédent sont
0,171160 s, 0,298018 s et 0,819090 s. À 512 éléments, LU corrigée vaut
0,695890 s : l'écart total observé est faible. Les plages de toutes les
mesures sont conservées ; trois répétitions ne démontrent pas une avance
universelle ni sa stabilité sur d'autres machines.

## Coût compté et portée

La sélection d'une direction intérieure, sa contrainte B, le certificat
frais d'inertie, quatre blocs de Krylov (128 directions au maximum), le
contrôle de profondeur 8 et ses facteurs éventuels sont inclus dans la
préparation. Les réponses comprennent la reconstruction des champs et
leurs normes physiques. Les processus sont frais, sur le CPU 8, avec un
seul fil numérique ; l'ordre des variantes est inversé à chaque passage.

Les imports, lectures, empreintes de provenance, sauvegardes et juges
contre les oracles sont hors chronomètre. Les empreintes produites par
l'appel du certificat sont incluses dans le chronomètre extérieur de
certification du pilote ; son temps interne s'arrête avant leur émission.
Les médianes de phases distinctes ne doivent pas être additionnées comme
si elles constituaient la médiane d'un même essai.

Le certificat dirigé vise `D_i.T D_i − gamma M_ii` sur `ker(B.T)`, les
entrées binary64 étant interprétées comme exactes. Le gamma est
`252661.87266788757`, soit `0x1.ed7aefb394d2bp+17`, issu de `(2*pi*80)**2`.
La preuve à 32 chiffres dirigés établit la positivité de M et l'inertie
`(n,1,0)` du KKT de dimension `n+1`. Les 24 certificats correspondent à
trois identités D/M/B distinctes. La marge `Fraction(gamma) >
Fraction(omega)**2` est vérifiée pour les fréquences NPZ.

Le seuil de 80 Hz concerne le complément contraint ; les réponses restent
qualifiées sur 0–40 Hz. Les ports et la direction conservée restent dans
le système global, dont les résonances nécessitent une marge distincte.
Le [carnet de preuves](../../CONTROLE_FACTEURS_PREUVES.md) décrit les
transferts de normes et leurs contre-exemples. Les contrôles de réponse
et les audits de facteurs sont évalués en binary64 : **aucune certification
machine des réponses** n'est revendiquée. La confrontation aux oracles
emploie la marge déclarée `64*epsilon*norme_reference`.

Les témoins HCB emploient le calcul modal officiel Exudyn et le pont de
réponse local corrigé. Il ne s'agit pas de simulations FFRF ou de systèmes
multicorps complets. Les mémoires ne sont pas comparables entre familles.
Ces consoles ne constituent pas un classement général d'Exudyn et MBDyn.

## Représentation compacte et provenance

- `rapport.json.gz` conserve les octets du rapport original, avec les
  60 résultats complets, retours par fréquence, erreurs et bornes par
  charge, refus, environnements, chronos et états terminaux.
- `provenance-fichiers.json` relie chaque fichier brut à ses SHA et taille
  d'origine. Les 60 `resultat.json` sont extraits du rapport avec
  `json.dumps(..., ensure_ascii=False, indent=2, allow_nan=False)+'\n'`.
  Leurs octets reconstruits ont été vérifiés avant suppression des doublons.
- `entrees/` conserve les trois modèles NPZ et leurs métadonnées, dont la
  qualification des références Decimal70/90 réutilisées de
  `confrontation-ports-exudyn-0.10.0`.
- `essais/` conserve B, Phi, la base physique, K/M réduits et les indices
  dans les NPZ, ainsi que les journaux, terminaux et résultats bruts des
  témoins. Les facteurs de réparation ne sont pas archivés séparément :
  leur régénération nécessite les sources figées et les modèles.
- `sources/` conserve exactement les **16 sources mesurées**. Le pilote
  est épinglé au SHA `508efe7f731cfd3c368229b62706b16e32456ae22ae881571e0ad6a3a1f30a7b`,
  le contrôleur au SHA `fed375d6b120df072d056a7029e84916dd666ab8bbfc926b3ccae849dd769bbe`.
- `manifeste-campagne-brute.json` est le manifeste original inchangé.
  `exclus.json` conserve les SHA, tailles, formes et chemins des **63 NPY
  exclus** : 60 champs et trois oracles. Les indices de fréquences valides
  des champs sont conservés. Ces fichiers restent dans `/tmp`.
- `manifest.json` couvre par SHA et taille tous les fichiers compacts,
  sauf lui-même. Le journal du parent est également conservé en gzip.

Les roues, interpréteurs et bibliothèques partagées ne sont pas copiés ;
leurs versions et empreintes figurent dans les environnements. Les profils,
pilotes exploratoires et réfutations précédant le gel sont dans
l'[archive des sondes](../controle-facteurs-sondes-2026/README.md). Leurs
chronométrages ne sont pas mélangés à ceux de cette campagne.

## Vérification et reconstitution

Depuis la racine du dépôt et avec les dépendances Python installées :

```sh
python ci/archive_controle_facteurs.py --verifier docs/bancs/controle-facteurs-2026
python ci/test_archive_controle_facteurs.py
```

Le validateur n'accède à aucun NPY exclu et n'exécute aucune source
capturée. Il recalcule la couverture, les qualifications depuis les erreurs
et bornes détaillées, les maxima, marges, rapports absolus/relatifs et
chronos. Il contrôle les diagnostics de facteurs et leurs dimensions.
Il ne reconstruit pas les erreurs physiques ni les facteurs de réparation.

Les signatures et déterminants dirigés des pivots sont vérifiés. Les
congruences sont ensuite **rejouées depuis les modèles D/M et B archivés**,
une fois par identité distincte, avec le module local connu dont le SHA
doit être identique au module mesuré. Tous les champs déterministes des
certificats doivent concorder ; seuls les temps sont exclus de l'égalité.
Les tests Fraction apportent un contrôle algorithmique indépendant sur
de petites matrices, distinct de ce rejeu du même algorithme dirigé.

Les commandes d'origine et chemins des interpréteurs figurent dans les
terminaux. Une nouvelle campagne exige des environnements correspondant
aux empreintes, un dossier de sortie absent, les modèles et les références
complètes. Les références peuvent être réutilisées ou régénérées par les
sources capturées `confronte_ports_exudyn.py`, `oracle_champs_ports.py` et
`reference_ports_precision.py`, avec vérification de leurs SHA.
Le pilote `experience_controle_facteurs.py --campagne` expose les arguments
d'entrée et d'interpréteurs ; un autre emplacement nécessite de préciser
le bilan HCB historique. Aucun chronométrage futur identique n'est promis.

Cette archive a été créée après confirmation de la terminaison des
60 processus, depuis `/tmp/vinkulum-campagne-controle-facteurs-2026`.
Les données mesurées, leurs instantanés de sources et le dossier brut
n'ont pas été modifiés pendant l'archivage.
