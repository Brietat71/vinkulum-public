# Campagne inertie contrainte — 8 septembre 2026

Archive compacte de **84 processus terminés** : trois maillages natifs
(32, 128 et 512 éléments), sept variantes, une chauffe puis trois mesures
par configuration. Chaque réponse porte sur 257 fréquences de 0 à 40 Hz et
six charges. Tous les résultats, y compris les refus de précision, sont
conservés. Aucun essai n'a été relancé pour remplacer un résultat défavorable.

Les sept variantes sont `inertie_champs`, `inertie_controle`,
`krylov_champs`, `krylov_controle`, `lu_corrigee`, `hcb_standard` et
`hcb_energie`. La sélection d'une direction intérieure, la construction
de quatre blocs de Krylov (128 directions au maximum), le certificat
frais et le contrôle éventuel de profondeur 8 sont inclus dans la
préparation chronométrée. Les imports, lectures, empreintes de provenance,
sauvegardes et confrontations avec les oracles sont hors chronomètre.
Les empreintes produites par l'appel du certificat lui-même sont incluses
dans le chronomètre extérieur de certification du pilote ; son compteur
interne `preparation_s` s'arrête avant leur émission.

## Qualifications conservées

La tolérance du juge est `1e-6`, sur les colonnes et opérateurs dans les
normes masse, déformation et port. **75/84 essais** passent ce juge :
56/63 mesures et 19/21 chauffes. **18/21 configurations** passent leurs
quatre essais.

- Les 48 essais Krylov, avec certificat d'inertie ou de trace, passent.
- Les 12 essais LU corrigée passent.
- Les deux variantes HCB à 512 éléments échouent dans leurs quatre essais.
- HCB standard à 128 éléments échoue dans une des trois mesures
  (`1.0339530569901925e-6` au pire) ; sa chauffe et les deux autres mesures
  passent. Sa configuration complète reste donc non admissible.

Les **24 essais munis du contrôle** (`inertie_controle` et
`krylov_controle`) passent également les tests de marge, de majorants
relatifs et de confrontation des majorants avec les oracles. Le
[bilan](bilan.json) contient toutes les médianes, plages, qualifications
et paramètres. La [note de travail](../../INERTIE_DIRIGEE_PREUVES.md)
décrit le fondement mathématique du nouveau certificat.

## Portée des preuves

L'inertie vise exactement `D_i.T D_i − gamma M_ii` sur `ker(B.T)`, où les
entrées binary64 de D, M et B sont interprétées comme des nombres exacts.
Le paramètre est le binary64 `252661.87266788757`, soit
`0x1.ed7aefb394d2bp+17`, obtenu par `(2*pi*80)**2`. Les congruences à
32 chiffres dirigés établissent la positivité de M et l'inertie `(n,1,0)`
du KKT de dimension `n+1`. La marge entre gamma et chaque carré exact
des fréquences archivées est vérifiée en `Fraction`.

Le seuil de 80 Hz concerne le **complément contraint**. Les réponses
archivées et leur qualification portent sur **0–40 Hz**. La direction
retenue et les ports restent dans le système réduit global ; le
certificat spectral seul ne contrôle pas les résonances de ce système.

Les contrôles de réponse sont évalués en binary64. Leur confrontation
utilise la marge déclarée `64*epsilon*norme_reference`. Ils ne constituent
pas une certification machine des réponses. Les certificats d'inertie
et de trace portent sur les matrices discrètes exactes, pas sur le
modèle physique continu.

Les témoins HCB emploient le calcul modal officiel Exudyn et le pont de
réponse local corrigé archivé. Le lot ne comprend aucune simulation FFRF
ou comparaison de systèmes multicorps complets. Les mesures mémoire ne
sont pas comparables entre familles. Les résultats sur ces consoles
n'établissent pas une supériorité généraliste sur Exudyn ou MBDyn.

## Contenu et intégrité

- `rapport.json.gz` conserve les octets originaux du rapport complet :
  les 84 résultats, les retours à chaque fréquence, les erreurs et les
  majorants par charge, les refus, les environnements et les terminaux.
- `provenance-fichiers.json` relie chaque fichier brut à sa représentation
  compacte et à son SHA/sa taille d'origine. Les 84 `resultat.json` sont
  extraits du rapport avec `json.dumps(..., ensure_ascii=False, indent=2,
  allow_nan=False)+'\n'` ; leurs octets reconstruits ont été comparés aux
  originaux avant suppression des doublons.
- `entrees/` contient les trois modèles NPZ, leurs métadonnées et la
  qualification des références Decimal70/90 réutilisées du lot historique
  `confrontation-ports-exudyn-0.10.0`.
- `essais/` conserve les bases NPZ (B, Phi, base physique, K/M réduits et
  indices), les journaux, les terminaux et les résultats historiques
  bruts des témoins. Les journaux et JSON volumineux sont comprimés.
- `sources/` conserve les 14 fichiers exactement mesurés. Les empreintes
  du paquet Vinkulum et des dépendances figurent dans les environnements.
- `manifeste-campagne-brute.json` est le manifeste original de campagne,
  sans modification. `exclus.json` décrit les **87 NPY exclus** : 84 champs
  et trois oracles, avec chemins d'origine, SHA, tailles, formes et indices
  de fréquences valides pour les champs. Les originaux restent dans `/tmp`.
- `manifest.json` couvre tous les fichiers compacts, sauf lui-même, par
  SHA-256 et taille. Le journal du parent est aussi conservé en gzip.

Les roues, interpréteurs et bibliothèques partagées ne sont pas copiés dans
l'archive ; leurs versions et empreintes sont conservées. Les caches des
sondes antérieures sont dans une
[archive exploratoire séparée](../inertie-contrainte-sondes-2026/README.md).

## Vérification autonome

Depuis la racine du dépôt, avec les dépendances Python installées :

```sh
python ci/archive_inertie_contrainte.py --verifier docs/bancs/inertie-contraint-2026
python ci/test_archive_inertie_contrainte.py
```

Le validateur n'accède à aucun champ ou oracle NPY exclu et n'exécute
aucune source capturée dans l'archive. Il contrôle la couverture et les
empreintes, requalifie les erreurs et majorants enregistrés, recoupe les
normes absolues/relatives, les marges et les chronos, puis réanalyse le bilan.

Pour l'inertie, il vérifie indépendamment les signatures et les
déterminants dirigés des pivots. Il **recalcule aussi les congruences
depuis les D/M/B archivés** avec le module local connu dont le SHA est
comparé au module mesuré. Les 24 certificats d'inertie possèdent trois
identités distinctes, une par maillage : trois recalculs suffisent. Tous
les champs déterministes des certificats doivent concorder, y compris
les pivots, l'ordre, les empreintes, la preuve de masse et les compteurs.
Seuls les temps de préparation et de phases sont exclus de cette égalité.

Ce rejeu reconstitue le même algorithme dirigé. Les contre-épreuves
Fraction des carnets de preuve et d'implémentation apportent un contrôle
algorithmique distinct sur de petites matrices. Le validateur de trace
contrôle l'arithmétique finale des intervalles archivés sans reconstruire
U/H/J. Aucun de ces contrôles d'archive ne reconstitue les erreurs
physiques à partir des grands champs absents.

## Reconstitution d'une campagne

Les commandes originales, chemins des interpréteurs, paramètres et SHA
figurent dans chaque terminal du rapport. Pour lancer un nouveau lot,
il faut un dossier de sortie absent, les environnements dont les
identités correspondent à celles enregistrées, les trois modèles et les
trois références complètes dont les SHA sont dans `exclus.json`. Ces
références peuvent être réutilisées depuis les originaux ou régénérées
avec les sources `confronte_ports_exudyn.py`, `oracle_champs_ports.py` et
`reference_ports_precision.py` capturées ; leur SHA doit être vérifié.

Le pilote `ci/experience_inertie_contrainte.py --campagne` expose les
arguments d'entrée et d'interpréteurs. La reproduction doit employer ses
sources figées et indiquer explicitement le bilan HCB historique si elles
sont exécutées depuis un autre emplacement. De nouveaux chronométrages
ne seront pas identiques à ceux archivés.

L'archive présente a été créée après confirmation de la terminaison des
84 processus, depuis `/tmp/vinkulum-campagne-inertie-contrainte-2026`.
Les sources mesurées et le dossier brut n'ont pas été modifiés pendant
l'archivage.

Après la campagne, le pilote courant a reçu un correctif de sérialisation
des refus : il conserve désormais le champ `bilan` des exceptions de budget
lorsque `diagnostic` est absent. Le pilote mesuré archivé, d'empreinte
`26c08cdd54ff3ee338147a7e4620ba155f11f096b92c467727c8115cc0d89ee1`,
est conservé inchangé. Le correctif ultérieur n'a pas été appliqué aux
sources capturées ni aux résultats des 84 essais.
