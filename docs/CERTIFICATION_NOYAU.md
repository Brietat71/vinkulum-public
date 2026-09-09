# Certification du noyau : contrat, preuves et obligations

État de la démarche, 9 septembre 2026. **Le noyau entier n'est pas certifié.**
Ce document engage une certification mathématique et numérique par propriétés
précises. Il ne constitue ni une attestation industrielle, ni une conformité
déclarée à un référentiel externe. La vocation généraliste du noyau est
conservée ; chaque domaine mécanique doit avoir son périmètre de preuve.

La **0.14.0** ajoute la [preuve des quotients exactement redondants](CERTIFICATION_QUOTIENT.md) :
rang et compatibilité exacts, accélération et réaction généralisée uniques,
bornes sur le représentant natif des multiplicateurs. Cette propriété ne
déclare pas uniques des multiplicateurs qui ne le sont pas.

La **0.13.0** ajoute des [certificats d'unicité et d'erreur en avant](CERTIFICATION_LINEAIRE.md)
pour les systèmes linéaires carrés admis, avec export du système d'accélération
initiale réellement résolu et vérificateur autonome. Ce périmètre est
facultatif et limité ; les obligations générales ci-dessous restent ouvertes.

## 1. Ce qui change dans le code

Depuis 0.12.2, un succès de `Noyau.assemble(vitesses=True)` passe par un
**garde en arithmétique entière exacte** avant l'écriture des vitesses.
Chaque équation originale active est contrôlée, y compris les redondances
retirées de la factorisation dynamique et les contraintes non holonomes.
Les éléments explicitement gelés sont masqués conformément au contrat.

Pour les coefficients binary64 effectivement calculés, le garde décide :

    r_i = d_i + Σ_j G_ij u_j
    s_i = |d_i| + Σ_j |G_ij u_j|
    |r_i| ≤ τ + 2⁻⁴⁶ s_i

Ici `d = ∂Φ/∂t`, `τ = tol` et `2⁻⁴⁶ = 64 ε` en binary64. Les produits,
sommes et comparaison **de ce contrôle** sont exacts. La marge relative
existante n'est pas augmentée. Le calcul compensé précédent reste un
préfiltre : il peut proposer un candidat au garde, jamais autoriser seul
l'écriture. Un échec déclenche le raffinement existant ou le refus atomique.

**Portée :** certificat a posteriori de l'inégalité sur les données stockées.
Ce n'est pas un certificat des fonctions géométriques ayant produit `G` et
`d`, ni de la distance à la variété admissible, ni de la norme minimale de
la projection, ni des accélérations ou trajectoires. En particulier, un
grand `s_i` autorise un résidu absolu plus grand que `tol`. Un petit résidu
ne borne pas l'erreur en avant sans information de conditionnement.
`vitesses=False` et l'API Rust `Modele::assemble`, limitée aux positions,
n'acquièrent pas cette garantie. Les réassemblages lors d'un changement du
masque de contraintes utilisent aussi ce garde. Les autres projections des
pas de temps et leurs erreurs de trajectoire ne sont pas couvertes par cette
seule propriété.

## 2. Démonstration du contrôle entier

Le vérificateur est dans `src/certificat_vitesse.rs`. Il ne fait aucune
opération flottante après la lecture des bits.

Un binary64 fini s'écrit exactement `(-1)^s m 2^e`. Pour un champ d'exposant
`E` et une fraction `F` :

| Cas | Mantisse entière m | Exposant e |
|---|---:|---:|
| `E=0` (zéro ou sous-normal) | `F` | `−1074` |
| `1 ≤ E ≤ 2046` | `2⁵² + F` | `E−1075` |
| `E=2047` (infini ou NaN) | refus | — |

Ainsi `0 ≤ m < 2⁵³` et `−1074 ≤ e ≤ 971`. Le produit de deux mantisses
tient dans un `u128` car il est strictement inférieur à `2¹⁰⁶`. Tout produit
de deux entrées est un entier en unités de `q=2⁻²¹⁴⁸` ; son décalage
`e₁+e₂+2148` est compris entre 0 et 4090. Les termes simples `d` et `τ`
sont représentables dans la même unité. Aucun décalage négatif n'est possible.

Le programme accumule `R=r/q`, `S=s/q` et `T=τ/q` avec `num-bigint`.
Par induction sur les termes, l'addition de leur entier signé conserve
exactement `R`, et celle de leurs magnitudes conserve exactement `S`.
Les zéros signés ont une magnitude nulle et ne modifient aucune somme.
Puisque `q>0`, l'inégalité cible est équivalente à :

    2⁴⁶ |R| ≤ 2⁴⁶ T + S

C'est exactement la comparaison effectuée par `CertificatLigne::accepte`.
Elle reste définie lorsque des produits flottants déborderaient ou
sous-déborderaient. Le solveur peut cependant refuser ces situations en
amont : la preuve d'acceptation n'est pas une promesse de terminaison réussie.

La mémoire des entiers dépend logarithmiquement du nombre de termes ; le
calcul parcourt toutes les colonnes de chaque ligne contrôlée. Une panne
mémoire ou matérielle n'entre pas dans cette preuve d'acceptation. Aucune
garantie de temps réel n'est revendiquée.

**Niveau de preuve :** démonstration mathématique inspectable et calcul
exact par exécution, sous les hypothèses d'implémentation ci-dessous.
Le [dossier Lean 4](../preuves/README.md), ajouté après la 0.14.1,
formalise le décodage binary64 et l'équivalence du garde entier au critère
rationnel : 13 théorèmes audités et 13 314 contre-épreuves Lean/Rust/Fraction.
Le code Rust compilé et `num-bigint` ne sont pas formellement vérifiés ici.
Les contre-épreuves n'établissent pas à elles seules une preuve universelle.

## 3. Base de confiance et vérification indépendante

La garantie suppose la fidélité de `f64::to_bits` au format binary64,
les opérations entières Rust, `num-bigint` 0.4.8, le compilateur, le matériel,
et la connexion du garde à toutes les sorties de succès de cette projection.
Elle ne suppose pas que QR/SVD calcule une bonne correction : une correction
incorrecte doit être rejetée si elle viole l'inégalité.

`ci/certifie_noyau.py` appelle le **même vérificateur compilé** via une entrée
privée. Son oracle utilise `fractions.Fraction` et la conversion rationnelle
de CPython, indépendamment du décodage de bits et des entiers Rust. Il compare
les trois entiers témoins et la décision, sans tolérance :

- tous les 2047 champs d'exposant finis, deux signes, deux mantisses extrêmes ;
- produits avec exposants variés et 2048 vecteurs pseudo-aléatoires reproductibles ;
- 3072 cas autour des frontières d'acceptation ;
- sous-normaux, débordements annulés, annulation de 1024 grands produits ;
- refus des dimensions incohérentes, NaN, infinis et tolérances négatives.

Un test du **chemin public** construit une translation dont l'ancien budget
flottant arrondit vers le résidu. La 0.12.1 laisse sa vitesse inchangée ;
la 0.12.2 doit effectuer la correction. Ce test détecte notamment le
contournement du garde exact dans le chemin de succès. Les autres tests
d'assemblage vérifient redondances, gel, non-holonomie et restitution d'état.

Le contre-exemple prend `u=1−2⁻⁵³`, `G=1`, `d=0` et
`tol=1−2⁻⁴⁶−2⁻⁵³`. Le budget exact vaut `u−2⁻⁹⁹` : le résidu est donc
hors budget de `2⁻⁹⁹`, alors que l'addition flottante arrondit le budget à `u`.
La taille de cet écart ne doit pas être confondue avec celle du défaut
d'initialisation dynamique corrigé en 0.12.1.

Le dossier JSON contient l'empreinte de l'extension effectivement chargée,
les empreintes des sources présentes, la graine, l'empreinte de la campagne
et des témoins exacts. Les empreintes identifient les objets examinés :
elles ne prouvent pas, à elles seules, que l'extension a été construite
depuis ces sources. Le dossier de livraison doit ajouter construction,
roue installée et journaux de validation. Une autre plateforme requiert
sa propre exécution ; les résultats Linux ne lui sont pas attribués.

```bash
python ci/certifie_noyau.py --sortie certification.json
python ci/certifie_noyau.py --sortie certification.json --exiger-couverture-complete
```

La première commande échoue à la moindre divergence avec l'oracle et
énonce les obligations ouvertes. La seconde rend **2** tant que la couverture
complète n'est pas établie. Un code 0 de la campagne ordinaire n'est donc
jamais présenté comme une certification globale. La campagne fait partie
de la CI locale obligatoire avant push et du modèle de CI GitHub.

## 4. Registre des garanties et des lacunes

| Propriété | Preuve ou mesure disponible | Ce qui manque pour la certifier |
|---|---|---|
| Résidu de vitesse à la sortie de l'assemblage Python | Calcul entier exact ci-dessus ; contre-épreuves rationnelles | Preuve formelle de l'implémentation et revue indépendante |
| Géométrie et lois temporelles, SO(3), G, ∂Φ/∂t | 0.15.0 : géométrie holonome admise et Jacobien encadrés, norme des quaternions incluse ; tests de tangentes et covariance | Autres fonctions mécaniques, lois imposées et incertitudes de construction des coefficients |
| Assemblage des positions | 0.15.0 : existence et unicité locales par Krawczyk, bornes de distance, jauges explicites | Autres mécanismes, grandes dimensions et garantie de convergence de l'algorithme natif |
| Rang, contraintes compatibles, systèmes KKT | Quotient exact dense ; 0.16.0 : facteurs creux et enveloppes, quotient structurel perturbé avec rang préservé | Perturbations sans structure de rang, systèmes hors budgets et justification mécanique des enveloppes |
| Initialisation dynamique | Certificat facultatif des accélérations, du représentant des multiplicateurs et de la réaction généralisée | Contacts, partitions et erreurs des coefficients mécaniques |
| Intégration α-généralisée/GGL | Premier pendule α-généralisé sur la roue 0.15.0 : bornes rationnelles aux nœuds sur 1 s, trois grilles, incertitudes initiales et arrondis de référence inclus | Autres mécanismes/horizons, GGL, interpolations entre nœuds, contraintes et bilans |
| Contacts et frottement | Bancs de contact | Garanties aux transitions, complémentarité et dissipation discrète |
| Réduction flexible et masse | Certificats d'inertie/comparaison sur matrices stockées | Propagation des arrondis jusqu'aux champs et réponses physiques |
| Refus et restitution, interfaces Rust/C++/Python | Sauvegardes et tests de régression, validations de buffers | Preuve de toutes les sorties d'échec et sûreté des interfaces |
| Validité physique des modèles généralistes | Campagnes référencées, couverture inégale selon les modèles | Domaines d'usage, incertitudes et validation expérimentale indépendante |

Les `certification_machine=True` existants des modules spectraux conservent
leur portée locale. Les réponses physiques restent explicitement non
certifiées contre les arrondis. Les rapports historiques ne sont pas des
validations automatiques de la version courante. Voir notamment
`INVERSE_SELECTIONNEE_PREUVES.md`, `MASSE_COMPAREE_PREUVES.md`,
`ASSEMBLAGE_COMPLEMENTS_PREUVES.md` et `VERSION_0.12.1.md`.

En particulier, le rapport historique `VALIDATION.md`, introduit au commit
`afb277b` du 6 septembre 2026, annonce **44/62 lignes dans leur tolérance**.
Ses 18 écarts, notamment aérodynamiques, restent des écarts publiés. Ce
rapport n'identifie pas une roue 0.12.2 et n'est pas réattribué à cette
livraison. Ce rapport historique seul ne qualifie pas une nouvelle roue ;
le garde cinématique ne résout pas les écarts physiques.

La [qualification 0.14.1](VERSION_0.14.1.md) rejoue désormais cette campagne
sur les roues identifiées 0.14.0 et 0.14.1 : 62 lignes reproduites,
44 dans leur tolérance et 18 écarts, avec les mêmes valeurs et critères.
Cette reprise ferme la lacune de traçabilité pour ces expériences ; elle
ne résout pas les désaccords physiques. La correction du Newton des petites
inerties et les contre-épreuves lisses n'ajoutent pas de preuve temporelle.

Les dossiers [d'assemblage géométrique](CERTIFICATION_ASSEMBLAGE.md) et de
[garantie temporelle du pendule](GARANTIE_TEMPORELLE_PENDULE.md) décrivent
les garanties ajoutées depuis, avec leurs hypothèses et leur base de
confiance. La preuve temporelle est un argument mathématique accompagné
d'un vérificateur rationnel et de contre-épreuves analytiques ; elle n'est
pas formalisée dans Lean et ne ferme pas la ligne d'intégration générale.

## 5. Ordre des obligations à fermer

1. Faire relire indépendamment le garde exact et sa formalisation Lean
   désormais disponible ; renforcer le lien entre cette preuve et le code
   compilé, au-delà des contre-épreuves de correspondance.
2. Étendre les certificats linéaires et les quotients exacts de la 0.14.0
   aux opérateurs creux et aux dépendances avec perturbations des coefficients.
   Les quotients denses admis sont désormais équipés d'une borne en avant.
3. Encadrer géométrie et dérivées ; certifier localement Newton lorsque les
   hypothèses d'existence/unicité sont vérifiées. Refuser la certification
   si ces hypothèses ne sont pas établies, sans prétendre que le modèle est faux.
4. Construire un contrôle a posteriori des trajectoires et des événements,
   incluant erreurs de discrétisation et d'arrondi. Les lois dissipatives et
   les contacts nécessitent leurs propres théorèmes et domaines.
5. Qualifier les usages retenus avec données expérimentales, incertitudes,
   traçabilité des défauts et évaluation indépendante du dossier de livraison.

Chaque fermeture devra ajouter une affirmation précise, ses hypothèses,
un artefact contrôlable, un refus en cas d'hypothèse non satisfaite et une
contre-épreuve capable de détecter une violation. Aucun nombre de tests
verts ne ferme automatiquement une obligation mathématique.

Le [plan de fiabilité](PLAN_FIABILITE.md) distingue les premiers domaines
livrés des obligations générales de ce registre. Le
[diagnostic énergie–moment](ATTRIBUTION_REPERES_EM.md) ajoute un travail
d'attribution des différences entre trajectoires discrètes : neuf comparaisons
admises sur 420 000 pas. Cette attribution ne borne pas l'erreur de
discrétisation par rapport à la solution continue.

## 6. Référentiels et portée industrielle

Le [NASA-STD-7009B, approuvé le 5 mars 2024](https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf)
organise la crédibilité des modèles et simulations sur leur cycle de vie,
avec des critères d'acceptation rattachés au projet. Il distingue notamment
vérification, validation, incertitudes et appréciation des résultats.
Il sert ici à structurer les questions du dossier, sans déclaration de conformité.

La présentation officielle de l'[ASME V&V 10-2019 (R2025)](https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-solid-mechanics)
fournit également un cadre de vérification, validation et quantification des
incertitudes en mécanique des solides. Le texte intégral payant n'a pas
été audité dans ce travail ; aucune conformité article par article n'est établie.

Une certification industrielle exige de convenir de l'usage, du périmètre,
du référentiel applicable et de l'autorité d'évaluation. Ces décisions
ne se déduisent pas de la seule vocation généraliste de Vinkulum.
