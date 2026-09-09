# Vinkulum 0.12.2 — garde exact des contraintes de vitesse

La 0.12.2 corrige une acceptation à tort, due à l'arrondi du budget de
tolérance, dans l'assemblage des vitesses. Elle ajoute un garde entier exact
au chemin de succès et un [dossier de certification du noyau](CERTIFICATION_NOYAU.md).

**Le noyau entier n'est pas certifié.** Cette version fournit une garantie
locale sur un résidu, sous une base de confiance explicitée. Elle ne fournit
ni une preuve formelle du programme complet ni une certification industrielle.

## Défaut reproduit sur la roue 0.12.1

Une translation impose `u=0`. On initialise :

    u = 1 − 2⁻⁵³
    tol = 1 − 2⁻⁴⁶ − 2⁻⁵³

Le budget exact `tol + 64 ε |u|` vaut `u − 2⁻⁹⁹`. Le résidu est donc
strictement hors budget. L'addition flottante arrondit pourtant ce budget
vers `u` : la roue 0.12.1 laisse la vitesse inchangée.

La 0.12.2 refuse ce candidat et corrige la vitesse à zéro. Le test de
régression appelle l'API publique et échoue effectivement sur la roue
0.12.1 figée. Il ne s'agit pas d'une nouvelle erreur macroscopique de
dynamique : l'écart au budget de ce contre-exemple est `2⁻⁹⁹`.

## Garantie et intégration

Pour chaque ligne active, le code vérifie exactement :

    |d_i + Σ_j G_ij u_j| ≤ tol + 2⁻⁴⁶ (|d_i| + Σ_j |G_ij u_j|)

Les binary64 sont décodés en entiers dyadiques ; produits et sommes sont
accumulés dans une unité commune `2⁻²¹⁴⁸` avec `num-bigint` 0.4.8,
déjà présent transitivement dans le verrou Cargo. La décision finale
n'effectue aucun calcul flottant. La démonstration détaille conversion,
absence de débordement du produit des mantisses et comparaison entière.

Le résidu compensé propose les candidats ; chaque succès passe ensuite
par ce garde avant l'écriture des vitesses. Les lignes redondantes et non
holonomes restent incluses. La restitution de l'état lors d'un refus
conserve le contrat de la 0.12.1.

Le budget de tolérance, les signatures publiques et les lois mécaniques
ne changent pas. Une entrée auparavant acceptée à la faveur d'un arrondi
peut désormais demander une correction supplémentaire ou être refusée.
C'est une correction du contrat existant, d'où le numéro **0.12.2**.
L'entrée Python privée de contre-épreuve n'est pas une nouvelle API publique.

## Preuves et limites

L'oracle indépendant utilise les rationnels de CPython, pas le décodage
Rust ni ses entiers. Il vérifie les entiers témoins et la décision sans
tolérance. La campagne couvre les exposants finis, des vecteurs
pseudo-aléatoires reproductibles, des frontières, des sous-normaux et
des annulations de produits qui déborderaient en flottant.

`ci/certifie_noyau.py` est obligatoire dans la CI locale et présent dans
le modèle GitHub. Son option `--exiger-couverture-complete` rend le code
**2** : les neuf obligations ouvertes interdisent le label global.

La preuve ne couvre pas les erreurs de calcul de `G` ou de `∂Φ/∂t`,
l'erreur en avant de la projection, la position, les accélérations,
l'intégration ou les contacts. La démonstration n'est pas mécanisée dans
un assistant de preuve ; le code et ses dépendances restent dans la base
de confiance. Les certificats spectraux existants gardent leur portée propre.

Le [manifeste de livraison](bancs/version-0.12.2.json) identifie sources,
extension installée, roue, environnement, résultats et limites. L'[archive
des preuves](bancs/version-0.12.2-preuves.json.gz) contient les journaux,
les témoins exacts, la campagne physique et le générateur du manifeste.

## Validation de la livraison

- 13 314 contre-épreuves exactes : 6 666 acceptations et 6 648 refus,
  aucune divergence avec l'oracle ; 7 entrées invalides refusées.
- 87 tests Rust, 8 tests du prototype de poutre et 145 tests Python sur
  la roue installée dans un environnement neuf, hors du dépôt.
- 41 cas de vérification, 46 bancs étendus et 9 cas de contact réussis.
- 2 304 projections et six trajectoires de pendule : résultats identiques
  à ceux archivés pour la 0.12.1, avec les mêmes critères. Les deux
  explorations de conditionnement extrême restent des refus atomiques.
- Trois tests d'intégration Exudyn ignorés, le binaire étant absent.
  Les contrôles algébriques correspondants sont exécutés.

La CI étendue termine avec le code 0. Le refus du label global termine
avec le code 2 attendu. La roue validée est
`vinkulum-0.12.2-cp314-cp314-linux_x86_64.whl`, SHA-256 :
`fb78e96446f04a69e4e22f95228ba0ec18e51b70d638b36ec7b0f5e2a41daaa6`.
Environnement : CPython 3.14.7, NumPy 2.5.3, SciPy 1.18.1, Rust 1.98.0,
Linux x86_64. Les autres plateformes ne sont pas qualifiées par cette livraison.

## Coût mesuré, sans extrapolation aux trajectoires

`ci/mesure_certification_vitesse.py` mesure des appels répétés à `assemble`
sur un état déjà convergé, construction et première projection exclues.
Deux familles distinguent une vitesse normale nulle et un levier mobile
avec deux produits non nuls qui se compensent. La seconde empêche de
mesurer seulement le raccourci des facteurs nuls.

Mesures après la CI, une CPU imposée, ordre des roues 0.12.1 / 0.12.2 /
0.12.2 / 0.12.1, médiane de 18 lots par taille et version :

| Famille | Corps | 0.12.1, µs/appel | 0.12.2, µs/appel | Rapport |
|---|---:|---:|---:|---:|
| Translation, vitesse normale nulle | 1 | 0,870 | 1,153 | 1,326 |
| Translation, vitesse normale nulle | 10 | 10,757 | 14,585 | 1,356 |
| Translation, vitesse normale nulle | 50 | 200,280 | 253,741 | 1,267 |
| Translation, vitesse normale nulle | 100 | 785,519 | 979,709 | 1,247 |
| Levier mobile | 1 | 0,868 | 1,482 | 1,708 |
| Levier mobile | 10 | 10,756 | 17,500 | 1,627 |
| Levier mobile | 50 | 200,208 | 268,372 | 1,340 |
| Levier mobile | 100 | 785,293 | 1007,993 | 1,284 |

Le surcoût existe : **environ 25 à 71 % sur ces appels déjà convergés**.
Ces pourcentages ne caractérisent ni un assemblage nécessitant Newton,
ni une simulation complète, ni tous les motifs de matrices possibles.
Le garde est ajouté à l'assemblage des vitesses ; il n'est pas ajouté à
chaque pas de temps. Aucun gain de performance n'est revendiqué pour ce lot.
