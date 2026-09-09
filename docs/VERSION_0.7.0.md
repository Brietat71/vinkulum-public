# Vinkulum 0.7.0 — statique contrôlable et noyau numérique renforcé

Cette version regroupe les sept commits de développement entre `v0.6.2`
et `24ba6a7`. Elle ajoute des capacités publiques : le mode statique
strict, son bilan, l'état cinématique précis et une formulation de poutre
optionnelle. Cela justifie une **mineure**, conformément aux
[règles du projet](VERSIONNEMENT.md). Des adaptations des appels Rust
directs sont également indiquées ci-dessous.

## Capacités ajoutées

- `statique(strict=True)` exige la tolérance demandée et refuse le secours
  sur stagnation. `statique_info()` distingue `tolerance`, `stagnation` et
  `echec`, avec résidus, tentatives et état restauré. Le retour historique
  `(résidu, compteur)` de `statique()` est conservé.
- `etat_precis()` conserve les petites composantes des translations.
  `pose_etat(*n.etat_precis())` permet de les restaurer. Cela corrige les
  planchers de résidu provoqués par l'arrondi des positions sur les
  chargements faibles, notamment les paliers Princeton.
- `poutre(..., formulation="integree")` propose une formulation
  expérimentale à interpolation des rotations et flexibilités condensées.
  Le défaut reste `formulation="milieu"`.

Les descriptions et domaines détaillés sont dans les rapports de
[statique stricte](STATUT_STATIQUE.md), de
[positions compensées](POSITIONS_COMPENSEES.md) et de
[poutre intégrée](POUTRE_INTEGREE.md).

## Changements du noyau

Les tangentes de force sont assemblées localement ; le chemin statique
creux évite les réévaluations et matrices globales correspondantes. Les
projecteurs et réactions emploient le QR creux lorsque leur domaine le
permet. Une réduction orthonormale traite les liaisons équivalentes ; une
factorisation orthogonale complète traite le repli des dépendances générales
avec contrôles de rang et de stationnarité.

Les corrections de SVD sont vérifiées contre des références fabriquées
et LAPACK. La normalisation par puissance de deux traite les échelles
uniformes extrêmes. La projection polaire des rotations remplace les SVD
répétées près de SO(3), ainsi que Gram–Schmidt dans le domaine des corps
libres du schéma énergie–moment.

Les références, coûts et limites sont publiés séparément :
[assemblage local](ASSEMBLAGE_LOCAL.md),
[contraintes creuses](CONTRAINTES_CREUSES.md),
[réactions redondantes](REDONDANCES_CREUSES.md),
[Newton redondant](NEWTON_REDONDANT.md),
[factorisations orthogonales](FACTORISATIONS_ORTHOGONALES.md) et
[rotations polaires](ROTATIONS_POLAIRES.md).

## Compatibilité et migration

Python reste requis en version **3.14 ou plus**. Aucune méthode Python
historique n'est supprimée ; les nouveaux arguments sont optionnels.
La référence d'API compte désormais **185 entrées**, contre 183 en 0.6.2.

En mode statique habituel (`strict=False`), une acceptation sur stagnation
émet désormais `RuntimeWarning`. Un environnement qui transforme les
avertissements en exceptions obtient une exception et une restauration.
Employer `strict=True` quand un arrêt de secours doit être refusé, et
consulter `statique_info()` pour distinguer les verdicts.

`etat()` conserve ses six champs arrondis. `etat_precis()` en ajoute un
septième, `r_bas` ; les deux parties de la position doivent rester séparées
pour préserver la précision. Cette API porte la cinématique et les états
`v_i`, sans constituer une sauvegarde complète des historiques de contact
et de toutes les variables internes.

Les utilisateurs Rust qui construisent les structures publiques directement
doivent ajouter `Corps::r_bas`, initialisé à `V3::zeros()` pour une position
ordinaire, et `Poutre::integree`, initialisé à `false` pour l'élément
historique. Les méthodes de tangente concernées par les positions compensées
prennent désormais des poses `(r, R, r_bas)` au lieu de `(r, R)` ; adapter
les appels directs aux signatures de [tangent.rs](../src/tangent.rs).
Les clients Rust doivent être recompilés.

Des configurations auparavant acceptées peuvent désormais échouer après
correction d'un faux succès numérique : déséquilibre libre masqué par de
grandes réactions, stagnation refusée en mode strict, ou matrice dont
l'étendue interne perdrait des coefficients pendant la normalisation SVD.
Une valeur singulière non représentable est également refusée. Les résultats
numériques ne sont pas promis identiques bit à bit à ceux de 0.6.2.

## Validation et identification de la livraison

Le noyau de cette livraison est celui de `24ba6a7` : la préparation de version
modifie les métadonnées et la documentation. Les mesures historiques gardent
leurs commits, empreintes et numéros réellement utilisés. Elles ne sont pas
renommées en mesures 0.7.0.

La CI étendue a réussi sur **0.7.0** : formatage, Clippy sans avertissement,
**52 tests Rust**, **66 tests Python**, vérification **41/41**, bancs
mécaniques **46/46**, contacts **9/9** et référence d'API à jour.

La roue CPython 3.14 / Linux x86_64 (`manylinux_2_39`) a été installée hors
du dépôt dans un venv neuf. Les imports isolés, les versions, le mode
statique strict, la restauration de l'état précis et le chargement d'une
poutre intégrée passent. Les données embarquées comprennent 21 modèles
MJCF et 9 URDF ; le modèle Panda se charge avec ses 11 corps. L'extension
de cette roue a la même empreinte que celle utilisée par la CI.

Les 75 empreintes de fichiers inscrites dans `RECORD` sont vérifiées. Les résultats,
versions synchronisées, sources et empreintes du paquet sont consignés dans
[le relevé 0.7.0](bancs/version-0.7.0.json).

## Limites conservées

Le repli des dépendances générales reste dense. La poutre intégrée reste
expérimentale et son câble condensé extrême conserve son échec connu.
Des surcoûts sont publiés, dont environ 2 % sur le corps libre dans la
dernière comparaison ; les écarts entre repères à long terme restent à
diagnostiquer. Aucun gain universel n'est annoncé.

La confrontation MBDyn archivée reste celle de 0.6.2. Il n'existe toujours
aucune nouvelle exécution Simpack dans ce dépôt. L'
[objectif généraliste](OBJECTIF_MBDYN.md) reste ouvert.

Le tag annoté `v0.7.0` désigne le commit de cette livraison. Le paquet reste
privé et n'est pas publié sur un index de paquets.
