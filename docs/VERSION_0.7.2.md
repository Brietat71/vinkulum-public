# Vinkulum 0.7.2 — initialisation et analyse des contraintes

Cette corrective rétablit le calcul des accélérations initiales à petite
échelle, permet l'analyse de mécanismes redondants auparavant refusés et
libère les contacts en décollement lors du redémarrage après impact.
Elle conserve les signatures Rust et Python, les unités et les formats.

Le numéro **0.7.2** suit les [règles du projet](VERSIONNEMENT.md) : les
résultats corrigés changent, sans extension de l'API ni modification
volontaire du modèle physique. Python **3.14 ou plus** reste requis.

## Changements observables

- Une masse liée à une bielle de **1 µm**, avec une vitesse tangentielle de
  **10 m/s**, reçoit désormais l'accélération centripète attendue :
  **−1e8 m/s²**, contre **−7,0710678e7 m/s²** en 0.7.1.
- `k_c_m_z()` peut calculer les matrices de la cascade de parallélogrammes
  sans multiplicateurs préalablement initialisés. Le LU singulier est
  remplacé, quand nécessaire, par une projection dans la métrique de masse.
- Les commandes sont différenciées avec la géométrie. Le calcul de
  courbure n'utilise plus un déplacement d'essai fixé à `10⁻⁷ s`.
- Après un rebond, l'accélération lisse de la bille retrouve **−g** dès le
  décollement. La sélection des contacts traite aussi les réactions qui
  deviendraient attractives et les appuis couplés.

Le [rapport technique](INITIALISATION.md) détaille les équations, les
références indépendantes et les limites. Le LU ordinaire est conservé
quand il satisfait les contrôles ; ses facteurs servent aussi au
raffinement de résidus compensés. Le repli SVD travaille sur les
composantes indépendantes du gradient pondéré par la masse, sans matrice
normale globale. Il préserve les corps gelés et les lignes inactives.

## Compatibilité et limites

Aucune migration des appels n'est nécessaire depuis 0.7.1. L'API compte
toujours **185 entrées** et les dépendances sont conservées. La nouvelle
extension doit être compilée ou installée pour utiliser les corrections.

Une accélération de contrainte incompatible avec le rang numérique, une
masse non définie positive sur le repli, ou un échec du contrôle de
complémentarité produit une erreur explicite. Les mécanismes valides du
corpus existant doivent continuer à passer ; le premier refus injustifié
du pendule vertical URDF a été corrigé par raffinement.

Les matrices d'analyse et les décompositions des composantes connectées
restent denses. La norme minimale des réactions est une propriété du repli
orthogonal, sans certification indépendante de tous les LU acceptés.
Les cassures de commandes et de géométrie n'ont pas de dérivée seconde
classique. Le frottement du redémarrage dépend encore de la dernière
réaction connue. Ces limites ne sont pas supprimées par cette version.

## Validation et livraison

La CI étendue passe : **66 tests Rust**, **72 tests Python**, **41/41 groupes
de vérification**, **46/46 bancs mécaniques**, **9/9 contacts**, formatage,
Clippy et contrôle des **185 entrées d'API**.

Les **138 essais alternés** couvrent 23 configurations. Les 69 essais de
la nouvelle version passent ; la roue 0.7.1 échoue dans 30 essais par
singularité. La cascade de 144 corps est analysée en **0,83–0,85 s**, avec
**214–215 Mio** de pic RSS. Les cas déjà valides restent proches en temps :
variations des médianes de **−4,1 % à +1,5 %**. Les temps des échecs ne sont
pas utilisés pour annoncer un gain.

Le [relevé de livraison](bancs/version-0.7.2.json) lie les sources, mesures,
journaux et roue indépendante. Son installation hors du dépôt vérifie
l'exemple du README et rejoue les 41 groupes de vérification et 72 tests
Python. Le tag annoté `v0.7.2` identifie cette livraison. Le paquet reste
privé, sans publication sur un index.

La dernière confrontation MBDyn archivée reste celle de **0.6.2**. Aucune
exécution Simpack n'est ajoutée ici. L'objectif généraliste de dépasser ces
références reste actif.
