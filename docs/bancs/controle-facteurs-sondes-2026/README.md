# Sondes du contrôle par facteurs — 8 septembre 2026

Quatre groupes exploratoires, distinctes de la campagne comparative fraîche.

- `profil/` contient le profil du contrôleur historique à 512 poutres : quatre
  balayages identiques en champs et résultats, dont cProfile et attribution
  par lignes. Les directions, le certificat et la préparation commune sont
  offerts pour observer les réponses. Les coûts instrumentés ne constituent
  pas des gains comparatifs. Son README donne le détail.
- `pilote-v1/` conserve neuf essais : ancien contrôleur, extension Gram et
  extension QR sur 32, 128 et 512 poutres. Chaque maillage partage QR, sélection,
  certificat et Krylov ; les durées totales sont donc indicatives. Chaque
  contrôleur est reconstruit, et ses 257 fréquences traitent six charges.
  Les neuf champs sont identiques bit à bit au témoin de leur maillage et
  tous les contrôles passent. Le choix Gram est arrêté avant la campagne
  formelle ; l'option QR, légèrement plus coûteuse, reste conservée ici.

Les profils ont conduit à réutiliser les normes, majorer seulement la norme
*d'opérateur* de l'extension via le Gram déjà calculé, et factoriser le
correctif de contrainte en conservant le défaut des facteurs. Les normes des
petites combinaisons sont calculées sur les champs physiques retournés.
Aucune certification machine des réponses n'est revendiquée.

Les snapshots, scripts, environnements, résultats complets, journaux et
manifestes originaux sont conservés. `manifest.json` couvre les octets de
cette archive, les SHA/tailles décompressées et les chemins sources.
`exclus.json` décrit dix-huit grands champs NPY par SHA, taille, forme et chemin.
La vérification de cette archive ne rejoue pas le juge depuis ces NPY absents.
Les sources et résultats des sondes sont ceux de leurs captures ; une
correction ultérieure du prototype ne modifie pas ces octets.

Les JSON et le profil texte cProfile sont comprimés sans modification de leur contenu. Les manifestes
originaux peuvent donc être recoupés via les métadonnées de provenance,
y compris les grands fichiers exclus. Aucun essai n'est supprimé au motif
qu'il serait moins favorable. Les exécutables et roues ne sont pas copiés.

## Réfutations et correction avant la campagne

`refutations-v1/` conserve deux sous-majorations constatées sur un état
intermédiaire : norme d'une colonne de 128 lignes mal accumulée, et image
non nulle dont le Gram sous-flue à zéro (ou devient subnormal). Le snapshot,
les données et le script de reproduction sont conservés, même après correction.

`pilote-v2/` répète les neuf essais après protection des normes d'audit par
mise à l'échelle et budgets d'arrondi. Une diagonale Gram nulle/subnormale
provoque désormais un repli sur une norme de Frobenius physique robuste.
Tous les champs restent bit à bit identiques aux témoins et les majorants
restent acceptés. Cette vérification ne transforme pas le contrôleur en
certificat machine des réponses. La variante Gram reste choisie avant le
lot formel ; ses gains indicatifs ne remplacent pas les mesures fraîches.
