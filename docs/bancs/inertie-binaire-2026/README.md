# Inertie binary64 compilée : 36 essais frais

Campagne du 8 septembre 2026 : trois consoles, certificat compilé avec
séparateurs BFS, certificat Decimal témoin et LU corrigée ; une chauffe et
trois mesures par configuration, CPU 8 et bibliothèques à un fil.
Les 36 essais passent le juge physique ; les 24 contrôles passent leurs
majorants. La permutation est reconstruite et comptée dans chaque candidat
binaire. Aucun nouveau classement Exudyn, MBDyn ou Simpack dans ce lot.

`rapport.json.gz` contient tous les résultats et diagnostics. Les JSON,
journaux et la bibliothèque compilée sont compressés sans changer leurs
octets décompressés. Sources, bases B/Phi et modèles NPZ sont conservés.
Les 36 grands champs NPY sont exclus, avec SHA et taille dans le manifeste ;
ils restent dans `/tmp/vinkulum-campagne-inertie-binaire-2026`. Les trois
oracles externes restent identifiés dans le rapport et leurs métadonnées.
L'archive ne suffit donc pas à rejuger les champs sans reconstruction.

Le validateur local requalifie les diagnostics des 36 essais, vérifie les
empreintes et l'identité des champs/normes/majorants entre contrôleurs,
puis reconstruit les six identités de certificats depuis les modèles et B.
Il compile le C++ courant connu du dépôt ; il n'exécute ni source ni
bibliothèque provenant de l'archive. Les durées et les chemins/empreintes
de la bibliothèque reconstruite ne sont pas des éléments de la preuve.
Les pivots binary64 sont encodés en hexadécimal exact.

La compilation unique et le chargement du module sont hors chronomètre,
comme les imports des autres solveurs ; les commandes, le compilateur,
le source et le binaire mesuré sont identifiés dans `build/`.
La réponse et ses majorants restent non certifiés machine. La certification
porte uniquement sur la coercivité du complément discret exact.

Voir [protocole, preuves et résultats](../../INERTIE_BINAIRE_SEPARATEURS.md).
