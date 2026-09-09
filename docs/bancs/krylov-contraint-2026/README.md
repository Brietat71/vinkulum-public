# Campagne Krylov contraint, 8 septembre 2026

Cette archive conserve les **60 essais** de la campagne : trois maillages,
cinq variantes, une chauffe et trois mesures. Les 24 essais Krylov passent
le juge physique ; les 12 essais avec contrôle passent aussi leurs
majorants. Parmi les 36 témoins, dix essais HCB échouent en précision.
Tous les résultats sont conservés. Le [guide et les tableaux](../../KRYLOV_CONTRAINT_PROTOTYPE.md)
présentent les coûts et les limites de comparaison.

La campagne brute a terminé sans relance. Ses 12 sources ont été capturées
avant le premier essai, puis leurs empreintes contrôlées dans chaque
processus. La préparation est fraîche : QR, sélection de la direction,
certificat à 16 chiffres et Krylov sont chronométrés. Les versions,
roues, fichiers exécutables et bibliothèques numériques sont identifiés
dans chaque résultat. Les sources capturées n'inspectent pas
l'implémentation d'Exudyn ; elles utilisent son API officielle HCB.

## Contenu sans perte des diagnostics

- `rapport.json.gz` : rapport original intégral, avec les 60 résultats,
  erreurs par fréquence et charge, contrôles, temps et environnements.
- `bilan.json` : médianes, plages et qualifications par configuration.
- `entrees/` : les trois modèles NPZ et leurs métadonnées historiques.
- `sources/` : les 12 sources effectivement mesurées.
- `essais/` : bases physiques, B/Phi, petits K/M, témoins historiques
  bruts, journaux compressés et états terminaux des processus.
- `provenance-fichiers.json` : correspondance entre les fichiers bruts
  et leur représentation compacte, avec SHA-256 et tailles originales.
- `manifeste-campagne-brute.json` : empreintes de la campagne complète.
- `exclus.json` : SHA-256, taille, forme et chemin d'origine des
  **63 grands NPY exclus** : 60 champs et trois oracles.
- `manifest.json` : inventaire exact et empreintes de l'archive compacte.

Les 60 `resultat.json` ne sont pas dupliqués : ils se reconstruisent
exactement depuis les objets `essais[*].resultat` du rapport, avec
`json.dumps(..., ensure_ascii=False, indent=2, allow_nan=False) + "\n"`.
La reconstruction a été comparée aux octets originaux avant le retrait
des doublons et ses empreintes sont vérifiées à chaque audit. L'archive
occupe environ 31 Mo, contre environ 1,1 Go pour la campagne brute.
Les journaux et JSON compressés restituent également leurs octets originaux.

## Auditer

Depuis la racine du dépôt, avec NumPy/SciPy disponibles :

```bash
python ci/archive_krylov_contraint.py --verifier docs/bancs/krylov-contraint-2026
python ci/test_archive_krylov_contraint.py
```

Le vérificateur contrôle l'inventaire, la couverture des 60 essais,
les identités, les temps, les erreurs physiques, les décisions et
les bases. Il recalcule les empreintes de B/Phi, la séparation rationnelle
de bande, les inégalités finales des certificats, les identités
d'enveloppe, les marges, les majorants relatifs et leur confrontation
aux erreurs absolues enregistrées. Les huit tests incluent des
falsifications dont les empreintes ont été recalculées : la seule
vérification SHA ne suffirait pas à les détecter.

Il **ne rejoue pas** les solveurs, les intervalles intermédiaires U/H/J,
les grands champs ni les oracles. Il ne lit aucun chemin d'origine des
NPY exclus et n'exécute pas les sources archivées. L'audit vérifie les
diagnostics enregistrés ; il n'ajoute aucune certification machine aux
majorants de réponse flottants. L'oracle Decimal 70/90 reste une
contre-vérification de précision, sans encadrement par intervalles.

## Reproduire les calculs

Le [guide](../../KRYLOV_CONTRAINT_PROTOTYPE.md#reproduire-et-auditer)
donne la commande du pilote figé et les versions requises. Préparer un
dossier neuf contenant les modèles de `entrees/` et leurs métadonnées.
Les oracles peuvent être régénérés hors chronomètre, pour chaque N parmi
32, 128 et 512, par la source capturée :

```bash
"$PY_VINKULUM" docs/bancs/krylov-contraint-2026/sources/confronte_ports_exudyn.py \
  --reference "$ENTREES_AVEC_ORACLES/n${N}-f40.npz" \
  "$ENTREES_AVEC_ORACLES/n${N}-f40.ref.npy"
```

La commande écrit aussi les métadonnées de référence. Le pilote impose
les SHA historiques des modèles et oracles avant de mesurer ; une
reconstruction différente ne passe pas silencieusement. Il refuse un
dossier de sortie déjà présent. Les installations des moteurs ne sont
pas embarquées dans cette archive.

Les [sondes exploratoires](../krylov-contraint-sondes-2026/README.md)
restent séparées : leurs anciennes sources ne sont pas toutes figées,
et leurs temps partiels ne constituent pas des mesures comparatives.
