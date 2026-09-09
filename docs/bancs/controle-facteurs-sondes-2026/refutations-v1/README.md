# Deux réfutations avant correction du contrôle par facteurs

8 septembre 2026. Snapshot de `ci/controle_facteurs.py` avant correction :
SHA256 `861ecf586bc36344a4e5673c3c1450ab66bbea482f16886f20307ac7e99d99ed`.
Le script de reproduction charge ce fichier local, pas le module corrigé.
Les auxiliaires historiques importés proviennent du `ci` du clone ; aucun
solveur natif, champ mécanique ou campagne comparative n'est exécuté.

## 1. La norme d'audit perd de petites contributions normales

Soit `image` de taille 128×2 : première ligne `(1,1)`, puis 127 lignes
`(2^-27,2^-27)`. Les facteurs proposés `gauche` 128×1 et `droite` 1×2
sont nuls ; le contrat autorise des facteurs imparfaits puisque leur
défaut est censé rester dans la majoration. Appliquer à `v=(1,0)^T`.

Le helper retourne `b=0x1.0000000000002p+0 = 1+2^-51`, mais

```text
||image v||² = 1 + 127·2^-54,
||image v||² − b² = 111·2^-54 − 2^-102 > 0.
```

Tous les nombres et produits de ce cas sont normaux. La somme de carrés
stridée de `np.linalg.norm(..., axis=0)` perd les petites contributions.
L'inflation utilisant le nombre de colonnes de l'image protège un autre
produit ; elle ne contrôle pas cette réduction sur 128 lignes.

L'inégalité abstraite de compression reste correcte avec des bornes
effectives des normes du défaut. Une correction doit donc protéger aussi
ces normes, par accumulation stable et inflation dépendant de leur vraie
dimension, ou par arrondis dirigés. Sous le modèle relatif classique,
hors débordement/sous-flux, une norme obtenue par somme de carrés et racine
demande de tenir compte du défaut du produit scalaire et de la racine.
Le seul remplacement par une autre routine non encadrée ne constitue pas
automatiquement un certificat machine.

## 2. Un Gram sous-flué masque une image non nulle

La routine de norme est isolée avec les dimensions physiques 3×2 :

```text
DX = [[2^-600, 0], [0, 2^-600], [0, 0]],
fl(DX.T DX) = 0,     ||DX||₂ = 2^-600 > 0.
```

Le helper retourne zéro. L'hypothèse du modèle relatif sans sous-flux
n'est pas vérifiée : les entrées sont finies et normales, mais leurs
carrés disparaissent. Une détection des diagonales nulles ou subnormales
avec repli vers une norme supérieure évaluée par mise à l'échelle permet
de traiter ce cas. Un refus explicite serait également sûr. Détecter
uniquement un Gram entièrement nul ne couvre pas les produits subnormaux
non nuls dont l'erreur relative peut aussi être grande.

La régression ajoutée au dépôt couvre également une diagonale
`1.2 * 2^-537`, dont le carré est subnormal non nul. Ces tests de helper
ne démontrent pas qu'une réponse mécanique native aurait été mal bornée.
Le module déclarait déjà `certification_machine=False` ; les réfutations
précisent pourquoi ses audits ordinaires ne peuvent être lus comme des
majorations machine inconditionnelles.

## Conservation et reproduction

`donnees.npz` contient toutes les matrices des deux cas originaux.
`resultats.json` conserve les valeurs, fractions exactes, diagnostics,
versions et empreinte du snapshot ; `execution.log` conserve la sortie.
Le programme vérifie que les deux sous-majorations sont effectivement
présentes. Il s'exécute sur CPU0 avec les bibliothèques à un fil.

```sh
PYTHONDONTWRITEBYTECODE=1 /tmp/vinkulum-release-0.11.0-final/venv/bin/python /tmp/vinkulum-controle-facteurs-refutations-v1/reproduire.py
```

Les régressions du module corrigé sont ajoutées uniquement dans
`ci/test_controle_facteurs.py`. Aucune source ancienne ou mesure archivée
n'est modifiée par cette archive. `manifest.json` donne les tailles et
empreintes de tous les autres fichiers.
