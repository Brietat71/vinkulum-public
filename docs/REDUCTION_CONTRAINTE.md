# Réduction contrainte et assemblage — API 0.12.0

`vinkulum.reduction_contrainte` rend installables les développements sur le
complément spectral, les masses couplées et les assemblages. Le modèle est
**linéaire, réel, conservatif**, avec `K=D.T @ D`, masse `M`, pulsations en
rad/s et charges appliquées aux ports. La dynamique multicorps non linéaire,
l'amortissement, les charges intérieures et les impacts ne sont pas traités
par cette API fréquentielle.

## Exemple : traverser une résonance intérieure

```python
import numpy as np
from vinkulum.reduction_contrainte import ReductionContrainte

D = np.array([[1., 0., .5], [0., 4., .25], [0., 0., 2.]])
M = np.eye(3)
r = ReductionContrainte(
    D, M, interieur=[0, 1], interface=[2], metrique=[[1.]],
    omega_max=1.25, modes_retenus=1, gamma=4.,
)
rep = r.reponses(omega=1., forces=[[1.]])
print(rep['champ'])  # [[2.], [0.], [0.]], à l'arrondi près
print(rep['statut_controle'])  # borne_numerique_disponible
print(r.certificats['complement']['certification_machine'])  # True
print(rep['certification_machine'])  # False
```

À `omega=1`, le premier degré intérieur a un pôle lorsque le port est fixé.
La direction correspondante reste dans le petit système conservé. Son
couplage avec le port permet ici une réponse globale régulière. Le
certificat porte sur les directions éliminées, au seuil `gamma=4`, et
n'exige donc pas la coercivité de tout l'intérieur sur la bande.

## Contrat du constructeur

- `D` peut être rectangulaire ; `M` est carrée, exactement symétrique et de
  même dimension physique. Le stockage creux doit être canonique, sans
  doublons ni indices non triés. Aucune somme implicite ne modifie le modèle.
- `interieur` et `interface` constituent une partition complète, non vide
  dans chaque partie. La métrique des ports est exactement symétrique et
  numériquement définie positive.
- `omega_max >= 0`. `gamma > omega_max**2` est un **seuil à démontrer**, pas
  une valeur propre supposée. Par défaut : `4*omega_max**2`, ou `1` pour une
  bande limitée à zéro. Un seuil trop élevé peut être refusé.
- `modes_retenus=1` propose les directions basses par un problème propre
  généralisé utilisant la **masse complète intérieure**. Le rang est
  strictement compris entre zéro et le nombre d'intérieurs.
  `directions_retenues=Phi` fournit alternativement une matrice finie
  `(intérieurs, rang)` ; sa largeur détermine alors le rang et remplace
  la sélection automatique. `B=M_ii Phi` est calculé puis certifié tel que
  stocké. Ni `Phi` ni les valeurs d'eigsh ne sont des preuves.
- `blocs=4`, `max_directions=128`, `profondeur=8` règlent la réduction et
  l'enveloppe de contrôle. Aucun enrichissement ou assouplissement de seuil
  automatique ne masque un refus.
- `alpha_masse=.25`, `beta_masse=2.` proposent une comparaison stricte
  `alpha L < M_active < beta L`, avec `L=diag(M_active)`. Ces valeurs sont
  vérifiées, jamais présumées universelles. Une masse SPD peut être refusée
  si ces constantes ne conviennent pas.
- La masse complète doit être PSD : son bloc à diagonale positive est
  certifié SPD ; chaque diagonale nulle exige une ligne et une colonne
  **exactement nulles**. Les intérieurs sont strictement massiques. Des ports
  sans masse sont ainsi permis ; les noyaux massiques obliques restent hors
  domaine. La preuve de masse complète s'ajoute à celle de l'intérieur.
- `budget_qr=10_000_000` borne les mises à jour Givens ; les certificats
  possèdent leurs budgets `budget_operations=100_000_000`,
  `budget_coefficients=2_000_000`, `budget_rectangulaire=2_000_000`.
  Ces compteurs sont **par certificat**, et ne constituent pas un plafond
  global de mémoire ou de temps. Le graphe des séparateurs et les bases
  numériques ont aussi un coût.

Les entrées sont copiées. `certificats` et `diagnostic` rendent des copies ;
leur modification ne change pas le calcul. Les attributs préfixés `_`
sont privés. Pour changer le modèle ou les réglages, construire un nouvel
objet.

`ReductionContrainte.depuis_noyau(noyau, libres, interface, metrique,
omega_max, **options)` reprend l'extraction des facteurs matériels des
poutres. Elle exige des contraintes fixes holonomes satisfaites et
`G[:, libres] == 0` exactement. Une tangente sous précontrainte ou un état
en mouvement ne devient pas un modèle matériel équivalent par extraction.

## Résultats et refus

`reponses(omega, forces)` prend une matrice `(ports, charges)` et retourne
**un dictionnaire groupé**. Ce format est distinct de la liste par charge
de `ReductionMaterielle.reponses_groupees`, dont le contrat reste inchangé.

| Champ | Signification |
|---|---|
| `champ` | Déplacements physiques `(DDL, charges)`, dans l'ordre d'entrée |
| `coordonnees` | Coordonnées normalisées des ports et directions retenues |
| `bornes['masse'/'deformation']['normes']` | Normes des champs reconstruits |
| `bornes[...]['absolues']` | Majorants numériques d'erreur, ou `None` |
| `bornes[...]['relatives']` | Majorants relatifs numériques par charge, ou `None` |
| `marge` | Plus petite valeur singulière du bloc conservé, moins l'enveloppe d'erreur |
| `statut_controle` | `borne_numerique_disponible` ou `marge_non_positive` |
| `certification_machine` | Toujours `False` pour une réponse |

Une marge positive permet de calculer des bornes ; elle ne garantit pas
qu'elles satisfont une tolérance choisie par l'utilisateur. Une marge non
positive rend les bornes indisponibles ; un système conservé singulier
lève `numpy.linalg.LinAlgError`. Un certificat non démontré lève
`CertificatIndisponible`, avec `diagnostic`. Les données invalides lèvent
`ValueError` ; une sélection propre qui ne converge pas peut lever l'erreur
SciPy correspondante. Il n'y a pas de repli silencieux vers un autre modèle.

## Assemblage

```python
from vinkulum.reduction_contrainte import AssemblageContraint

a = AssemblageContraint(
    reductions=[r, r], applications=[[[1.]], [[1.]]], metrique=[[2.]],
    facteur_externe=[[1.]], masse_externe=[[0.]],
)
rep = a.reponses(omega=.5, forces=[[1.]])
# rep['champs'] : liste des champs locaux
# rep['champ_ports'] : déplacements physiques globaux
# rep['bornes'] : erreur dans la somme des énergies des pièces et ajouts
```

Chaque application transforme un déplacement global en ports physiques
locaux. Les directions retenues de chaque instance restent privées, même
lorsqu'un objet réduction est réutilisé. La métrique globale est explicite ;
les énergies externes s'ajoutent aux énergies déjà présentes dans les pièces.
Aucune masse d'interface n'est ajoutée automatiquement.

Le bloc conservé global est assemblé **avant** sa résolution ; aucun Schur
local n'est inversé. Pour les applications normalisées `E_j` et bornes
locales `b_j`, l'enveloppe transportée est
`G = sum(b_j E_j.T @ E_j)`, et non la simple somme des `b_j`.
Sa plus grande valeur propre entre dans la marge globale. Les bornes et
les égalités de jonction restent évaluées en flottants. Le
[raisonnement et ses contre-exemples](ASSEMBLAGE_COMPLEMENTS_PREUVES.md)
détaillent la portée exacte de cette composition.

## Construction et preuves

Le C++ d'intervalles mesuré dans les campagnes historiques est repris sans
modification et lié à l'extension Rust via `cc`. Compiler les sources exige
un compilateur **C++17**, en plus de Rust. Les options strictes désactivent
fast-math et la contraction FMA. Chaque appel vérifie l'arrondi au plus
proche et le sous-flux graduel ; un environnement incompatible est refusé.
La roue ne lance aucun compilateur à l'import ou à l'appel.

Le pont Rust possède ses tableaux et vérifie leurs dimensions, pointeurs
CSR, indices, symétrie, permutation, finitude et budgets avant l'ABI C++.
Les preuves retournent leurs pivots en intervalles hexadécimaux et
l'identité des sources compilées. Ces identités ne prouvent pas à elles
seules la correction du compilateur ou du matériel.

La plateforme de livraison testée est Linux x86_64 / CPython 3.14 / GCC 13.3.
Les chemins de construction Clang et MSVC ne valent pas qualification de
ces plateformes. Les anciens scripts `ci/` et leurs archives sont conservés
pour permettre les rejeux historiques. Ils ne sont pas requis à l'exécution
du paquet installé.
