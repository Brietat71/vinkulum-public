# Structure et pivots du KKT natif : sondes non certifiées

8 septembre 2026. Analyse des modèles natifs figés à 32, 128 et 512
poutres ; aucune source du dépôt n'a été modifiée. Le script lit D et M
dans `/tmp/vinkulum-confrontation-ports-0.10.0-corrigee`, et les directions
fraîches dans l'archive `docs/bancs/krylov-contraint-2026/essais` du clone
`/tmp/vinkulum-energie-native-travail`. Il n'inspecte aucun code concurrent.

**Toutes les signatures numériques ci-dessous sont des sondes binary64,
pas des certificats.** Le produit K = fl(DᵀD), les soustractions avec M et
les éliminations sont arrondis sans encadrement. Le futur certificat doit
viser le produit mathématique des valeurs stockées de D, le M original et
le B effectivement retenu, avec arrondis dirigés.

## Données et portée

Les sondes d'inertie utilisent le B exact stocké dans
`n{n}-f40-krylov_controle-passage0/bases.npz`, de taille n_i×1 avec
n_i = 186, 762, 3066. Le rapport contient les chemins, empreintes SHA256,
versions et nombres de fils. Les huit bases fraîches de chaque taille
sont lues pour vérifier leur forme et leur corrélation absolue avec la
base sondée : celle-ci vaut 1 à quelques unités d'arrondi près. Cela
**ne remplace pas** un test d'inertie pour chaque B distinct.

Le script fixe BLAS, OpenMP, MKL et Rayon à un fil. La durée totale des
sondes et les temps de LDL avec reconstruction sont conservés à titre
indicatif seulement : pas de répétition statistique ni de comparaison
de solveurs, et aucun coût de certification dirigée n'est mesuré.

## Structure : l'ordre naturel est déjà très favorable

| Poutres | n_i | nnz fl(DᵀD), zéros retirés | nnz du motif DᵀD sans annulation | Produits pour assembler son triangle supérieur |
|---:|---:|---:|---:|---:|
| 32 | 186 | 786 | 910 | 980 |
| 128 | 762 | 3 282 | 3 790 | 4 052 |
| 512 | 3 066 | 13 266 | 15 310 | 16 340 |

Chaque ligne de D_i possède au plus quatre nonzéros. M_ii est diagonale
sur ces trois modèles. Le graphe physique comporte quatre composantes
de tailles N−1, 2(N−1), 2(N−1), N−1. Sa largeur naturelle est 10 ; RCM
appliqué au graphe physique réduit la largeur numérique à 4.

Pour le KKT Cγ = [[D_iᵀD_i−γM_ii,B],[Bᵀ,0]], B possède
186, 761 et 3 065 valeurs stockées non nulles. Certaines sont de l'ordre
de 10⁻²⁴. **Aucun seuil ne les supprime.** Les nombres suivants utilisent
le motif de DᵀD sans exploiter les annulations, donc conviennent à un
budget de préparation de l'assemblage dirigé :

| Poutres | nnz KKT structurel | nnz L naturel, diagonale comprise | Voisins futurs maximaux | Mises à jour de coefficients supérieurs |
|---:|---:|---:|---:|---:|
| 32 | 1 282 | 735 | 4 | 1 150 |
| 128 | 5 312 | 3 039 | 4 | 4 798 |
| 512 | 21 440 | 12 255 | 4 | 19 390 |

Il ne se crée que 0, 1 et 1 arête de remplissage. Le comptage d'une
mise à jour correspond à `a_ij ← a_ij − a_ik a_jk / pivot` pour i≤j ;
ce n'est pas un comptage des opérations Decimal ou de leurs allocations.
Les divisions peuvent être mutualisées par voisin. Avec l'ordre naturel,
le front possède au plus trois voisins physiques et le multiplicateur.
La grande largeur scalaire apparente du KKT, égale à n_i à cause de sa
bordure, ne doit donc pas conduire à allouer une bande pleine de largeur
n_i : un stockage par voisins actifs conserve cette petite complexité.

RCM physique ne réduit pratiquement pas ces comptes de factorisation.
Le RCM **global** testé place certes le multiplicateur presque à la fin,
mais augmente le remplissage : pour 512, 17 340 coefficients L et
42 781 mises à jour sur le motif numérique, contre 12 253 et 19 382 en
ordre naturel sur ce même motif. Il n'apporte pas d'avantage ici.

Une permutation générale raisonnable est P = diag(P_phys,I_s) :
réordonner le graphe physique par une méthode adaptée, puis laisser les
s multiplicateurs en bordure. Pour ces chaînes, P_phys = I est déjà un
bon choix. Sur un graphe physique de largeur d'élimination w, la bordure
ajoute au plus s voisins au front ; le coût des mises à jour reste de
l'ordre de n_i(w+s)² tant que les pivots n'imposent pas de grands retards.
La validité de la permutation ne dépend pas de la qualité de la méthode
de réordonnancement ; le remplissage et les budgets, eux, en dépendent.

## Signature aux seuils demandés

Pour γ = fl((2πf)²), les ordres naturel et RCM physique donnent :

| f, Hz | 32 poutres : inertie Cγ | 128 poutres : inertie Cγ | 512 poutres : inertie Cγ |
|---:|:---:|:---:|:---:|
| 40 | (186,1,0) | (762,1,0) | (3066,1,0) |
| 50 | (186,1,0) | (762,1,0) | (3066,1,0) |
| 60 | (186,1,0) | (762,1,0) | (3066,1,0) |
| 70 | (186,1,0) | (762,1,0) | (3066,1,0) |
| 80 | (186,1,0) | (762,1,0) | (3066,1,0) |
| 100 | (185,2,0) | (761,2,0) | (3065,2,0) |

La convention est (positifs, négatifs, nuls). Une LDL dense avec pivotage
concorde aux six seuils pour 32 et 128. Une diagonalisation indépendante
de la matrice physique bandée indique un négatif jusqu'à 80 Hz, puis
deux à 100 Hz, pour les trois tailles. Le dernier pivot bordé est positif.
Cela rend **80 Hz plausible comme cible d'un certificat strict**, au-delà
du minorant de trace vers 65,7 Hz. Ce tableau ne prouve ni la position
exacte du seuil contraint ni l'absence de singularité entre deux sondes.

## Échelles et pivots

La sonde de base du rapport utilise une normalisation facultative de B
contre diag(K) ; elle est conservée comme contre-exemple de choix
d'échelle défavorable. Son facteur de bordure est
1/||diag(K)⁻¹ᐟ²B||. Il crée des entrées de Schur de l'ordre de 10¹⁰ et
une forte croissance intermédiaire, bien que les signatures concordent.
Ce résultat ne constitue pas une instabilité intrinsèque du modèle.

La section `scalings` compare les données brutes, l'équilibre physique
seul, et deux équilibres avec une échelle inertielle du multiplicateur.
Le choix le plus simple est de **laisser le multiplicateur inchangé** et,
si nécessaire, d'équilibrer les variables physiques par puissances de
deux proches de 1/sqrt(K_jj). À 512, en ordre naturel :

| f, Hz | Maximum absolu de L_ij | Croissance maximale des entrées actives | Minimum du rapport pivot absolu / maximum de sa ligne active | Dernier pivot |
|---:|---:|---:|---:|---:|
| 60 | 4,080 | 1,000 | 0,245 | 9,68906 × 10⁻⁶ |
| 80 | 2,441 | 1,306 | 0,410 | 4,67839 × 10⁻⁶ |

La reconstruction flottante LDLᵀ présente un résidu relatif en norme
infinie inférieur à 2 × 10⁻¹⁶ dans ces deux cas. Ce résidu calculé ne
borne pas ses propres arrondis et ne certifie rien. Les facteurs d'échelle
physiques vont de 2⁻¹⁷ à 2⁻⁸. Une congruence par de telles puissances
préserve l'inertie en arithmétique exacte ; son application machine doit
encore éviter débordements et sous-flux et préserver les entrées visées.

Les six seuils ne réclament aucun pivot 2×2 dans la LDL scalaire, et les
contrôles denses n'en choisissent aucun. Le changement de signe d'un
pivot physique est normal : Aγ est déjà indéfinie avant le test du
complément. Il faut compter les signes, sans imposer tous les pivots
physiques positifs. L'intervalle d'un pivot doit exclure zéro avant division.

## Quand A est singulière : retarder, puis pivoter 2×2

La non-singularité de A n'est pas une hypothèse du critère KKT. Exemple
exact : A = diag(0,1), B = (1,0)ᵀ. Le complément contraint est positif,
alors qu'une LDL scalaire commençant par la première variable refuse.
Éliminer d'abord la seconde variable laisse [[0,1],[1,0]], de signature
(1,1,0). Le KKT entier a bien la signature (2,1,0).

La stratégie générale peut chercher d'abord un pivot 2×2 physique local,
ou retarder une variable physique ambiguë jusqu'à la bordure. Pour un bloc
symétrique [[a,b],[b,c]], un encadrement du déterminant strictement négatif
établit directement un positif et un négatif. Le cas [[0,b],[b,0]], b≠0,
est particulièrement simple ; il n'exige pas de borne de valeurs propres.
Si le déterminant est positif, le signe d'une diagonale non nulle sépare
les cas définis positif et négatif. Les divisions par le déterminant et
les mises à jour du Schur doivent être encadrées elles aussi.

Pivoter trop tôt avec un multiplicateur dense peut rendre dense le graphe
restant : son élimination relie ses voisins. Il vaut mieux le réserver à
une petite bordure après élimination des variables physiques sûres.
Même si la nullité de A est au plus s lorsque A est coercive sur ker(Bᵀ),
cela ne borne pas automatiquement le nombre de retards d'une stratégie
de pivots donnée. Les budgets doivent donc couvrir taille de front,
nonzéros, opérations et retards, avec refus sûr si ces limites sont atteintes.

## Reproduire

```sh
/tmp/vinkulum-release-0.11.0-final/venv/bin/python /tmp/vinkulum-sonde-inertie-2026/sonde_inertie.py > /tmp/vinkulum-sonde-inertie-2026/execution.log
```

`rapport.json` est le résultat principal complet, y compris les échelles
défavorables. `echelles.json` provient de la première sonde ad hoc ; ses
cas sont repris et complétés dans `rapport.json`. Aucun fichier NPZ n'est
copié. La prochaine étape utile est de confronter ces signatures à un
certificat dirigé de DᵀD−γM ; les comptes creux ne préjugent pas du nombre
de chiffres requis ni du coût réel de cette preuve.
