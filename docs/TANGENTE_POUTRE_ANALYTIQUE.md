# Efforts analytiques et tangente des poutres

Livraison 0.8.2 du 7 septembre 2026. Les formules ci-dessous différencient les deux
énergies existantes de Vinkulum. Elles ne changent ni les inconnues nodales,
ni la loi constitutive, ni l'interpolation en une autre famille d'éléments.

## Dérivation par travail virtuel

Soient les repères matériels A et B, la corde d = x_B − x_A et la longueur
au repos L. On pose θ = log(AᵀB), Q = A exp([θ]/2), c = Qᵀd/L et κ = θ/L.
Le crochet [θ] représente le produit vectoriel par θ. Les perturbations
de rotation sont à gauche, dans le repère monde.

Pour l'option milieu, S = I. Pour l'option intégrée :

\[
S=I+a(q)[\theta]^2,\quad q=\theta^T\theta,\qquad
a(q)=\frac{1-(\sqrt q/2)/\sin(\sqrt q/2)}q.
\]

La déformation γ = Sc − e₁ donne l'énergie
U = L(γᵀCₙγ + κᵀCₘκ)/2. Les rigidités de cisaillement Cₙ sont celles
effectives de l'option intégrée, incluant la flexibilité interne condensée.
Avec n = Cₙγ, le résultat en translation est F = QSn : F_A = F, F_B = −F.

Le gradient partiel en θ, à Q et c fixés, vaut

\[
w=C_m\kappa+L\{a[c(\theta\cdot n)+n(\theta\cdot c)
-2\theta(n\cdot c)]
+2a'\theta[(\theta\cdot n)(\theta\cdot c)-q(n\cdot c)]\}.
\]

Pour l'option milieu, w = Cₘκ. La variation de θ est
δθ = H Aᵀ(δω_B − δω_A), avec H = J_l(θ)⁻¹. La variation de la rotation
moyenne utilise D = (I + exp([θ]/2))⁻¹ = I/2 − c_m[θ], où
c_m = tan(|θ|/4)/(2|θ|). Ces transports proviennent de la composition des
rotations et du travail conjugué ; ils ne se remplacent pas par une moyenne
arithmétique des gradients.

En posant g = d × F et

\[
k=A H^T w-c_m(A\theta)\times g,
\qquad M_A=g/2+k,\qquad M_B=g/2-k,
\]

on obtient les douze efforts du gradient à gauche −dU. Le bilan
M_A + M_B = d × F est satisfait par construction. Pour évaluer Hᵀw :

\[
H^Tw=w+\tfrac12\theta\times w+b\theta\times(\theta\times w),
\qquad b=\frac{1-(|\theta|/2)\cot(|\theta|/2)}q.
\]

La tangente est la dérivée de ce **champ d'efforts**, avec un passage de
duaux à douze directions. Ce n'est pas la Hessienne symétrique de l'énergie
dans une carte exponentielle unique : le bloc rotation/rotation peut être
non symétrique hors équilibre. L'assemblage utilise la convention opposée
K = −∂f/∂q, augmentée des termes de réactions des contraintes.

## Calcul flottant et limites

- Pour q < 1/4, huit coefficients de série évitent les annulations de a,
  a', b et c_m. La dérivée a' est évaluée par Horner simultanément avec a,
  donc elle est exactement celle du polynôme employé par l'énergie.
- L'origine des deux corps est commune et locale. La partie basse de la
  corde compensée est conservée séparément dans les déformations et dans
  les contributions linéaires aux moments.
- La coupure du logarithme principal à π reste exclue du gradient. Les
  efforts et la tangente y signalent une valeur non définie ; fournir un
  résultat fini artificiel ne rendrait pas le modèle différentiable.
- Les rigidités dérivées peuvent être signées. La même expression sert
  aux sensibilités constitutives, avec la règle de chaîne des coefficients
  effectifs de l'option intégrée.

## Contrôles indépendants

Le calcul de référence différencie directement le scalaire d'énergie.
Pour la tangente, les rotations successives exp(ε) exp(δ) R garantissent
que la dérivation emboîtée porte sur le même champ à gauche.

Les tests Rust couvrent les deux formulations, les rotations de zéro à
π − 10⁻⁴, les deux côtés du raccord des séries, le bilan des moments,
les repères matériels distincts, les unités de longueur 10⁻⁶ à 10⁶,
les coefficients signés et l'objectivité sous rotation et translation
communes jusqu'à une origine distante de 10¹² longueurs. Les contrôles
préexistants d'allongement inférieur à l'ULP et de modes rigides restent
appliqués. La coupure π est vérifiée séparément.

Les identités s'inscrivent dans le calcul variationnel sur groupes de Lie
employé notamment par [Sonneville, Cardona et Brüls, 2014](https://orbi.uliege.be/bitstream/2268/159471/1/pa_SonnevilleCardonaBruls2014_GeometricallyExactBeamFEOnSE3.pdf).
La dérivation présentée ici est propre à nos deux énergies ; elle ne
revendique ni une nouvelle théorie de poutre, ni une reproduction de
l'interpolation SE(3) de cet article.

## Mesures entre roues figées

La campagne comporte 51 configurations et trois répétitions par version,
soit 306 résolutions chronométrées. Les versions alternent, sur CPU 0,
avec un fil demandé aux bibliothèques. Chaque processus chauffe un modèle,
puis mesure une résolution neuve ; construction, import et calcul des
références indépendantes sont exclus. Les sorties détaillées des
échauffements ne sont pas conservées. Python 3.14.7, NumPy 2.5.3 et
SciPy 1.18.1 sont identiques dans les deux environnements isolés.

Les six rampes Princeton de cinquante paliers donnent :

| Formulation | Éléments | 0.8.1 (ms) | 0.8.2 (ms) | Temps divisé par |
|---|---:|---:|---:|---:|

| milieu | 8 | 75.741 | 31.178 | 2.43 |
| milieu | 20 | 228.835 | 115.660 | 1.98 |
| milieu | 60 | 611.377 | 280.493 | 2.18 |
| integree | 8 | 86.385 | 33.807 | 2.56 |
| integree | 20 | 253.685 | 123.208 | 2.06 |
| integree | 60 | 694.692 | 295.301 | 2.35 |

Les rampes milieu conservent 334 évaluations principales. Les rampes
intégrées passent de 334 à 335 ; elles demandent toujours une tentative
par palier. Le gain vient surtout du coût de chaque évaluation et tangente.

Les 49 contrôles précédemment réussis passent toujours, dont le refus
attendu du corps sans équilibre ; cela représente 48 équilibres obtenus.
L'écart maximal des positions convergées est 2,151 × 10⁻¹² m et celui des
composantes de rotation 6,916 × 10⁻¹². Les deux cas à petite échelle restent
non résolus au seuil resserré. Leurs refus coûtent 2,79 et 2,64 fois moins
cher, sans que cela les transforme en solutions acceptables.

Sur les médianes de cette campagne, les chaînettes varient de −0,6 % à
+0,2 %, les boucles de −1,4 % à +6,1 % et les contacts de moins de 0,1 %.
Trois répétitions ne suffisent pas à attribuer au changement de code les
petites variations de temps des familles sans poutre.

Le diagnostic distinct `k_c_m_z` mesure l'appel complet d'analyse, incluant
matrices globales et base admissible : ce n'est pas le coût isolé d'une
poutre. Le gain vaut 1,42–1,53 fois à 10 éléments, mais seulement 1,04–1,09
fois à 120 éléments. Les produits K par quatre vecteurs fixés diffèrent
de moins de 7 × 10⁻¹⁵ en norme relative globale ; ceux de C sont identiques.
Ce contrôle global ne remplace pas les tests locaux composante par
composante. Le diagnostic conserve trois mesures après un échauffement
par taille, dans quatre processus successifs ; son ordre est archivé.

[Résolutions brutes](bancs/tangente-poutre-0.8.2.json),
[bilan recalculé](bancs/tangente-poutre-0.8.2-bilan.json),
[diagnostic des matrices](bancs/tangente-poutre-0.8.2-matrices.json).
Pour rejouer le contrôle sans relancer les solveurs :

```bash
python ci/bilan_globalisation.py docs/bancs/tangente-poutre-0.8.2.json
```

Le processus de mesure rend volontairement un code non nul à cause des
deux échecs numériques. Le vérificateur d'archive exige leur présence et
leur restauration tout en vérifiant les autres cas. Aucun gain face à
MBDyn ou Simpack n'est déduit directement de ces mesures internes.
