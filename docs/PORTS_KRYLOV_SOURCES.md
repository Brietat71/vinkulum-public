# Réduction par ports et Krylov : sources, antériorité et cible scientifique

Lecture ciblée du 7 septembre 2026. Cinq articles primaires et trois
documentations officielles ont été consultés. Ce complément au
[carnet de preuves](QUOTIENT_LOCAL_PREUVES.md) et à la
[recherche fondamentale](FONDEMENTS_MATHEMATIQUES_VINKULUM_2026.pdf)
porte sur les transferts linéaires conservatifs. Aucun code source de
solveur concurrent n'a été consulté.

**La cible est une réduction dont la dimension est choisie par l'erreur
du transfert mécanique, avec un coût de construction amortissable.**
Les sous-espaces de Krylov, les approximants de Padé, les réalisations
de Stieltjes et la synthèse modale par interfaces ont des antériorités
établies. Leur emploi ne constitue ni une nouveauté mathématique démontrée
ni une avance mesurée sur les solveurs multicorps.

## 1. Le problème utile au noyau

Après traitement des compatibilités et choix d'une métrique, considérons
un opérateur intérieur symétrique défini positif A et un bloc B de sollicitations
aux ports. Une partie du transfert dynamique s'écrit :

\[
F(z)=B^T(I-zA)^{-1}B,\qquad
0\le z<1/\lambda_{\max}(A).
\]

Dans la normalisation mécanique correspondante, z représente le carré
de la pulsation. Le domaine ci-dessus précède le premier pôle intérieur.
L'identification de A et B avec les matrices mécaniques et la composition
avec les contraintes doivent être dérivées pour chaque formulation ; une
borne sur F seule ne borne pas toutes les sorties du mécanisme assemblé.

Un petit sous-espace peut approcher ce transfert sans calculer tous les
modes intérieurs. Il faut comparer son coût à une méthode qui ne calcule
elle aussi que les modes nécessaires. Le nombre de ports, les résolutions
avec les matrices originales, la réorthogonalisation, le stockage des
bases et le nombre de configurations réutilisables comptent dans ce coût.
Un port de grande dimension ou une fréquence proche d'un pôle peut
annuler l'intérêt de la réduction.

Le transfert actuellement visé emploie **un contrôle par résidu général**.
Les articles sur Gauss–Radau précisent une possibilité d'affinement et
ses hypothèses ; ils ne décrivent pas un algorithme Radau livré dans
Vinkulum.

## 2. Cinq sources primaires et ce qu'elles établissent

### S1 — Encadrement matriciel par Gauss et Gauss–Radau

Jörn Zimmerling, Vladimir Druskin, Valeria Simoncini, *Monotonicity,
Bounds and Acceleration of Block Gauss and Gauss–Radau Quadrature for
Computing Bᵀφ(A)B*. **Journal of Scientific Computing 103, article 5,
2025**, publication le 12 février 2025.
[Texte intégral et DOI 10.1007/s10915-025-02799-z](https://link.springer.com/article/10.1007/s10915-025-02799-z).

**Lecture :** sections 2–5, théorème 1, proposition 3, remarques 1–3,
note 1 et périmètre de l'appendice.

Pour A définie positive et s réel positif, les approximations de
Bᵀ(A+sI)⁻¹B par blocs Gauss et Radau donnent des bornes monotones en
**ordre de Loewner sur la matrice complète**. Leur différence fournit
un majorant calculable. Les deux approximations réutilisent les
coefficients de la récurrence de Lanczos.

**Limites :** blocs de rang plein, absence de déflation et de rupture
de récurrence. L'argument sur l'arrondi contient une extrapolation aux
blocs confirmée expérimentalement. Les essais à fréquence imaginaire
ne constituent pas un théorème d'encadrement dans ce domaine. La
justification asymptotique de la moyenne en appendice est scalaire,
avec hypothèses de spectre continu.

### S2 — Erreur énergétique du gradient conjugué par blocs

Gérard Meurant, Petr Tichý, *Error Norm Estimates for the Block Conjugate
Gradient Algorithms*. **SIAM Journal on Matrix Analysis and Applications
46(4), 2025**.
[Publication, DOI 10.1137/25M1735408](https://epubs.siam.org/doi/10.1137/25M1735408) ;
[manuscrit auteur lu, arXiv:2502.14979](https://arxiv.org/pdf/2502.14979).

**Lecture :** théorèmes 5.1 et 7.6, lemme 7.5, discussion de l'arrondi
après le théorème 7.6 et conclusion.

Le théorème 5.1 décompose la matrice d'erreur énergétique entre
itérations en une somme positive. Le théorème 7.6 établit le signe
du reste Radau à partir d'un nombre \(0<\mu<\lambda_{\min}(A)\),
sous hypothèses de rang plein et de terminaison exacte. Il donne
ainsi des bornes matricielles, dont les diagonales contrôlent les
erreurs énergétiques des différents seconds membres.

**Limites :** une estimation de Ritz ne remplace pas une sous-estimation
spectrale établie. Les résidus et coefficients par blocs doivent avoir
les rangs requis. La discussion en précision finie explique pourquoi
l'argument de terminaison exacte ne s'applique pas automatiquement.
La déflation adaptative n'est pas couverte par ces hypothèses.

### S3 — Bornes a posteriori pour fonctions matricielles par blocs

Qichen Xu, Tyler Chen, *A posteriori error bounds for the block-Lanczos
method for matrix function approximation*. **Numerical Algorithms 98,
903–927, 2025**, publication en ligne en 2024.
[Publication, DOI 10.1007/s11075-024-01819-7](https://doi.org/10.1007/s11075-024-01819-7) ;
[manuscrit v4 du 14 avril 2024 lu](https://arxiv.org/pdf/2211.15643v4).

**Lecture :** note 1, théorème 2.1, section 2.1 et équation (11),
sections 2.2 et 3.5 sur l'arrondi.

Le théorème 2.1 borne l'erreur de f(H)V par une erreur de résolution
décalée multipliée par une intégrale de contour. L'équation (11) traite
V*f(H)V en norme opérateur, avec un terme quadratique en résidu. Les
quantités du contour s'évaluent sur les petites matrices réduites.

**Hypothèses et limites :** H hermitien, ensemble contenant son spectre
connu, analyticité de f dans les contours et facteurs réduits inversibles.
La note 1 suppose des espaces Krylov non dégénérés. La borne est une
borne en norme, pas une paire monotone de matrices en ordre de Loewner.
L'intégrale numérique et les erreurs d'arrondi doivent être contrôlées
pour obtenir une certification machine.

### S4 — Antériorité de la quadrature matricielle du résolvant

C. Fenu, D. Martin, L. Reichel, G. Rodriguez, *Block Gauss and Anti-Gauss
Quadrature with Application to Networks*. **SIAM Journal on Matrix
Analysis and Applications 34(4), 1655–1684, 2013**,
DOI 10.1137/120886261.
[Texte intégral auteur lu](https://www.math.kent.edu/~reichel/publications/blkagauss.pdf).

**Lecture :** sections 3.1, 5.1 et 6.1, théorème 5, équations
(3.7)–(3.10), conditions de l'estimateur au début de la section 8.

Le théorème 5 établit l'exactitude de Gauss matriciel jusqu'au degré
2N−1. Les équations (3.8)–(3.10) caractérisent l'opposition des erreurs
Gauss/anti-Gauss sur les polynômes jusqu'au degré 2N+1 et l'exactitude
de leur moyenne. La fonction 1/(1−cx) est explicitement traitée.

**Limites :** les encadrements discutés sont entrée par entrée, sous
conditions sur le développement de la fonction. Une différence
Gauss/anti-Gauss ne devient donc pas automatiquement un certificat
matriciel universel. L'antériorité concerne directement le calcul d'un
résolvant réduit sans diagonalisation de l'opérateur original.

### S5 — Sous-domaines, ports et réalisations de Stieltjes

Vladimir Druskin, Alexander V. Mamonov, Mikhail Zaslavsky, *Multiscale
S-Fraction Reduced-Order Models for Massive Wavefield Simulations*.
**Multiscale Modeling & Simulation 15(1), 445–475, 2017**.
[Publication, DOI 10.1137/16M1072103](https://doi.org/10.1137/16M1072103) ;
[manuscrit v2 du 16 octobre 2016 lu](https://arxiv.org/pdf/1604.06750v2).

**Lecture :** sections 3–5, hypothèses 1–2, lemme 4.1, proposition
8.1, remarques de coût et introduction de l'appendice C.

La méthode réduit les transferts Neumann–Dirichlet de sous-domaines
par Krylov rationnel, puis les transforme en réalisations équivalentes
par blocs tridiagonales. La proposition 8.1 préserve le caractère de
Stieltjes après assemblage, avec conséquences de stabilité et de
conservation pour la réalisation considérée.

**Limites :** dynamique linéaire conservative, masse positive,
entrées/sorties sur le squelette, rangs Krylov pleins et hypothèse de
degré de McMillan maximal. Les résolutions locales et leur factorisation
ont un coût de préparation à amortir. La proposition ne donne pas une
borne bilatérale générale d'erreur du mécanisme. Ports, Krylov rationnel
et réalisation mécanique creuse ont donc déjà été composés avant ce
travail sur Vinkulum.

## 3. Ce qui se transfère et ce qui reste à démontrer

La connexion de S1 avec F(z) demande un changement de variable explicite.
Pour z>0, si l'on connaît un nombre
\(\lambda_{\max}(A)<\alpha<1/z\), poser
\(C=\alpha I-A\) et \(s=1/z-\alpha\) donne :

\[
C\succ0,\qquad s>0,\qquad
F(z)=z^{-1}B^T(C+sI)^{-1}B.
\]

C'est une déduction algébrique de notre part. Elle conserve l'espace
Krylov polynomial par translation et réflexion de l'opérateur.
Elle ne supprime ni le besoin d'un majorant spectral ni les conditions
de rang de S1. Le cas z=0 se calcule directement.

La voie par résidu utilise une identité plus générale. Si
H=I−zA est définie positive, pour toute approximation X, poser
R=B−HX et \(F_c=B^TX+X^TB-X^THX\). Alors :

\[
B^TH^{-1}B-F_c=R^TH^{-1}R.
\]

Sous des bornes établies \(\delta I\preceq H\preceq\Delta I\),
avec δ>0, l'inversion de l'ordre donne :

\[
F_c+\Delta^{-1}R^TR
\preceq B^TH^{-1}B
\preceq F_c+\delta^{-1}R^TR.
\]

Cette identité variationnelle, redérivée sans revendication de nouveauté,
ne requiert pas que X provienne d'une récurrence de rang plein. Elle
permet de contrôler une approximation après déflation ou résolution
réduite imparfaite. Elle ne valide pas une suppression de contrainte
mécanique : celle-ci concerne le problème physique original.

Les inégalités sont mathématiques en arithmétique exacte. Si δ, Δ,
les produits ou le résidu ne sont qu'estimés, leur résultat ne doit pas
être présenté comme un encadrement rigoureux en précision finie.
Il reste également à propager l'erreur locale dans le système assemblé,
en tenant compte de son conditionnement et des sorties demandées.

## 4. Références concurrentes à inclure

Les trois documents ci-dessous ont été relus le 7 septembre 2026.
Ils décrivent des fonctionnalités et des formulations, pas des mesures
comparatives exécutées dans ce travail.

| Source officielle | Capacité directement documentée | Conséquence pour la confrontation |
| --- | --- | --- |
| **S6 — Exudyn 1.11.0**, [Model order reduction and component mode synthesis](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/ModelOrderReductionAndComponentModeSynthesis.html) | FFRF avec réduction, modes propres creux, HCB combinant modes statiques et modes internes, interfaces RBE2/RBE3 ; factorisation creuse réutilisée pour construire les modes statiques | Comparer à HCB avec modes sélectifs et interfaces pertinentes. Réduire l'interface et éviter une inverse dense sont déjà des pratiques disponibles. |
| **S7 — MBDyn**, [FAQ, composants déformables](https://www.mbdyn.org/Documentation/FAQ.html) | Synthèse modale de composants généraliste par élément « modal », poutres géométriquement exactes | Inclure une référence CMS, au-delà d'une comparaison à un maillage flexible complet. La FAQ lue ne précise pas une base HCB particulière. |
| **S8 — Dassault Systèmes / SIMULIA**, [Flexible Body Simulation Modules](https://www.3ds.com/products/simulia/simpack/flexible-body) | Approches modales, réduction linéaire et non linéaire, SIMBEAM et couplage Abaqus | La réduction dynamique et non linéaire est déjà proposée. Cette page ne permet pas d'identifier toutes les bases, bornes d'erreur ou méthodes internes. |

La dénomination HCB est donc vérifiée explicitement ici pour Exudyn.
Pour MBDyn et Simpack, les sources choisies établissent la CMS ou les
approches modales sans justifier d'attribuer une même réalisation HCB
aux trois logiciels. Une fonctionnalité non décrite par ces pages
n'est pas démontrée absente.

## 5. Expérience susceptible de décider de l'intérêt

La question scientifique est : **à transferts de ports et tolérance
physique identiques, un espace enrichi par Krylov avec contrôle résiduel
réduit-il le coût total par rapport à des modes HCB bien sélectionnés ?**

Le protocole doit faire varier séparément la taille intérieure, le nombre
et le rang des ports, la proximité d'un pôle et le nombre de requêtes
réutilisant la réduction. Il doit compter construction, factorisations,
produits, orthogonalisation, mémoire, évaluation et reconstruction des
sorties. Une économie sur l'évaluation seule ne suffit pas.

Les contre-épreuves doivent inclure ports collinéaires ou presque
colinéaires, multiplicités modales, couplages très anisotropes, fréquences
proches de la limite et charges hors de l'espace utilisé pour construire
la réduction. Pour une nouvelle charge, il faut établir son appartenance
aux sollicitations couvertes ou enrichir la construction.

Une référence directe indépendante contrôle le transfert complet et ses
directions d'erreur sur les petits cas. Les grands cas nécessitent un
résidu dans les équations originales et des bornes spectrales établies.
Les témoins réduits doivent conserver les mêmes interfaces et les mêmes
sorties physiques. Les cas où la réduction coûte davantage restent
dans les résultats.

La réussite de cette expérience établirait une brique utile dans son
domaine linéaire conservatif. Le passage aux grandes rotations, aux
contacts, aux changements de topologie, aux trajectoires et aux gradients
exige des démonstrations et des validations distinctes. Aucun classement
généraliste d'Exudyn, MBDyn ou Simpack n'est déduit de cette recherche.
