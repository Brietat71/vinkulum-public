# Relèvement statique, facteurs d'énergie et précision des ports

État de lecture : **7 septembre 2026**. Recherche bornée à six travaux
primaires, dont deux prépublications récentes identifiées comme telles.
Les articles ont été consultés ; aucun code de solveur concurrent n'a été
consulté. Ce document établit des conditions mathématiques et des
antériorités. Il ne constitue ni un benchmark ni une certification du
prototype Vinkulum.

La piste concrète consiste à conserver un **facteur de l'énergie avant
formation de la raideur**, puis à éliminer les variables intérieures par
transformations orthogonales. Le complément de Schur statique devient le
Gram d'un résidu énergétique de port. Cela évite une soustraction
potentiellement catastrophique, sans supprimer le conditionnement du
problème physique ni garantir à lui seul une erreur relative proche de
l'arrondi machine.

## 1. Six sources et ce qu'elles permettent d'affirmer

### S1 — Équivalence exacte entre Schur et élimination du facteur

**Nikolaus Demmel, Christiane Sommer, Daniel Cremers, Vladyslav Usenko,
2021.** *Square Root Bundle Adjustment for Large-Scale Reconstruction*,
CVPR 2021. [Manuscrit auteur](https://www.usenko.net/pdf/demmel2021rootba.pdf).

Lecture : **§4.3–4.4, équations (16)–(25)**. Pour un jacobien partitionné
\([J_l,J_p]\), un QR du bloc éliminé \(J_l=Q_1R_1\), avec \(R_1\)
inversible et \([Q_1,Q_2]\) orthogonal, donne exactement

\[
H_{pp}-H_{pl}H_{ll}^{-1}H_{lp}
  =(Q_2^TJ_p)^T(Q_2^TJ_p).
\]

L'équivalence porte aussi sur le second membre réduit et la reconstruction.
Il s'agit d'une démonstration algébrique, pas d'un théorème numéroté de
stabilité flottante. Le cas amorti est traité par augmentation du jacobien.

**Apport :** antériorité directe du calcul d'un Schur en conservant sa
racine, transposable à une énergie mécanique quadratique. La structure
par observations du problème de reconstruction explique son organisation
particulière ; ses gains expérimentaux ne se transfèrent pas aux ports
mécaniques. Ni l'identité ni le principe QR ne sont une invention Vinkulum.

### S2 — QR stable et petites valeurs singulières relativement précises

**Zlatko Drmač, Krešimir Veselić, 2005/2008.** *New Fast and Accurate
Jacobi SVD Algorithm. I*. Lecture : [LAPACK Working Note 169,
2005](https://netlib.org/lapack/lawnspdf/lawn169.pdf) ; publication
[SIAM J. Matrix Anal. Appl., 29(4), 1322–1342,
2008](https://epubs.siam.org/doi/10.1137/050639193).

Les numéros suivants sont ceux du rapport lu. **Proposition 2.2** : QR
Householder/Givens admet une perturbation colonne par colonne
\(\|\delta A_{:j}\|_2\le\epsilon_{qr}\|A_{:j}\|_2\), avec facteur
dépendant des dimensions et de l'algorithme. Pour \(A\) de rang colonne
plein, poser \(A_c=A\operatorname{diag}(\|A_{:j}\|_2)^{-1}\).
**Proposition 2.4** : les valeurs singulières exactes du facteur triangulaire
calculé vérifient

\[
\max_i\frac{|\widetilde\sigma_i-\sigma_i|}{\sigma_i}
 \le\sqrt n\,\epsilon_{qr}\|A_c^\dagger\|_2.
\]

Le **corollaire 2.5** contrôle l'étape Jacobi sur le transposé de ce facteur,
sous \(\sqrt n\epsilon_{qr}\|A_c^\dagger\|_2<1\).

**Limites :** modèle flottant sans dépassement ni sous-dépassement ; une
bonne mise à l'échelle des colonnes doit rendre le facteur suffisamment
conditionné. Ce résultat n'accorde pas une grande précision relative à
toute petite valeur propre d'une raideur arbitraire. Former préalablement
le Gram peut avoir détruit l'information recherchée.

### S3 — Former un Gram pour orthogonaliser exige des conditions explicites

**Takeshi Fukaya, Ramaseshan Kannan, Yuji Nakatsukasa, Yusaku Yamamoto,
Yuka Yanagisawa, 2020.** *Shifted Cholesky QR for Computing the QR
Factorization of Ill-Conditioned Matrices*, SIAM J. Sci. Comput., 42(1),
A477–A503. [Publication](https://doi.org/10.1137/18M1218212),
[prépublication 2018 lue](https://arxiv.org/pdf/1809.11085).

**Théorème 4.1**, avec les hypothèses (3.1)–(3.3), (4.12) : pour
\(X\in\mathbb R^{m\times n}\) de rang plein, \(m\ge n\),
\(h=mn+n(n+1)\),

\[
6n^2u\kappa_2(X)<1,\quad mnu\le1/64,\quad n(n+1)u\le1/64,
\quad\kappa_2(X)\le\frac1{96hu},
\]

`shiftedCholeskyQR3` satisfait

\[
\|\widehat Q^T\widehat Q-I\|_F\le6hu,\qquad
\frac{\|\widehat Q\widehat R-X\|_F}{\|X\|_2}\le15n^2u.
\]

Le décalage initial prescrit est \(11hu\|X\|_2^2\), suivi de deux étapes
CholeskyQR sans décalage. L'analyse étend la plage d'orthogonalisation
au-delà de la barrière usuelle de CholeskyQR autour de \(u^{-1/2}\), sous
ces restrictions dimensionnelles.

**Limite :** c'est un décalage interne à une factorisation QR, pas une
autorisation de modifier la raideur physique. Rang déficient, précision
relative du Schur et conservation mécanique nécessitent leurs propres
analyses. Householder reste une référence distincte pour éviter le Gram.

### S4 — Erreur arrière d'une élimination par blocs et croissance

**Neil Lindquist, Piotr Luszczek, Jack Dongarra, 2025.** *The Stability of
Block Eliminations and Additive Modifications*. [Prépublication arXiv
2509.07305](https://arxiv.org/pdf/2509.07305), version du 9 septembre 2025 ;
aucune validation éditoriale supposée ici.

Lecture : **§2.2, théorème 2.3 et corollaire 2.5**. Avec les erreurs locales
des factorisations, résolutions et mises à jour satisfaisant (2.4)–(2.7),
et les compléments de Schur effectués dans l'ordre considéré, la
factorisation calculée vérifie

\[
\|A-\widehat L\widehat R\|_\alpha
 \le C_AuP_\alpha\|A\|_\alpha
       +C_{LU}u\|\widehat L\|_\alpha\|\widehat R\|_\alpha.
\]

Les constantes cumulent les erreurs locales et \(P_\alpha\) mesure la
croissance des compléments successifs. Le corollaire ajoute les erreurs
de résolution au contrôle arrière du système complet. Le travail prolonge
l'analyse classique de l'élimination par blocs.

**Limite :** une erreur arrière petite relativement à \(A\) n'est pas une
erreur avant petite relativement à un petit Schur. Ce théorème ne certifie
pas automatiquement une chaîne d'assemblages différente de ses hypothèses.
Les modifications additives étudiées ailleurs dans l'article ne doivent
pas devenir une régularisation physique silencieuse de Vinkulum.

### S5 — QR creux : rang numérique et coût d'une suppression

**Leslie V. Foster, Timothy A. Davis, 2013.** *Algorithm 933: Reliable
Calculation of Numerical Rank, Null Space Bases, Pseudoinverse Solutions,
and Basic Solutions using SuiteSparseQR*, ACM TOMS, 40(1), article 7.
[Manuscrit auteur lu](https://people.engr.tamu.edu/davis/publications_files/Reliable_Calculation_of_Numerical_Rank_Null_Space_Bases_Pseudoinverse_Solutions_and_Basic_Solutions_using_SuiteSparseQR.pdf),
[DOI](https://doi.org/10.1145/2513109.2513116).

Lecture : **§2.1, théorème 5**. Lorsqu'une triangularisation supprime des
diagonales de valeur absolue au plus \(\tau\), le vecteur \(w\) des valeurs
supprimées donne, pour la perturbation correspondante \(E_1\),

\[
\|E_1\|_F=\|w\|_2,\qquad
\|E_1\|_2\le\sqrt{n-\ell}\,\tau.
\]

L'identité repose sur les transformations orthogonales du raisonnement ;
elle ne comprend pas à elle seule toutes leurs erreurs d'arrondi. Le
papier utilise des réflecteurs implicites et des estimations de valeurs
singulières pour signaler un rang incertain. Son introduction signale
explicitement une exception expérimentale à la fiabilité annoncée.

**Apport :** conserver des transformations creuses et comptabiliser une
compression est une antériorité établie. Un pivot petit ne prouve pas
qu'un degré de liberté mécanique est nul. Rang numérique, tolérance de
compression et erreur flottante doivent rester trois objets distincts.

### S6 — Formulation mixte mécanique et condensation récente

**Alexander Humer, Ivo Steinbrecher, Astrid Pechstein, 2026.** *Mixed Finite
Elements for Geometrically Exact Beams using Discontinuous Rotations and
Discrete Curvature*. [Prépublication arXiv 2605.04573v1, 6 mai
2026](https://arxiv.org/html/2605.04573v1).

Lecture : **§3, équations (30)–(31), puis §3.2**. Une transformation de
Legendre introduit le moment indépendant :

\[
\overline\varphi_\kappa(M,\kappa)
   =-\tfrac12M^TC_\kappa^{-1}M+\kappa^TM.
\]

La stationnarité en \(M\) retrouve la loi constitutive. Le dispositif
combine champs locaux discontinus et rotations nodales ; ces dernières
hybrident la continuité des moments. Les auteurs décrivent explicitement
la condensation des variables locales, laissant déplacements et rotations
nodales dans le système global. Ils précisent que, pour des poutres déjà
reliées seulement par leurs nœuds, l'avantage de condensation est limité
comparé aux éléments volumiques.

**Limites :** ce travail concerne une formulation de poutre et ses
propriétés mécaniques ; ce n'est pas un théorème d'erreur flottante ni une
preuve de meilleur conditionnement de tout système multicorps. Son statut
est celui d'une prépublication. Il offre une direction généraliste par
variables physiques et interfaces, sans justifier une supériorité actuelle.

## 2. Dérivation pour les ports linéaires

Les calculs suivants sont une spécialisation algébrique explicite du
principe de S1. On suppose un facteur réel
\(D\in\mathbb R^{m\times(n_I+n_S)}\), de rang colonne plein, tel que
\(K=D^TD\). Les colonnes sont partitionnées
\(D=[D_I,D_S]\) entre intérieur et ports. Une énergie élémentaire positive
\(\tfrac12(B_ex_e)^TC_e(B_ex_e)\) fournit par exemple des lignes
\(C_e^{1/2}B_e\), assemblées avec les applications locales/globales.
Ce facteur doit être conservé au niveau des données constitutives ; une
Cholesky de la raideur déjà arrondie ne peut récupérer l'information perdue.

Un QR complet de \(D_I\), appliqué implicitement aux colonnes de port,
donne

\[
Q^T[D_I,D_S]
  =\begin{bmatrix}R_I&C\\0&F\end{bmatrix},\qquad R_I\text{ inversible}.
\]

Pour un déplacement de port \(s\), minimiser l'énergie sur l'intérieur
annule le premier bloc :

\[
x_I^*(s)=-R_I^{-1}Cs=:Ts,\qquad
\min_{x_I}\|D_Ix_I+D_Ss\|_2^2=\|Fs\|_2^2.
\]

Il en résulte exactement

\[
S_0=K_{SS}-K_{SI}K_{II}^{-1}K_{IS}=F^TF.
\]

L'application pratique utile est de conserver \(F\), éventuellement
compressé par un second QR en un facteur triangulaire de port. Il n'est
pas nécessaire de former \(Q\) complet. Un QR creux peut réduire le
stockage, mais son remplissage dépend du graphe et de l'ordre d'élimination ;
le seul fait que \(D\) soit creux ne garantit pas un coût linéaire.
Une suppression de rang ajoute une perturbation, à borner séparément
selon S5.

### Une perte d'information que cette représentation peut éviter

Considérons, en arithmétique exacte,

\[
D_\varepsilon=\begin{bmatrix}1&1\\0&\varepsilon\end{bmatrix},\qquad
K_\varepsilon=\begin{bmatrix}1&1\\1&1+\varepsilon^2\end{bmatrix}.
\]

En éliminant la première colonne, le vrai Schur est
\(S_0=\varepsilon^2\). Si \(1+\varepsilon^2\) s'arrondit à \(1\), la
raideur stockée donne un Schur nul. Le facteur contient encore
\(F=\varepsilon\), tant que cette donnée reste représentable. Cet exemple
algébrique isole une perte due à la représentation, sans constituer une
mesure du prototype. Une perturbation absolue de \(D\) de taille \(u\)
peut néanmoins dominer \(\varepsilon\) lorsqu'il est lui-même trop petit.

### Pourquoi le conditionnement ne disparaît pas

Pour le facteur de rang plein,
\(\kappa_2(K)=\kappa_2(D)^2\). Éviter le Gram évite donc d'exposer certaines
opérations à ce carré, mais ne transforme pas la sensibilité intrinsèque
de \(Kx=f\) pour tout second membre en \(\kappa_2(D)\).

Une autre dérivation donne la sensibilité directe du Schur. Posons
\(X=[T^T,I]^T\), c'est-à-dire \(X=\begin{bmatrix}T\\I\end{bmatrix}\).
Alors \(S_0=X^TKX\), et la stationnarité intérieure annule les termes en
\(dX\). La différentielle est

\[
dS_0=X^T(dK)X,\qquad
\|dS_0\|_2\le\|X\|_2^2\|dK\|_2.
\]

Le facteur \(\|X\|_2^2\|K\|_2/\|S_0\|_2\) expose ainsi une possible
amplification relative. Une résolution intérieure avec petite erreur
arrière relativement à \(K\) ne garantit pas un petit défaut relativement
à \(S_0\), conformément à la distinction rendue nécessaire par S4.

Pour un facteur de port perturbé de \(\Delta F\), l'identité propre à ce
facteur donne

\[
\|(F+\Delta F)^T(F+\Delta F)-F^TF\|_2
 \le2\|F\|_2\|\Delta F\|_2+\|\Delta F\|_2^2.
\]

Cette borne suppose que \(\Delta F\) est déjà contrôlé ; elle ne constitue
pas une borne sur sa production par QR. Une erreur arrière petite dans
\(D\) peut donner une erreur relative importante dans \(F\) si les
colonnes de port sont presque contenues dans l'espace intérieur. Les
angles entre sous-espaces et le conditionnement de \(R_I\) interviennent.
Il faut également ajouter l'arrondi de la formation finale de \(F^TF\),
si elle est effectuée.

## 3. Passage exact à la dynamique et portée généraliste

Le changement de coordonnées statique est

\[
x=P\begin{bmatrix}y\\s\end{bmatrix},\qquad
P=\begin{bmatrix}I&T\\0&I\end{bmatrix},\qquad
P^TKP=\begin{bmatrix}K_{II}&0\\0&S_0\end{bmatrix}.
\]

La masse doit subir **la même congruence** : \(\widetilde M=P^TMP\).
Même lorsque \(M_{IS}=0\), le nouveau bloc \(\widetilde M_{IS}\) n'est
généralement pas nul. Pour \(z=\omega^2\) et
\(K_{II}-z\widetilde M_{II}\) inversible, le Schur dynamique exact est

\[
S(z)=S_0-z\widetilde M_{SS}
 -z^2\widetilde M_{SI}
  (K_{II}-z\widetilde M_{II})^{-1}\widetilde M_{IS}.
\]

Cette écriture est une déduction par blocs, pas une nouvelle méthode de
réduction. Elle sépare l'énergie statique de la correction dynamique et
définit un objet compatible avec une approximation de résolvante par
Krylov. Les preuves de troncature de cette approximation restent distinctes
des erreurs de calcul de \(T\), de \(F\), de la masse transformée et des
produits réduits. Des forces intérieures exigent aussi la transformation
et la condensation cohérentes du second membre.

Pour les basses valeurs propres, si \(M=LL^T\), alors

\[
L^{-1}KL^{-T}=(DL^{-T})^T(DL^{-T}),\qquad
\lambda_j=\sigma_j(DL^{-T})^2.
\]

Le résultat de S2 devient pertinent si les données du facteur normalisé
en masse et leur mise à l'échelle respectent ses hypothèses. La précision
relative des valeurs singulières peut alors protéger celle des valeurs
propres. Le calcul de \(L\), son application et les écarts spectraux pour
les vecteurs propres restent à analyser. Une SVD seulement stable en
norme ne garantit pas cette propriété pour les plus petites valeurs.

La formulation à variables physiques indépendantes de S6 est une seconde
voie pour conserver l'information avant élimination. Déjà dans le cas
linéaire, \(e=Dx\), \(D^Te=f\) remplace le système de raideur par des
relations du premier ordre. Leur résolution mixte demande un choix de
métriques, de préconditionnement et de stabilité des espaces ; compter
moins de dérivées dans chaque bloc ne prouve pas une meilleure précision.

La vocation généraliste impose enfin des distinctions : une tangente
avec précontrainte ou flambement peut être indéfinie et ne possède alors
pas de représentation globale réelle \(D^TD\). Une formulation signée ou
mixte et une nouvelle analyse sont nécessaires. Les contraintes, les
modes rigides, le contact et les mises à jour géométriques ne sont pas
couverts par l'hypothèse SPD du cas initial.

## 4. Objectif scientifique pour Vinkulum

L'objectif défendable est de **préserver une réponse de port faible au
milieu de contributions de raideur fortes**, en transportant l'information
constitutive et les facteurs orthogonaux jusqu'à l'interface. Une
revendication technique devra établir conjointement :

- la fidélité du facteur local à l'énergie discrétisée et à ses unités ;
- une borne d'erreur du relèvement et du facteur de port, tenant compte
  du conditionnement, du rang et des éventuelles compressions ;
- la cohérence des masses, charges et champs reconstruits ;
- une séparation explicite entre erreur de modèle, erreur de réduction
  sur la bande et erreur d'arrondi ;
- le coût du remplissage, du stockage et des reconstructions pour les
  graphes mécaniques visés.

Le QR énergétique, le Schur, la précision relative Jacobi et la
condensation mixte ont des antériorités précises. L'apport éventuel de
Vinkulum se situerait dans leur composition contrôlée pour des ports
mécaniques généralistes, avec des critères de refus quand les hypothèses
échouent. Les six sources ne suffisent pas à démontrer cette composition
en arithmétique flottante, et ne permettent pas de proclamer un avantage
sur les solveurs multicorps existants.
