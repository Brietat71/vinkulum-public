# Relèvement, condensation en énergie et assemblage des ports

Carnet mathématique du 7 septembre 2026. Ces preuves prolongent les
[bornes des ports Krylov](PORTS_KRYLOV_PREUVES.md). Elles traitent le
relèvement statique, les couplages de masse, la condensation par QR et
l'assemblage. Les identités sont exactes en arithmétique exacte. Aucune
certification machine ni performance mesurée n'est revendiquée ici.

Le problème numérique visé est précis : soustraire deux grandes matrices
K_SS et un transfert statique presque égal peut perdre la précision du
petit Schur résultant. Une représentation par énergie permet de calculer
ce Schur statique comme un Gram résiduel. Les erreurs du modèle assemblé,
de la factorisation et des évaluations restent à distinguer.

## 1. Un relèvement arbitraire donne une identité exacte

Partitionnons les coordonnées en intérieur I et interface S. K et M sont
symétriques ; M peut avoir un bloc M_IS non nul. Supposons

\[
K_{II}\succeq\lambda_*M_{II}\succ0,\qquad
0\le z=\omega^2\le\Omega^2<\lambda_*.
\]

La matrice dynamique intérieure est A_z = K_II−zM_II. La matrice
dynamique complète K−zM n'est pas nécessairement définie positive.
Soit W une transformation inversible des p coordonnées d'interface et
soit Ψ une matrice intérieure de relèvement quelconque, même inexacte.
Notons

\[
J_I=\binom I0,\qquad L=\binom\Psi I W,\qquad
K_0=L^TKL,\quad M_0=L^TML,
\]

\[
R_0=(K_{II}\Psi+K_{IS})W,\qquad
C=(M_{II}\Psi+M_{IS})W.
\]

**Proposition 1.** Le Schur exact dans les coordonnées normalisées du
port est

\[
\boxed{
S(z)=K_0-zM_0-(R_0-zC)^TA_z^{-1}(R_0-zC).}
\]

Il est égal à Wᵀ[(K−zM)_SS−(K−zM)_SI A_z⁻¹(K−zM)_IS]W,
indépendamment du relèvement choisi.

**Preuve.** Le changement de coordonnées u = J_Iv+Ls est inversible
puisque W l'est. Dans ces coordonnées, la matrice dynamique vaut

\[
\begin{pmatrix}
A_z&R_0-zC\\
(R_0-zC)^T&K_0-zM_0
\end{pmatrix}.
\]

L'élimination de v donne la formule. Le même champ physique résulte
de l'élimination dans les coordonnées originales, ce qui prouve
l'indépendance vis-à-vis de Ψ.

En particulier, si K_IIΨ+K_IS = 0 exactement, alors R₀ = 0 et

\[
\boxed{S(z)=K_0-zM_0-z^2C^TA_z^{-1}C.}
\]

Le transfert restant porte une échelle z². Cela supprime la soustraction
du grand transfert statique à K_SS ; cela ne supprime pas une éventuelle
cancellation physique entre les termes dynamiques près d'une résonance.

La formule de M₀ conserve intégralement les termes croisés :

\[
M_0=W^T(M_{SS}+M_{SI}\Psi+\Psi^TM_{IS}
+\Psi^TM_{II}\Psi)W.
\]

Poser M_IS = 0 sans hypothèse de modèle changerait le problème. Une
erreur de relèvement ne permet pas non plus de poser artificiellement
R₀ = 0 tout en conservant le même K₀ et le même M₀.

## 2. Un seul espace pour les deux semences

Posons μ = z/λ_*, ρ = Ω²/λ_* < 1 et M̂ = λ_*M_II. Le bloc des
semences et sa restriction fréquentielle sont

\[
\mathcal B=[R_0,\lambda_*C],\qquad
P(\mu)=\binom I{-\mu I},\qquad
B_\mu=\mathcal BP(\mu)=R_0-zC.
\]

Une base V commune peut être construite par résolutions dans K_II à
partir de toutes les colonnes de 𝓑. En notant H(μ) = 𝓑ᵀA_z⁻¹𝓑,
on a S = K₀−zM₀−PᵀHP.

Supposons qu'un transfert réduit énergétique H̃ vérifie

\[
0\preceq H(\mu)-\widetilde H(\mu)\preceq bI
\quad(0\le\mu\le\rho).
\]

La congruence conserve l'ordre et PᵀP = (1+μ²)I. Par conséquent

\[
\boxed{
0\preceq S_V-S
=P^T(H-\widetilde H)P
\preceq(1+\rho^2)bI.}
\]

Le facteur 1+ρ² est une majoration exacte pour cette norme et ce bloc.
Il n'est pas une marge empirique. Si les coordonnées des deux semences
ont été redimensionnées séparément, P et sa norme doivent être modifiés
en conséquence.

Avec Y(μ) les coefficients réduits associés à 𝓑, la reconstruction est

\[
u_S=Ws,\qquad
u_I=\Psi Ws-VY(\mu)P(\mu)s.
\]

Il faut utiliser les mêmes coefficients et les mêmes matrices physiques
pour le transfert, la reconstruction et les résidus.

En particulier, avec U_c = L−J_IVYP, le transfert énergétique corrigé
donne exactement S_V = U_cᵀ(K−zM)U_c. L'erreur de cette énergie est
le carré du résidu intérieur dans la métrique A_z⁻¹. Cette propriété
utilise la coercivité intérieure, sans exiger celle de toute la matrice
dynamique assemblée.

### Une enveloppe pondérée peut être moins pessimiste

Considérons les matrices réellement utilisées

\[
T=V^T\widehat MV,\qquad
[D_0,D_1]=V^T\mathcal B,\qquad
Y=(I-\mu T)^{-1}[D_0,D_1],
\]

et le transfert énergétique corrigé utilisant le vrai Gram VᵀK_IIV.
L'identité résiduelle reste valable sans supposer ce Gram égal à I.
Définir

\[
E_0=R_0-K_{II}VD_0,\qquad
E_1=\lambda_*C-K_{II}VD_1,\qquad
F=\widehat MV-K_{II}VT.
\]

Le résidu après congruence est exactement

\[
\mathcal R(\mu)=E_0-\mu E_1
+\mu F(I-\mu T)^{-1}(D_0-\mu D_1).
\]

Pour q ≥ 1, développer la résolvante jusqu'à q−2 comme dans le carnet
Krylov. Avec t = ‖T‖₂, ρt < 1 et la norme
‖X‖_(K⁻¹,2) = ‖K_II⁻¹/²X‖₂, on obtient

\[
\begin{split}
\mathcal E_q={}&\|E_0\|_{K^{-1},2}
+\rho\|E_1\|_{K^{-1},2}\\
&+\sum_{j=0}^{q-2}\left[
\rho^{j+1}\|FT^jD_0\|_{K^{-1},2}
+\rho^{j+2}\|FT^jD_1\|_{K^{-1},2}\right]\\
&+\frac{\rho^q\|FT^{q-1}\|_{K^{-1},2}
(\|D_0\|_2+\rho\|D_1\|_2)}{1-\rho t}.
\end{split}
\]

La somme est vide pour q = 1. L'identité d'erreur en énergie donne

\[
\boxed{\sup_{0\le\mu\le\rho}\|S_V-S\|_2
\le\frac{\mathcal E_q^2}{1-\rho}.}
\]

Cette enveloppe conserve séparément le défaut statique et la semence
massique, multipliée par μ. Elle peut être plus fine que l'enveloppe du
bloc entier suivie du facteur 1+ρ² ; aucune domination universelle entre
ces deux évaluations n'est supposée. Prendre le minimum de deux bornes
établies reste valide.

Pour un relèvement statique exact et une base Krylov exacte de q blocs
sur λ_*C, le transfert de C conserve 2q moments. Le facteur μ² ajoute
deux puissances : S_V−S = O(μ^(2q+2)) = O(ω^(4q+4)) à modèle fixé.
Comme auparavant, cet ordre local ne prouve pas une convergence uniforme
sur tous les maillages et toutes les bandes.

## 3. Condenser la racine carrée de l'énergie

Supposons maintenant que la raideur complète du modèle considéré possède
une représentation

\[
K=D^TD.
\]

Elle est donc positive semi-définie. Cette hypothèse convient notamment
à une assemblée d'éléments linéaires coercifs sans précontrainte
déstabilisante ; elle ne vaut pas pour toute tangente multicorps.
Les lignes de D peuvent être les déformations élémentaires pondérées
par les racines des rigidités. Posons D_I = DJ_I et D_P = D_S W.

Une élimination orthogonale des seules colonnes intérieures donne

\[
Q^T[D_I,D_P]=
\begin{pmatrix}R&P\\0&E\end{pmatrix},\qquad Q^TQ=I.
\]

La matrice R est carrée triangulaire inversible puisque K_II ≻ 0.
Aucune base admissible globale n'est nécessaire pour écrire ou appliquer
ces facteurs. Les colonnes de port restent présentes pendant l'élimination.

**Proposition 2.** Le relèvement statique normalisé et son énergie sont

\[
\boxed{\Psi W=-R^{-1}P,\qquad K_0=E^TE,\qquad K_{II}=R^TR.}
\]

**Preuve.** Pour un port s et un intérieur x, l'énergie est

\[
\|D_Ix+D_Ps\|_2^2
=\|Rx+Ps\|_2^2+\|Es\|_2^2.
\]

Le premier terme s'annule pour x = −R⁻¹Ps. Le second est exactement
l'énergie minimale. Il fournit le Schur statique comme un Gram, sans
soustraire K_SS et un grand transfert. Enfin D_IᵀD_I = RᵀR.

Les résolutions statiques peuvent réutiliser R :
K_II⁻¹B = R⁻¹(R⁻ᵀB). Mieux, pour un résidu de forces F_r,

\[
\boxed{F_r^TK_{II}^{-1}F_r
=(R^{-T}F_r)^T(R^{-T}F_r).}
\]

Une seule résolution triangulaire transposée suffit au Gram résiduel.
Cette écriture évite le produit entre un résidu et un déplacement résolu,
et construit explicitement un Gram positif en arithmétique exacte.

### Le relèvement calculé peut encore être inexact

Pour une matrice stockée X ≈ −R⁻¹P, poser h = RX+P. Dans le modèle
défini par les facteurs QR, le relèvement [X;W] a exactement

\[
\boxed{K_0=h^Th+E^TE,\qquad R_0=R^Th.}
\]

Son erreur statique d'énergie est quadratique en h. Employer K₀ = EᵀE
et R₀ = 0 tout en utilisant X inexact dans M₀ et C mélangerait deux
relèvements. On peut conserver h, ou représenter le relèvement exact
implicitement par la résolution triangulaire et contrôler les erreurs
de chacune de ses applications.

### Ce que garantit la structure creuse

Un QR par rotations de Givens appliquées aux lignes peut conserver les
facteurs et les colonnes de port dans des structures creuses. Il faut
rapporter le nombre de rotations, de mises à jour et les coefficients
stockés : R peut subir du remplissage. Une complexité linéaire sur une
chaîne exige effectivement un nombre linéaire de lignes actives et une
largeur de front bornée, ports compris. Elle ne découle pas du seul nom
« QR creux ». Un budget dépassé doit arrêter ou reporter l'élimination ;
supprimer des coefficients pour tenir ce budget constitue une perturbation
du modèle qui doit être comptée.

### Mise à l'échelle des colonnes dans le prototype

Les formules précédentes décrivent le QR sans changement des coordonnées
intérieures. Le code utilise une diagonale positive D_e pour équilibrer
les colonnes, puis factorise [D_I D_e, D_P]. On a alors exactement

    K_II = D_e⁻¹ RᵀR D_e⁻¹,
    ΨW = −D_e R⁻¹P,
    FᵀK_II⁻¹F = (R⁻ᵀD_e F)ᵀ(R⁻ᵀD_e F).

Pour le relèvement physique stocké X_I, le défaut devient
h = R D_e⁻¹ X_I+P. Ses contributions sont K0 = EᵀE+hᵀh et
R0 = D_e⁻¹Rᵀh. La masse et les applications de port restent exprimées
dans les coordonnées physiques. Ces congruences expliquent les facteurs
de mise à l'échelle du code ; elles ne sont pas des approximations
supplémentaires. Les calculs flottants de D_e et de ses inverses doivent
toujours être inclus dans une éventuelle certification machine.

## 4. Efforts équilibrés : contrôler l'erreur sans Schur exact dense

L'identité suivante donne une autre lecture du Gram précédent. Pour un
résidu intérieur F_r, considérons toute matrice d'efforts τ telle que

\[
D_I^T\tau=F_r.
\]

**Proposition 3.** Alors

\[
\boxed{F_r^TK_{II}^{-1}F_r\preceq\tau^T\tau.}
\]

**Preuve.** L'effort minimal est τ_* = D_IK_II⁻¹F_r. Tout autre effort
admissible s'écrit τ = τ_*+τ₀ avec D_Iᵀτ₀ = 0. Les colonnes de τ_*
et de τ₀ sont orthogonales ; ainsi τᵀτ = τ_*ᵀτ_*+τ₀ᵀτ₀.

Dans la factorisation QR exacte, τ = Q[R⁻ᵀF_r;0] est cet effort
minimal. Sa norme se calcule sans construire Q. Sur une chaîne ou un
arbre, une solution d'équilibre peut aussi se construire par accumulation
des efforts, indépendamment d'un solveur dense.

Si l'équilibre n'est qu'approximatif, définir δ = F_r−D_Iᵀτ. Lorsque
K_II ≽ κN avec N ≻ 0 et κ > 0 établis, on a encore

\[
\|F_r\|_{K^{-1},2}
\le\|\tau\|_2+\frac{\|N^{-1/2}\delta\|_2}{\sqrt\kappa}.
\]

En effet, D_IK_II⁻¹D_Iᵀ est un projecteur orthogonal, donc
‖D_Iᵀτ‖_(K⁻¹,2) ≤ ‖τ‖₂ ; ajouter la correction δ donne la formule.
Le carré de cette majoration, divisé par 1−ρ, contrôle l'erreur du Schur
énergétique. Le défaut d'équilibre et les erreurs d'évaluation doivent
rester présents ; un effort presque équilibré n'est pas automatiquement
un certificat.

Cette approche permet d'auditer une approximation à partir des seuls
champs de port et résidus, sans résoudre une référence dense de précision
supérieure pour tout l'intérieur. Elle ne fournit une borne machine que
si les constantes et les erreurs du petit audit sont elles-mêmes établies.

## 5. Quelle énergie les nombres flottants représentent-ils ?

Une matrice stockée D_f peut être considérée comme une matrice exacte
de nombres binaires réels. Elle définit alors le modèle K_D = D_fᵀD_f,
positif semi-défini en arithmétique réelle. La matrice assemblée par
K_a = fl(D_fᵀD_f) est en général différente. Une tangente native calculée
par une autre suite d'opérations peut encore donner un troisième résultat.
Ces différences ne doivent pas être effacées en déclarant les matrices
« identiques à l'arrondi près » lorsque la tolérance visée est comparable
à leur effet mécanique.

De même, les facteurs QR calculés sont exacts pour un modèle factorisé
qui peut différer de D_f. Si un opérateur orthogonal Q_* et une perturbation
ΔD satisfont

\[
[D_I,D_P]+\Delta D
=Q_*\begin{pmatrix}R_f&P_f\\0&E_f\end{pmatrix},
\]

alors les identités de la section 3 sont exactes pour ce modèle perturbé.
Pour les appliquer au modèle D_f original, il faut contrôler ΔD. Un
produit de rotations calculées n'est pas supposé exactement orthogonal
sans justification ; leurs défauts d'orthogonalité et les suppressions
éventuelles de coefficients entrent dans cet audit.

Pour deux modèles D et D+ΔD, la perturbation de raideur est

\[
\Delta K=D^T\Delta D+\Delta D^TD+\Delta D^T\Delta D.
\]

Sur des champs de port X, une borne plus informative que la seule norme
globale de ΔK est

\[
\boxed{
\|X^T\Delta KX\|_2
\le2\|DX\|_2\|\Delta DX\|_2+\|\Delta DX\|_2^2.}
\]

Elle mesure la perturbation de l'énergie effectivement portée par les
champs. Elle ne suffit pas, seule, à comparer les Schur exacts de deux
modèles : leurs champs minimisants peuvent différer. Une comparaison
rigoureuse utilise soit les principes de minimum dans les deux modèles,
soit un contrôle résiduel dans un modèle de référence explicitement fixé.

### Utiliser les facteurs comme outil, garder le modèle de référence

Si D_f est la référence, un champ candidat X_c à trace imposée définit
son Schur énergétique par l'énergie de D_fX_c et par la masse de X_c.
Son défaut par rapport au Schur exact de D_f est donné par son résidu
intérieur dans **ce même modèle**. Les facteurs QR peuvent donc servir
au calcul du candidat et au préconditionnement sans devenir implicitement
une nouvelle définition de K.

Dans cette démarche, utiliser K₀ = E_fᵀE_f puis vérifier seulement les
résidus dynamiques avec D_f ne contrôle pas le défaut statique entre les
facteurs et D_f. Il faut également évaluer l'énergie du relèvement dans
D_f, ou encadrer l'écart entre cette énergie et E_fᵀE_f.

Une équivalence énergétique établie permet de réutiliser le Gram
triangulaire pour le modèle original. Si

\[
cR_f^TR_f\preceq K_{II},\qquad c>0,
\]

alors

\[
F_r^TK_{II}^{-1}F_r
\preceq\frac1c(R_f^{-T}F_r)^T(R_f^{-T}F_r).
\]

Par exemple, si D_I = Q_I R_f+ΔD_I avec Q_IᵀQ_I = I et
‖ΔD_I R_f⁻¹‖₂ ≤ η < 1, l'inégalité triangulaire donne
c = (1−η)². Une estimation non garantie de η ne certifie pas c.
Ajouter une marge arbitraire de 10⁻³ à une valeur de Ritz ne remplace
aucun de ces arguments.

### Audit limité aux opérations utiles

La précision requise pour un contrôle fin peut être concentrée sur les
produits locaux D_fX_c, les petits Gram, les résidus et quelques résolutions
triangulaires, plutôt que sur une diagonalisation dense globale. Par
exemple, si Ŷ approche Y = D_fX_c avec une borne établie
‖Y−Ŷ‖₂ ≤ δ_Y, alors

\[
\|Y^TY-\widehat Y^T\widehat Y\|_2
\le2\|\widehat Y\|_2\delta_Y+\delta_Y^2.
\]

Il faut encore ajouter l'erreur du produit calculé ŶᵀŶ. Sous le modèle
usuel d'opérations correctement arrondies, sans dépassement ni
sous-dépassement invalidant ce modèle, un produit scalaire de h termes
admet une majoration composante par composante avec
γ_h = hu/(1−hu), où u est l'unité d'arrondi. Appliquée à chaque ligne
locale, elle donne |fl(D_fX_c)−D_fX_c| ≤ γ_h|D_f||X_c| si h majore
la longueur de chaque produit. Le calcul de cette enveloppe doit lui
aussi être arrondi vers l'extérieur pour constituer une preuve machine.

La même démarche vaut pour une masse factorisée M = F_MᵀF_M, lorsqu'une
telle représentation locale est disponible. Les erreurs de K₀, zM₀,
du transfert et de leur soustraction finale s'additionnent dans une
borne absolue. Une sommation compensée ou une précision accrue sur ces
petits audits peut réduire leur coût d'erreur ; leur seul emploi n'est
pas une certification sans enveloppe associée.

Viser une erreur de 10⁻¹⁰ requiert donc un budget dans la même métrique :
erreur de modèle, erreur de réduction, erreur d'évaluation et erreur
d'assemblage. Leur somme doit être établie sous 10⁻¹⁰ pour une affirmation
certifiée. Un accord mesuré avec un Schur analytique valide un essai
numérique, pas tous les modèles. Près d'un zéro physique du Schur, une
tolérance absolue peut rester pertinente quand une tolérance relative
devient impossible.

### Enveloppe uniforme dans D d'entrée, sans équivalence avec le QR

Une autre construction évite de devoir établir une équivalence entre
K_QR et K = D_IᵀD_I. Elle utilise directement la coercivité physique
K ≽ λ_*M_II et la norme duale de masse. Soient X₀ un relèvement normalisé,
de trace W, et V la base intérieure effectivement retenue. Les
coefficients calculés ont la forme

\[
Y(\mu)=(I-\mu T)^{-1}(a_0+\mu a_1),\qquad
X_c=X_0-J_IVY.
\]

Dans le prototype, a₀ provient du couplage statique QR projeté et a₁
du couplage massique projeté avec le signe négatif. L'identité ci-dessous
reste valable pour ces matrices de coefficients quelconques. Posons

\[
R_0^D=D_I^TDX_0,\qquad
\widehat C=\lambda_*(MX_0)_I,\qquad
\widehat M=\lambda_*M_{II},
\]

\[
E_0=R_0^D-KVa_0,\qquad
E_1=-\widehat C-KVa_1,\qquad
F=\widehat MV-KVT.
\]

Le résidu intérieur dans le modèle D vaut exactement

\[
\mathcal R_D=E_0+\mu E_1+\mu FY.
\]

En effet, Y−a₀ = μ(a₁+TY). Substituer cette identité dans
R₀^D−μĈ−(K−μM̂)VY donne la formule. Pour tenir compte de la
compensation entre les premiers coefficients, définir

\[
d=Ta_0+a_1,\qquad C_1=E_1+Fa_0.
\]

Comme Y = a₀+μ(I−μT)⁻¹d, on obtient

\[
\boxed{\mathcal R_D
=E_0+\mu C_1+\mu^2F(I-\mu T)^{-1}d.}
\]

Pour q ≥ 2, l'expansion finie exacte est

\[
\mathcal R_D=E_0+\mu C_1
+\sum_{j=0}^{q-3}\mu^{j+2}FT^jd
+\mu^qFT^{q-2}(I-\mu T)^{-1}d.
\]

La somme est vide pour q = 2. Avec t = ‖T‖₂, ρt < 1 et
‖Z‖_(M⁻¹,2) = ‖M_II⁻¹/²Z‖₂, une enveloppe uniforme est

\[
\begin{split}
\mathcal R_{D,q}^{\max}={}&\|E_0\|_{M^{-1},2}
+\rho\|C_1\|_{M^{-1},2}\\
&+\sum_{j=0}^{q-3}\rho^{j+2}\|FT^jd\|_{M^{-1},2}\\
&+\frac{\rho^q\|FT^{q-2}\|_{M^{-1},2}\|d\|_2}{1-\rho t}.
\end{split}
\]

La coercivité donne A_z ≽ λ_*(1−ρ)M_II. Si S_c,D est le Schur
énergétique du candidat X_c évalué dans D, alors

\[
0\preceq S_{c,D}-S_{\mathrm{exact},D}
=\mathcal R_D^TA_z^{-1}\mathcal R_D,
\qquad
\boxed{\sup\|S_{c,D}-S_{\mathrm{exact},D}\|_2
\le\frac{(\mathcal R_{D,q}^{\max})^2}{\lambda_*(1-\rho)}.}
\]

Aucune annulation Krylov ni égalité entre K_QR et D_IᵀD_I n'est
supposée. Les normes duales demandent des résolutions dans M_II, souvent
diagonale ou par petits blocs, au lieu de résolutions dans K. Une masse
générale peut demander une factorisation supplémentaire dont le coût doit
être inclus. Prendre le minimum des enveloppes pour plusieurs q est
valide en arithmétique exacte.

Il reste à comparer ce fonctionnel avec le Schur S_rapporté réellement
renvoyé. Notons ses matrices K₀^Q, M₀^Q, G_Q, T, a₀ et a₁ :

\[
S_{\mathrm{rapport\acute e}}
=K_0^Q-zM_0^Q
-2\operatorname{sym}[(a_0+\mu a_1)^TY]
+Y^T(G_Q-\mu T)Y.
\]

Définir les écarts de petites projections

\[
\begin{array}{ll}
\Delta K_0=(DX_0)^T(DX_0)-K_0^Q,&
\Delta\widehat M_0=\lambda_*[X_0^TMX_0-M_0^Q],\\
\Delta K_c=(D_IV)^TDX_0-a_0,&
\Delta M_c=V^T\widehat C+a_1,\\
\Delta K_V=(D_IV)^TD_IV-G_Q,&
\Delta M_V=V^T\widehat MV-T.
\end{array}
\]

Leur effet exact est

\[
S_{c,D}-S_{\mathrm{rapport\acute e}}
=\Delta K_0-\mu\Delta\widehat M_0
-2\operatorname{sym}[Y^T(\Delta K_c-\mu\Delta M_c)]
+Y^T(\Delta K_V-\mu\Delta M_V)Y.
\]

Avec Y_max = (‖a₀‖₂+ρ‖a₁‖₂)/(1−ρt), on en déduit

\[
\begin{split}
\delta_{\mathrm{fonctionnel}}={}&\|\Delta K_0\|_2
+\rho\|\Delta\widehat M_0\|_2\\
&+2Y_{\max}(\|\Delta K_c\|_2+\rho\|\Delta M_c\|_2)\\
&+Y_{\max}^2(\|\Delta K_V\|_2+\rho\|\Delta M_V\|_2).
\end{split}
\]

Quand les projections massiques sont exactement communes, les trois
écarts massiques s'annulent. Les conserver permet aussi de comptabiliser
des divergences entre deux suites de calcul des mêmes projections, tout
en gardant à contrôler les arrondis de cet audit lui-même.

Si ε_res majore l'erreur résiduelle précédente, le Schur renvoyé vérifie

\[
-\delta_{\mathrm{fonctionnel}}I
\preceq S_{\mathrm{rapport\acute e}}-S_{\mathrm{exact},D}
\preceq(\delta_{\mathrm{fonctionnel}}+\varepsilon_{\mathrm{res}})I.
\]

Sa norme d'erreur est donc majorée par la somme de ces deux termes.
Le signe positif de l'erreur ne se transporte pas automatiquement au
fonctionnel QR renvoyé. Cette enveloppe porte explicitement sur le modèle
D d'entrée, avec sa constante de coercivité. Son évaluation flottante
reste non certifiée tant que les produits D, les résolutions massiques,
les petites projections et la constante λ_* ne sont pas encadrés.

Cette preuve ne dépend pas de la méthode ayant produit V ou X₀. Un
préconditionneur LU du Gram assemblé K_a peut donc construire des
candidats par raffinement x ← x+K_a⁻¹[b−D_Iᵀ(D_Ix)], en évaluant
le résidu dans D d'entrée. Une baisse du résidu équilibré ou un nombre
fixe de corrections ne prouve pas une contraction ni la précision de
l'inverse. De même, τ = D_I solve(b) fournit seulement des coordonnées
duales approchées tant que le défaut D_Iᵀτ−b n'est pas encadré ; son
Gram ne doit pas être déclaré égal à bᵀK⁻¹b en machine. L'audit uniforme
indépendant dans D, avec la norme duale de masse et les écarts du
fonctionnel, reste valable pour n'importe quel candidat ainsi obtenu
et décide séparément si la tolérance estimée est atteinte. Un Gram qui
perd le rang lors de son assemblage peut rendre la LU inutilisable même
quand D_I garde son rang : ce refus ne doit pas être transformé en
acceptation implicite ni en bascule non annoncée vers le QR.

## 6. Assemblage et erreur de réponse

Pour chaque sous-structure s, soit E_s l'application des coordonnées
globales vers ses ports normalisés. Elle inclut les changements de base,
les orientations et les normalisations nécessaires. À configuration
fixée, le Schur assemblé est

\[
S_G(z)=A_{\mathrm{support}}(z)+\sum_sE_s^TS_s(z)E_s.
\]

A_support porte les coordonnées rigides retenues et les contributions
qui ne sont pas condensées. Chaque énergie élémentaire et chaque masse
doit être comptée une seule fois, ou selon une partition explicitement
définie ; dupliquer les masses de port changerait le modèle.

Si 0 ≼ Ŝ_s−S_s ≼ b_sI, et si les supports sont conservés exactement,

\[
\boxed{
0\preceq\widehat S_G-S_G
\preceq\sum_sb_sE_s^TE_s.}
\]

Ainsi, une borne de la norme de la petite matrice ou de l'opérateur
Σ_s b_sE_sᵀE_s contrôle l'erreur globale. La majoration plus grossière
Σ_s b_s‖E_s‖₂² reste valide. Imposer des contraintes linéaires globales
par restriction aux coordonnées compatibles conserve ces inégalités ;
cela ne fournit pas une paramétrisation des contraintes non linéaires.
Les arrondis d'assemblage doivent être ajoutés séparément et peuvent
perdre le signe positif du défaut calculé.

Si b_s borne seulement la norme d'erreur du Schur effectivement renvoyé
dans D d'entrée, comme dans l'audit précédent, l'encadrement assemblé
devient bilatéral :

\[
-\sum_sb_sE_s^TE_s
\preceq\widehat S_G-S_G
\preceq\sum_sb_sE_s^TE_s.
\]

La majoration de norme et l'argument de réponse ci-dessous restent
valides ; seule la conclusion de signe positif doit être retirée.

Pour une réponse harmonique, notons ΔS = Ŝ_G−S_G, et soit x̂ un candidat
avec résidu r = f−Ŝ_Gx̂. Si σ_min(S_G) ≥ β > 0 est établi, alors

\[
S_G(x-\widehat x)=r+\Delta S\widehat x,
\qquad
\boxed{\|x-\widehat x\|_2
\le\frac{\|r\|_2+\|\Delta S\|_2\|\widehat x\|_2}{\beta}.}
\]

Une constante β peut parfois être obtenue à partir du seul problème
assemblé réduit : si ‖ΔS‖₂ ≤ b et σ_min(Ŝ_G) > b sont établis, alors
β = σ_min(Ŝ_G)−b convient. Cette voie ne demande pas le Schur exact
de tous les intérieurs. Une valeur singulière seulement estimée fournit
une conclusion seulement conditionnelle. Des modes rigides non contraints
ou une résonance globale peuvent annuler β malgré la coercivité de tous
les intérieurs.

Toutes ces normes supposent des coordonnées physiques mises à l'échelle
de manière cohérente. Pour une métrique de déplacement différente,
les forces, opérateurs et bornes doivent être transportés ensemble.
Une petite erreur de Schur local n'autorise donc pas une affirmation de
précision de réponse sans contrôle de l'assemblage et de sa stabilité.
