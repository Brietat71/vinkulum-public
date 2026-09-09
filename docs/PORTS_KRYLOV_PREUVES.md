# Ports dynamiques : preuves et limites de la réduction par Krylov statique

Carnet mathématique du 7 septembre 2026. Les propositions ci-dessous sont
redérivées, sans revendication de nouveauté. Elles prolongent le
[carnet sur les quotients locaux](QUOTIENT_LOCAL_PREUVES.md), dans le domaine
plus restreint d'un intérieur linéaire conservatif et coercif. Ce document
ne rapporte aucun résultat chronométré.

La construction vise à préparer un transfert d'interface par résolutions
répétées avec une même raideur locale, puis à contrôler son approximation
sur une bande entière. Les preuves sont en arithmétique exacte. Leur
évaluation en nombres flottants constitue un diagnostic numérique tant que
les erreurs de calcul et la constante de coercivité ne sont pas encadrées.

## 1. Domaine et normalisation

Soient K et M deux matrices symétriques définies positives de taille n,
respectivement raideur et masse d'un intérieur. Les mouvements rigides ou
les contraintes restantes doivent avoir été traités avant de se placer dans
ce cadre. Soit B une matrice de ports de taille n×p, sans hypothèse de rang
plein. On considère le transfert

\[
H(z)=B^T(K-zM)^{-1}B,\qquad z=\omega^2.
\]

On suppose disposer d'une constante établie λ_* > 0 telle que

\[
K\succeq\lambda_*M.
\]

La bande physique est 0 ≤ ω ≤ Ω, avec Ω² < λ_*. Posons

\[
\widehat M=\lambda_*M,\qquad
\mu=z/\lambda_*,\qquad
\rho=\Omega^2/\lambda_*<1,\qquad
A_\mu=K-\mu\widehat M.
\]

Alors, sur toute la bande,

\[
(1-\rho)K\preceq(1-\mu)K\preceq A_\mu\preceq K.
\]

Nous écrirons H(μ) pour BᵀA_μ⁻¹B dans les démonstrations. Cette notation
change le paramètre, pas le transfert physique.

Toutes les inverses écrites ici désignent des opérateurs mathématiques :
leur mise en œuvre peut utiliser des résolutions et des factorisations.
Les racines K^(±1/2) servent aux preuves et aux définitions de normes ;
elles n'imposent pas la construction d'une matrice dense.

## 2. Blanchir le vrai Gram de la base

Soit V une matrice n×r de colonnes indépendantes. Elle peut provenir d'une
construction Krylov avec des blocs de tailles variables. Son Gram exact
est G = VᵀKV ≻ 0. Si G = LLᵀ est son Cholesky, définir

\[
\overline V=VL^{-T},\qquad
\overline V^TK\overline V=I_r.
\]

Les petites matrices et la solution réduite sont

\[
T=\overline V^T\widehat M\overline V,\qquad
D=\overline V^TB,\qquad
Y(\mu)=(I_r-\mu T)^{-1}D,\qquad X_V=\overline VY.
\]

Pour r > 0, 0 ≺ T ≼ I_r ; en particulier t = λ_max(T) ≤ 1. Le transfert
réduit est

\[
H_V(\mu)=D^T(I_r-\mu T)^{-1}D.
\]

Si r = 0, on définit H_V = 0 et t = 0 ; les termes impliquant la base sont
vides. Le résidu est alors simplement B. Le blanchiment suppose que les
colonnes retenues soient indépendantes. Un Gram singulier ne permet pas
de continuer cette formule sans décider quelles directions conserver.

En calcul flottant, l'identité du Gram ne doit pas être supposée à partir
du seul historique d'orthogonalisation : il faut former et contrôler le
Gram courant. Même après blanchiment numérique, son égalité à I reste
approximative. Une formulation équivalente conserve explicitement
VᵀKV dans le système réduit.

## 3. Identité d'erreur valable pour toute base

Définissons les défauts de semence et d'invariance

\[
E=B-K\overline VD,\qquad
F=\widehat M\overline V-K\overline VT.
\]

Ils satisfont exactement V̄ᵀE = 0 et V̄ᵀF = 0. Le résidu physique est

\[
R(\mu)=B-A_\mu X_V.
\]

**Proposition 1.** Sans supposer que B appartient à l'espace des semences
retenues, on a

\[
\boxed{R(\mu)=E+\mu F(I_r-\mu T)^{-1}D,}
\]

et

\[
\boxed{H(\mu)-H_V(\mu)=R(\mu)^TA_\mu^{-1}R(\mu)\succeq0.}
\]

**Preuve du résidu.** L'équation réduite donne Y−D = μTY. Ainsi

\[
B-(K-\mu\widehat M)\overline VY
=B-K\overline VD+\mu(\widehat M\overline V-K\overline VT)Y.
\]

**Preuve de l'identité d'erreur.** Poser X = A_μ⁻¹B. La condition de
Galerkin V̄ᵀA_μ(X−X_V) = 0 annule les termes croisés. Par conséquent

\[
(X-X_V)^TA_\mu(X-X_V)=H-H_V.
\]

Puis X−X_V = A_μ⁻¹R donne la formule annoncée. La preuve utilise seulement
la symétrie, la coercivité et la résolution exacte du problème réduit.
Elle ne nécessite ni Krylov exact, ni annulation de E, ni rang plein de B.

L'ordre des inverses donne également

\[
\boxed{
R^TK^{-1}R\preceq H-H_V
\preceq\frac{R^TK^{-1}R}{1-\mu}.}
\]

En particulier, le défaut de semence a un effet statique explicite :

\[
H(0)-H_V(0)=E^TK^{-1}E.
\]

Supprimer une direction de port faible n'est donc pas une opération
exacte ; son effet entre dans E et doit rester dans la borne.

## 4. Borne uniforme sans annulations Krylov supposées

Pour une matrice W à n lignes, définir la norme duale de K par

\[
\|W\|_{K^{-1},2}=\|K^{-1/2}W\|_2
=\sqrt{\lambda_{\max}(W^TK^{-1}W)}.
\]

La norme est spectrale : elle contrôle simultanément toutes les
combinaisons linéaires unitaires des ports.

Pour un entier q ≥ 1, définissons

\[
\delta_{\mathrm{sem}}=\|E\|_{K^{-1},2},\qquad
d_j=\|FT^jD\|_{K^{-1},2}\quad(0\le j\le q-2),
\]

\[
\eta_q=\|FT^{q-1}\|_{K^{-1},2}.
\]

La somme qui suit est vide pour q = 1. Posons

\[
\mathcal R_q(\rho)=
\delta_{\mathrm{sem}}
+\sum_{j=0}^{q-2}\rho^{j+1}d_j
+\frac{\rho^q\eta_q\|D\|_2}{1-\rho t}.
\]

**Proposition 2.** Pour toute base de la section 2 et tout q ≥ 1,

\[
\boxed{
\sup_{0\le\mu\le\rho}\|H(\mu)-H_V(\mu)\|_2
\le\frac{\mathcal R_q(\rho)^2}{1-\rho}.}
\]

**Preuve.** L'identité géométrique finie donne, sans hypothèse sur F,

\[
R(\mu)=E+
\sum_{j=0}^{q-2}\mu^{j+1}FT^jD
+\mu^qFT^{q-1}(I_r-\mu T)^{-1}D.
\]

Pour q = 1, c'est directement la première formule de la proposition 1.
Comme T est symétrique positive et μt < 1,

\[
\|(I_r-\mu T)^{-1}\|_2=\frac1{1-\mu t}
\le\frac1{1-\rho t}.
\]

L'inégalité triangulaire et la sous-multiplicativité donnent
‖R(μ)‖_(K⁻¹,2) ≤ ℛ_q(ρ). Enfin, la proposition 1 donne

\[
\|H-H_V\|_2
\le\frac{\|R\|_{K^{-1},2}^2}{1-\mu}
\le\frac{\mathcal R_q(\rho)^2}{1-\rho}.
\]

Cette majoration vaut entre les fréquences échantillonnées, puisqu'elle
est établie pour tout μ de la bande. Elle peut être pessimiste : elle
ignore les compensations entre les termes du résidu et sépare certaines
normes de produits. Augmenter q dans cette seule formule, sans enrichir
V, ne constitue pas une amélioration de l'approximation H_V.

Les quantités peuvent se calculer avec des résolutions dans K et des
Gram de petite taille. Par exemple, si J = FᵀK⁻¹F, alors

\[
d_j^2=\lambda_{\max}\big(D^TT^jJT^jD\big),\qquad
\eta_q^2=\lambda_{\max}\big(T^{q-1}JT^{q-1}\big).
\]

Une évaluation directe des vecteurs FT^jD et de leurs Gram permet aussi
de contrôler les défauts. Il faut éviter de décréter J positif parce que
sa formule théorique l'est, ou de remplacer silencieusement une valeur
propre négative calculée par zéro. La différence algébriquement équivalente
V̄ᵀM̂K⁻¹M̂V̄−T² peut perdre beaucoup de précision par cancellation.

La borne est absolue dans les coordonnées de ports choisies. Si B a été
normalisé par une transformation de colonnes P, elle porte d'abord sur
le transfert PᵀHP. Le retour aux unités d'origine exige de transporter
la congruence et la borne. Une affirmation relative demande en plus de
préciser le dénominateur et les éventuelles directions de port nulles.

## 5. Moments exacts, dimension variable et déflation

Supposons maintenant l'inclusion exacte

\[
\operatorname{range}\overline V\supseteq
\mathcal K_q(K^{-1}\widehat M,K^{-1}B)
=\operatorname{range}
\left[K^{-1}B,\ (K^{-1}\widehat M)K^{-1}B,\ldots,
(K^{-1}\widehat M)^{q-1}K^{-1}B\right].
\]

Il n'est pas nécessaire que la dimension soit qp. Une dépendance exacte
entre colonnes peut réduire la dimension sans perdre cette inclusion.
En revanche, une déflation à seuil peut la perdre ; le nombre de passes
effectuées ne suffit alors pas à revendiquer les moments ci-dessous.

Pour la preuve, introduire

\[
\mathcal A=K^{-1/2}\widehat MK^{-1/2},\qquad
C=K^{-1/2}B,\qquad Q=K^{1/2}\overline V.
\]

Alors QᵀQ = I, T = Qᵀ𝒜Q, D = QᵀC et
range Q contient C, 𝒜C, …, 𝒜^(q−1)C. Par récurrence,

\[
\mathcal A^jC=QT^jD\quad(0\le j\le q-1).
\]

**Proposition 3.** Les 2q premiers moments sont identiques :

\[
\boxed{C^T\mathcal A^jC=D^TT^jD,
\qquad 0\le j\le2q-1.}
\]

**Preuve.** Pour j ≤ 2q−2, écrire j = a+b avec a,b ≤ q−1 et appliquer
l'identité précédente aux deux facteurs de
(𝒜^aC)ᵀ𝒜^bC. Pour j = 2q−1, écrire

\[
C^T\mathcal A^{2q-1}C
=(\mathcal A^{q-1}C)^T\mathcal A(\mathcal A^{q-1}C)
=D^TT^{2q-1}D.
\]

Il en résulte E = 0 et FT^jD = 0 pour 0 ≤ j ≤ q−2. Les développements
des résolvantes autour de μ = 0 donnent

\[
H-H_V=\mu^{2q}D^TT^{q-1}JT^{q-1}D+O(\mu^{2q+1}),
\qquad J=F^TK^{-1}F.
\]

Le coefficient affiché est positif semi-défini ; il peut être nul si
l'espace est déjà invariant dans les directions concernées. À matrices
fixées, cela donne H−H_V = O(ω^(4q)) lorsque ω → 0. Pour q = 1, les
moments statique et inertiel sont reproduits, et l'erreur commence au plus
tôt en ω⁴. Ce résultat local en fréquence ne remplace pas la proposition 2
sur une bande finie, notamment à proximité d'une résonance intérieure.

Dans le cas exact, la proposition 2 se simplifie. Si
j_q = λ_max(T^(q−1)JT^(q−1)), alors

\[
\|H-H_V\|_2\le
\frac{\rho^{2q}j_q\|D\|_2^2}
{(1-\rho t)^2(1-\rho)}.
\]

Mieux, en ordre de Loewner,

\[
0\preceq H(\mu)-H_V(\mu)
\preceq\gamma_qH(0)\preceq\gamma_qH(\mu),\qquad
\gamma_q=\frac{\rho^{2q}j_q}{(1-\rho t)^2(1-\rho)}.
\]

Cette dernière forme utilise E = 0, donc H(0) = DᵀD. Elle n'est pas à
appliquer en effaçant un défaut de semence mesuré.

## 6. Complétion de type Gauss–Radau : résultat théorique

Cette section propose un encadrement plus précis, distinct de la borne
résiduelle précédente. La complétion n'est pas présentée comme une
fonctionnalité implémentée du prototype.

Supposons E = 0 et choisissons une constante établie a telle que
a > λ_max(𝒜) et ρa < 1. Une constante λ_* strictement inférieure au vrai
minimum généralisé permet de prendre a = 1. Dans une complétion orthogonale
de Q, écrire uniquement pour la preuve

\[
\mathcal A=
\begin{pmatrix}T&F_\perp^T\\F_\perp&A_{22}\end{pmatrix},
\qquad C=\binom D0.
\]

Le complément de Schur de aI−𝒜 donne

\[
A_{22}\preceq aI-F_\perp(aI-T)^{-1}F_\perp^T.
\]

Remplacer A₂₂ par cette majorante produit une matrice 𝒜_U telle que
0 ≺ 𝒜 ≼ 𝒜_U ≼ aI. Pour μ ≥ 0 et μa < 1, on en déduit

\[
(I-\mu\mathcal A)^{-1}
\preceq(I-\mu\mathcal A_U)^{-1}.
\]

Or F_⊥ᵀF_⊥ = J. Choisir une factorisation exacte J = L_RᵀL_R, où L_R
a s = rang J lignes, et former la petite matrice

\[
T_U=\begin{pmatrix}
T&L_R^T\\
L_R&aI_s-L_R(aI-T)^{-1}L_R^T
\end{pmatrix},\qquad D_U=\binom D0.
\]

Les directions de l'orthogonal de range F_⊥ sont découplées des ports
dans 𝒜_U. La compression aux directions restantes donne donc

\[
\boxed{
H_V(\mu)\preceq H(\mu)\preceq H_U(\mu),\qquad
H_U=D_U^T(I-\mu T_U)^{-1}D_U.}
\]

Si J = 0, l'espace est invariant et la borne supérieure se confond avec
H_V = H. Pour l'espace Krylov exact de q blocs sans enrichissement
supplémentaire, s ≤ p : l'extension ajoute au plus un bloc de ports.
Pour une base générale, s peut aller jusqu'à r. La preuve par Schur ne
demande pas une récurrence à blocs de taille constante.

L'extension reproduit le moment 2q également. En effet, elle représente
exactement l'action de 𝒜 sur range Q ; les vecteurs 𝒜^qC sont donc
représentés avec leur composante résiduelle. Le produit
(𝒜^qC)ᵀ𝒜^qC donne l'égalité du moment 2q. Ainsi
H_U−H = O(μ^(2q+1)) sous l'inclusion Krylov exacte.

H_V et H_U sont chacun croissants en μ, puisque leurs dérivées ont la
forme Dᵀ(I−μT)⁻¹T(I−μT)⁻¹D avec T positive. Il n'est pas nécessaire
de supposer leur différence croissante. Sur un intervalle [μ_l,μ_r],

\[
0\preceq H(\mu)-H_V(\mu)
\preceq H_U(\mu_r)-H_V(\mu_l).
\]

Une subdivision fournit donc un encadrement uniforme avec de petites
matrices. Échantillonner seulement le gap aux extrémités, sans cet
argument, ne fournit pas à lui seul la même garantie.

Si E ≠ 0, C a aussi une composante dans l'orthogonal de Q : la formule
utilisant D_U = [D;0] ne s'applique plus. Il faudrait notamment enrichir
par les directions manquantes de K⁻¹B, ou dériver un autre encadrement.
Le diagnostic résiduel de la proposition 2 couvre déjà ce cas.

## 7. Constantes établies et calcul flottant

La constante λ_* est une hypothèse substantielle. Un minimum de Ritz de
la paire (K,M), obtenu dans un sous-espace, est en général supérieur ou
égal au vrai minimum. De même, un maximum de Ritz de 𝒜 est inférieur ou
égal au vrai maximum. Substituer ces estimations aux bornes requises peut
placer une résonance intérieure dans la bande prétendument sûre.

Une inégalité de coercivité issue du modèle, ou une vérification rigoureuse
de K−λ_*M ≻ 0 avec contrôle des erreurs, peut établir la constante. Le
succès d'une factorisation Cholesky ordinaire et un petit résidu de vecteur
propre restent des diagnostics : ils ne prouvent pas, à eux seuls, la
position de tout le spectre en arithmétique exacte.

La proposition 2 garde explicitement les défauts E et FT^jD. Elle n'exige
donc aucune annulation Krylov supposée en machine. Elle ne borne cependant
pas automatiquement les erreurs produites en évaluant ces défauts, les
résolutions dans K, les puissances de T ou les valeurs propres des Gram.
Pour une certification, il faudrait aussi majorer ces erreurs et utiliser
des marges dirigées dans 1−ρt et 1−ρ. Une tolérance numérique ou un accord
avec une référence dense ne remplace pas ces majorations.

Une identité utile lorsque la résolution réduite est inexacte vaut pour
toute matrice candidate X̂, avec R = B−A_μX̂ :

\[
\widetilde H=\widehat X^TB+B^T\widehat X
-\widehat X^TA_\mu\widehat X,
\qquad
H-\widetilde H=R^TA_\mu^{-1}R\succeq0.
\]

Elle remplace le transfert brut BᵀX̂ par sa correction énergétique
symétrique. Son évaluation flottante demande encore un contrôle de ses
propres erreurs. Une vérification de résidu est nécessaire, mais un
résidu calculé très petit peut lui-même être affecté par cancellation.

**Proposition 4 : conserver le vrai Gram dans le transfert corrigé.**
Considérons une matrice V̄ quelconque, sans supposer V̄ᵀKV̄ = I. Définir
avec les matrices physiques

\[
G_b=\overline V^TK\overline V,\qquad
M_b=\overline V^T\widehat M\overline V,\qquad
T=M_b,\qquad D=\overline V^TB.
\]

Supposer ρt < 1, avec t = ‖T‖₂, et poser Y = (I−μT)⁻¹D. Si T est
positive semi-définie, t = λ_max(T). Il ne faut plus déduire t ≤ 1 du
seul cadre physique lorsque le Gram G_b n'est pas I. Utiliser le transfert
énergétique corrigé

\[
\boxed{
\widetilde H_V
=D^TY+Y^TD-Y^T(G_b-\mu M_b)Y
=2\,\operatorname{sym}(D^TY)-Y^T(G_b-\mu M_b)Y,}
\]

où sym(Z) = (Z+Zᵀ)/2. Avec E = B−KV̄D et F = M̂V̄−KV̄T, on a encore

\[
R=E+\mu F(I-\mu T)^{-1}D,
\qquad
H-\widetilde H_V=R^TA_\mu^{-1}R\succeq0.
\]

**Preuve.** L'identité du résidu utilise seulement Y−D = μTY, sans
intervention de G_b. La formule énergétique précédente appliquée à
X̂ = V̄Y donne exactement le transfert encadré. Enfin
‖(I−μT)⁻¹‖₂ ≤ 1/(1−ρt) permet de reprendre l'expansion finie et toute
la preuve de la proposition 2. Avec les défauts calculés pour cette V̄,

\[
\boxed{
\sup_{0\le\mu\le\rho}\|H-\widetilde H_V\|_2
\le\frac{\mathcal R_q(\rho)^2}{1-\rho},}
\]

où le dénominateur 1−ρ provient de la coercivité physique exacte
K ≽ M̂, et le dénominateur 1−ρt contrôle le système réduit choisi.
La preuve reste valable pour des colonnes dépendantes, dès lors que le
système réduit ainsi défini est inversible sur la bande.

Cette variante traite algébriquement un défaut d'orthogonalité ; elle ne
certifie toujours pas les arrondis qui interviennent dans les produits,
résolutions et évaluations de la borne. Si G_b = I, elle se réduit au
transfert Galerkin H_V. Pour une base ou un candidat arbitrairement
mauvais, H̃_V peut en revanche être indéfinie : elle est un minorant
énergétique, sans garantie automatique de positivité, de passivité ou
de conservation des moments. Les propositions 3 et la complétion de la
section 6 conservent leurs hypothèses exactes.

## 8. Passage au Schur et limites pour la réponse mécanique

Si le Schur exact de l'interface est S(z) = A_SS(z)−H(z), alors

\[
A_{SS}-H_U\preceq S\preceq A_{SS}-H_V.
\]

La réduction Galerkin sous-estime le transfert et surestime donc ce Schur.
Les signes ne doivent pas être inversés lors de l'assemblage. Des
contributions d'intérieurs indépendants s'additionnent ; les bornes
matricielles s'additionnent également, après transport aux mêmes ports.

Dans une partition physique, le couplage peut dépendre de la fréquence :
B(z) = B₀−zB₁ lorsque la masse possède des termes entre intérieur et
interface. Le traitement à B constant ne couvre pas silencieusement ce
cas. On peut construire les ports élargis B_aug = [B₀,B₁], puis utiliser
B(z) = B_aug P(z), avec P(z) = [I;−zI]. Les transferts et leurs
encadrements sont transportés par la congruence P(z)ᵀ(·)P(z). Pour une
borne uniforme en norme, la variation de P(z) doit aussi être incluse.

Un petit défaut de transfert ne garantit pas un petit défaut de
déplacement près d'une résonance du mécanisme assemblé. Par exemple, si
S est inversible, Ŝ = S+ΔS et Ŝû = f, alors

\[
u-\widehat u=S^{-1}\Delta S\widehat u,
\qquad
\|u-\widehat u\|\le\|S^{-1}\|\,\|\Delta S\|\,\|\widehat u\|.
\]

Le facteur de stabilité global ne disparaît pas parce que chaque
intérieur est sous sa première résonance. Une certification de réponse,
de mode propre, de trajectoire ou de gradient demande des arguments
supplémentaires adaptés à ces objets.

Enfin, le cadre présent ne traite pas directement une raideur intérieure
indéfinie, les mouvements rigides éliminés sans transport, les changements
de configuration, les impacts, l'amortissement général ou la gyroscopie.
Une base locale n×r peut déjà coûter cher si les ports ou r deviennent
nombreux. L'absence de base dense globale ne démontre ni complexité
linéaire générale, ni gain universel sur un solveur multicorps.

## 9. Même contrôle uniforme pour une référence HCB

Dans les mêmes coordonnées de ports, soit Ψ une extension statique
calculée, éventuellement inexacte, et Q une base de r modes d'intérieur
à interface fixe. Les p colonnes de Ψ reconstruisent l'intérieur depuis
les p coordonnées de port. On conserve tous les blocs physiques projetés,
y compris ceux qui résultent du résidu statique R₀ = KΨ+B.

Définir, avec le vrai Gram de Q,

\[
G=Q^TKQ,\quad T=Q^T\widehat MQ,\quad
C=\widehat M\Psi,\quad D=Q^TC,\quad J=Q^TR_0.
\]

L'élimination des coordonnées modales donne le champ intérieur candidat

\[
U(\mu)=\Psi+QY(\mu),\qquad
Y(\mu)=(G-\mu T)^{-1}(\mu D-J).
\]

Posons

\[
E_0=R_0-KQG^{-1}J,\qquad
E=C-KQG^{-1}D,\qquad
F=\widehat MQ-KQG^{-1}T.
\]

**Proposition 5.** Le résidu intérieur positif R = A_μU+B satisfait
exactement

\[
\boxed{R=E_0-\mu E-\mu FY.}
\]

**Preuve.** De GY−μTY = μD−J, on tire
Y = μG⁻¹D−G⁻¹J+μG⁻¹TY. Substituer cette expression au terme KQY dans
R = R₀−μC+(KQ−μM̂Q)Y donne l'identité. Elle ne suppose ni G = I,
ni résidu statique nul, ni diagonalité exacte des blocs modaux.

Avec t = ‖T‖₂, supposons la marge suffisante
η = λ_min(G)−ρt > 0. Définissons
δ₀ = ‖E₀‖_(K⁻¹,2), δ_E = ‖E‖_(K⁻¹,2) et δ_F = ‖F‖_(K⁻¹,2).
L'inégalité ‖Y‖₂ ≤ (ρ‖D‖₂+‖J‖₂)/η donne alors

\[
\mathcal R_{\mathrm{HCB}}=
\delta_0+\rho\delta_E
+\frac{\rho\delta_F(\rho\|D\|_2+\|J\|_2)}{\eta}.
\]

Le Schur énergétique du champ candidat est

\[
S_{\mathrm{HCB}}=A_{SS}+B^TU+U^TB+U^TA_\mu U.
\]

En complétant le carré autour de U_* = −A_μ⁻¹B, on obtient

\[
S_{\mathrm{HCB}}-S_{\mathrm{exact}}
=R^TA_\mu^{-1}R\succeq0,\qquad
\boxed{
\sup_{0\le\mu\le\rho}\|S_{\mathrm{HCB}}-S_{\mathrm{exact}}\|_2
\le\frac{\mathcal R_{\mathrm{HCB}}^2}{1-\rho}.}
\]

Pour r = 0, les termes modaux sont vides : R = R₀−μC et la borne
utilise δ₀+ρ‖C‖_(K⁻¹,2). L'argument contrôle donc aussi l'extension
statique seule. Toutes ces bornes ont le même statut conditionnel que
celles de Krylov : leur évaluation flottante ne certifie pas les arrondis.

Les deux familles peuvent ainsi être sélectionnées avec une tolérance
absolue commune sur le Schur, dans la même métrique et sur la même bande.
Cette comparaison reste distincte d'une tolérance de déplacement et ne
suppose pas que leurs majorants aient la même finesse. La référence HCB
doit conserver le vrai Schur statique et les couplages projetés ; remplacer
un de ces blocs par une identité de normalisation changerait le problème.

Le contrôle HCB réutilise la LU de K déjà employée pour l'extension
statique et les modes. Une résolution groupée du bloc [E₀,E,F], de
2p+r colonnes, suffit aux trois Gram résiduels ; s'y ajoutent les petits
calculs avec G et T. Ce coût appartient à la préparation. La dimension
totale du modèle HCB est **n_modes+p**, et celle du modèle Krylov est
**r_Krylov+p**. Le nombre de blocs Krylov n'est pas une dimension à
comparer directement au nombre de modes HCB. Ces décomptes supposent que
les p coordonnées d'interface sont effectivement conservées.

## 10. Pourquoi beaucoup de modes HCB peuvent être nécessaires : la chaîne

Cette section donne une explication spectrale exacte pour la chaîne
discrète de `ci/modeles_ports.py`. Elle concerne les modes à interface
fixe retenus par ordre de fréquence, et non toute méthode possible de
réduction modale. Les nombres ci-dessous proviennent des formules, sans
campagne ni comparaison chronométrée.

### 10.1 Erreur discrète exacte dans la métrique du port

Considérons n masses m, n ressorts k, u₀ = 0 et u_n comme unique port.
Les n−1 inconnues intérieures ont une masse mI et la raideur tridiagonale
k·tridiag(−1,2,−1). La raideur statique vue du port vaut k/n ; la
normalisation est donc W = √(n/k). Les modes intérieurs normalisés en
masse, leurs valeurs propres et leurs couplages au port normalisé sont

\[
\phi_j(i)=\sqrt{\frac{2}{mn}}\sin\frac{ij\pi}{n},\qquad
\lambda_j=\frac{4k}{m}\sin^2\frac{j\pi}{2n},\qquad
b_j^2=\frac{2k}{m}\sin^2\frac{j\pi}{n},
\]

pour 1 ≤ i,j ≤ n−1. En effet, B = K_IS W = −√(nk)e_(n−1) et
b_j = φ_jᵀB ; l'orthogonalité des sinus donne φ_jᵀMφ_l = δ_jl.

L'extension statique est u_i = (i/n)u_n. Sa masse projetée, y compris
la masse du port, vaut dans ces coordonnées

\[
d_n=mW^2\sum_{i=1}^n(i/n)^2
=\frac{m}{k}\frac{(n+1)(2n+1)}6.
\]

Le transfert intérieur exact est Σ_j b_j²/(λ_j−z). Développer chaque
terme jusqu'à l'ordre un avec un reste exact donne

\[
\frac1{\lambda_j-z}
=\frac1{\lambda_j}+\frac z{\lambda_j^2}
+\frac{z^2}{\lambda_j^2(\lambda_j-z)}.
\]

La raideur statique normalisée vaut 1, et le coefficient de z est d_n.
L'élimination HCB conserve les r premiers termes de correction dynamique.
Ainsi, pour 0 ≤ r ≤ n−1 et 0 ≤ z < λ₁,

\[
S_{n,r}^{\mathrm{HCB}}(z)
=1-zd_n-\sum_{j=1}^{r}
\frac{z^2b_j^2}{\lambda_j^2(\lambda_j-z)},
\]

\[
\boxed{
E_{n,r}(z)=S_{n,r}^{\mathrm{HCB}}(z)-S_n(z)
=\sum_{j=r+1}^{n-1}
\frac{z^2b_j^2}{\lambda_j^2(\lambda_j-z)}\ge0.}
\]

Cette formule est aussi une conséquence directe du couplage massique
entre l'extension statique et le mode j, égal à −b_j/λ_j. Son élimination
soustrait exactement le terme dynamique affiché. Pour r = n−1, la somme
restante est vide et HCB est exact en arithmétique exacte.

### 10.2 Limite de la somme, uniforme sous la première résonance intérieure

Écrivons z = a²λ₁ avec 0 ≤ a ≤ a_max < 1, et gardons r fixé pendant
que n tend vers l'infini. Il s'agit d'une limite après normalisation de
la fréquence et du port. En posant x = π/(2n), le terme j de l'erreur est

\[
e_{n,j}(a)
=\frac{2a^4\sin^4x\cos^2(jx)}
{\sin^2(jx)\,[\sin^2(jx)-a^2\sin^2x]}
=\frac{2a^4(\sin x/\sin jx)^4\cos^2(jx)}
{1-a^2(\sin x/\sin jx)^2}.
\]

Pour chaque j fixé, cette expression converge uniformément en a vers
2a⁴/[j²(j²−a²)]. Pour justifier le passage à la somme, prolonger
e_(n,j) par zéro lorsque j ≥ n. Pour 1 ≤ j < n, les inégalités
sin(jx) ≥ 2jx/π, sin x ≤ x et sin(jx) ≥ sin x donnent

\[
0\le e_{n,j}(a)
\le\frac{\pi^4a_{\max}^4}{8(1-a_{\max}^2)}\frac1{j^4}.
\]

Le majorant est sommable et indépendant de n et de a. On peut donc
couper la somme à un indice fixe : la partie finie converge uniformément,
et les deux queues sont arbitrairement petites grâce à Σj⁻⁴. Cela prouve,
sans interversion non justifiée des limites,

\[
\boxed{
E_{n,r}(a^2\lambda_1)\longrightarrow
E_r(a)=2a^4\sum_{j=r+1}^{\infty}\frac1{j^2(j^2-a^2)},}
\]

uniformément pour 0 ≤ a ≤ a_max. De plus,
a²λ₁d_n → π²a²/3 uniformément, si bien que

\[
S_r^{\mathrm{HCB}}(a)
=1-\frac{\pi^2a^2}{3}
-2a^4\sum_{j=1}^{r}\frac1{j^2(j^2-a^2)}.
\]

La limite exacte du Schur peut se retrouver indépendamment par la
récurrence de la chaîne. Posons θ = 2arcsin(a sin(π/(2n))). La solution
à déplacement de port imposé vérifie u_i/u_n = sin(iθ)/sin(nθ), et

\[
S_n(a^2\lambda_1)
=n[\cos\theta-1+\sin\theta\cot(n\theta)].
\]

Ici nθ → πa uniformément et 0 ≤ nθ ≤ πa_max. La fonction
f(t) = t cot t, prolongée par f(0) = 1, est continue sur cet intervalle.
Écrire n sin θ cot(nθ) = (sin θ/θ)f(nθ), et utiliser
|n(cos θ−1)| ≤ π²a_max²/(2n), prouve la convergence uniforme

\[
\boxed{S(a)=\pi a\cot(\pi a),\qquad S(0)=1.}
\]

On obtient donc aussi S_r^HCB−S = E_r par deux descriptions cohérentes
de la même limite. Les zéros de S, notamment a = 1/2, correspondent à
des résonances globales : la formule contrôle toujours l'erreur absolue
du Schur, sans garantir une erreur relative de déplacement à ces points.

### 10.3 Décroissance algébrique et dimension requise

Pour r ≥ 1 et 0 ≤ a < 1,

\[
\frac1{j^4}\le\frac1{j^2(j^2-a^2)}
\le\frac1{1-a^2/(r+1)^2}\frac1{j^4},\qquad j\ge r+1.
\]

Les comparaisons avec les intégrales de x⁻⁴ sur [r+1,∞[ et [r,∞[
donnent les bornes explicites

\[
\boxed{
\frac{2a^4}{3(r+1)^3}
\le E_r(a)
\le\frac{2a^4}{3r^3[1-a^2/(r+1)^2]}.}
\]

Ainsi E_r(a) ~ 2a⁴/(3r³) pour a fixé strictement positif. Le rapport
asymptotique est également uniforme pour 0 < a ≤ a_max < 1 ; à a = 0,
l'erreur est exactement nulle. Chaque terme de la série augmente avec a,
donc l'erreur maximale sur la bande est E_r(a_max).

Pour une tolérance absolue ε sur cette bande, avec a_max > 0 fixé,
le nombre minimal de modes dans la limite continue satisfait

\[
r_\varepsilon\sim
\left(\frac{2a_{\max}^4}{3\varepsilon}\right)^{1/3}
\qquad(\varepsilon\to0).
\]

La minoration impose déjà
r ≥ (2a_max⁴/(3ε))^(1/3)−1. À titre de lecture des formules,
a = 0,8 et r = 64 donnent un terme asymptotique d'environ 1,0417·10⁻⁶,
et la série limite environ 1,0176·10⁻⁶. Un nombre élevé de modes HCB peut
donc provenir d'une queue spectrale algébrique, même quand la fréquence
reste sous le premier mode intérieur.

L'ordre des limites est essentiel : les preuves ci-dessus prennent
d'abord n → ∞ à r fixé, puis r → ∞. À n fixé, la dimension n−1 suffit
pour l'exactitude ; une loi ε^(−1/3) ne peut pas s'étendre indéfiniment.
Un régime où r et n croissent ensemble demanderait un contrôle conjoint
supplémentaire. Ces résultats ne prouvent pas un taux global pour Krylov.
La conservation des moments O(ω^(4q)) à n fixé ne fournit, à elle seule,
ni une borne uniforme en n, ni une convergence exponentielle en q sur
toute la bande. Les constantes et les défauts de la section 4 doivent
encore être maîtrisés pour établir une telle comparaison.

## 11. Rattachement scientifique récent

Zimmerling, Druskin et Simoncini établissent en 2025 des encadrements
matriciels par Gauss et Gauss–Radau pour (A+sI)⁻¹ avec A définie positive
et s > 0. Voir leur théorème 1, proposition 3 et formule (23). Leur
analyse explique le lien entre Lanczos par blocs, quadratures et fractions
continues matricielles de Stieltjes ; l'extension Radau à q+1 blocs utilise
les informations de q étapes de Lanczos.

- [Article publié, Journal of Scientific Computing, 2025](https://link.springer.com/article/10.1007/s10915-025-02799-z).
- [Version intégrale des auteurs, janvier 2025](https://arxiv.org/pdf/2407.21505).

Notre résolvante sous résonance se ramène à ce domaine par
Ã = aI−𝒜 et s = 1/μ−a > 0 :
(I−μ𝒜)⁻¹ = μ⁻¹(Ã+sI)⁻¹ pour μ > 0. La valeur μ = 0 s'obtient par
continuité. La complétion de la section 6 est démontrée directement par
ordre de Schur, sans étendre une borne de résolvante positive à un
domaine indéfini. Aucun avantage d'une moyenne de quadratures n'est
supposé pour tous les mécanismes.
