# Quotients, stabilité et dimension : preuves pour le noyau généraliste

Carnet mathématique du 7 septembre 2026. Trois analyses indépendantes,
puis une confrontation des hypothèses, complètent ici la
[recherche fondamentale](FONDEMENTS_MATHEMATIQUES_VINKULUM_2026.pdf) et le
[sondage du quotient orthogonal](CONTRAINTES_ORTHOGONALES_PROTOTYPE.md).
Les propositions sont redérivées, sans revendication de nouveauté.
Elles définissent des conditions de conception ; elles ne sont pas toutes
implémentées dans le prototype.

La cible scientifique devient précise : **conserver exactement les
compatibilités et la métrique, maîtriser les transferts aux interfaces,
puis réduire les solutions utiles sous contrôle d'erreur**. Une faible
dimension d'interface ne garantit ni stabilité ni faible dimension dynamique.

## 1. Ce qui doit traverser une interface : un quotient exact

Soit G:V→C le jacobien de contraintes homogènes. Après choix des produits
scalaires, les mobilités infinitésimales forment ker G ; les multiplicateurs
qui ne produisent aucune force forment ker G*. Le théorème du rang donne :

\[
\dim\ker G-\dim\ker G^*=\dim V-\dim C.
\]

Une perte de rang crée simultanément des mobilités infinitésimales et des
auto-contraintes. Compter seulement les équations ne distingue pas ces effets.
Les constructions des sections 1–2 sont à configuration fixée, pour les
vitesses ou le problème linéarisé. Elles ne fournissent pas une
paramétrisation finie globale des contraintes non linéaires. Si une
injection E(q) sert à intégrer la dynamique, sa dérivée temporelle
contribue aussi à l'accélération et aux termes inertiels.

Partitionnons V=V_I⊕V_S et G=[A B]. L'intérieur x peut prolonger une valeur
s de l'interface si et seulement si Ax=−Bs est soluble. Soit
π:C→C/im A la projection canonique. L'opérateur d'interface exact est :

\[
R=\pi B:V_S\longrightarrow\operatorname{coker}A.
\]

**Proposition.** Il existe la suite exacte :

\[
0\to\ker A\to\ker G\xrightarrow{\mathrm{trace}}V_S
\xrightarrow{R}\operatorname{coker}A
\to\operatorname{coker}G\to0.
\]

**Preuve.** Une mobilité de trace nulle est exactement une mobilité
intérieure de ker A. L'image de la trace est ker R, puisque Bs doit
appartenir à im A. L'image de R est la partie de C/im A engendrée par
B ; quotienter par cette image donne C/im G. Cela vérifie toutes les
exactitudes. En particulier :

\[
\operatorname{rang}G=\operatorname{rang}A+\operatorname{rang}R,
\qquad
\dim\ker G=\dim\ker A+\dim\ker R.
\]

En identifiant coker A à ker A*, R peut se représenter par P₀B, où P₀
projette sur ker A*. Son adjoint a pour domaine ker A*, et
ker R*=ker G*. Les auto-contraintes sont donc conservées, avec les
mobilités internes et les compatibilités de l'interface.

Pour plusieurs intérieurs indépendants, A est diagonale par blocs : ces
quotients se calculent localement, puis se composent. **Les mouvements
internes de ker A restent présents.** Les éliminer en plus demanderait
une approximation dynamique distincte.

## 2. Transport exact de la masse et origine d'un mauvais pivot

Supposons M=diag(M_I,M_S) définie positive et M_I=L_I L_Iᵀ. Posons
F=A L_I^(−T). Une décomposition orthogonale complète de rang r donne :

\[
U^TFV=\begin{pmatrix}T&0\\0&0\end{pmatrix},
\quad U=[U_1,U_0],\quad V=[V_1,V_0],\quad T\text{ inversible}.
\]

Notons B₁=U₁ᵀB et C₀=U₀ᵀB. Les contraintes équivalent exactement à :

\[
x_I=L_I^{-T}V_0z-L_I^{-T}V_1T^{-1}B_1s,
\qquad C_0s=0.
\]

**Preuve.** Écrire x_I=L_I^(−T)V(y,z) et multiplier Ax_I+Bs=0 par Uᵀ.
On obtient Ty+B₁s=0 et C₀s=0. La première équation détermine y ; les
coordonnées z sont libres.

Les colonnes V₀ et V₁ étant orthogonales, les termes croisés de masse
disparaissent. Par conséquent :

\[
x^TMx=z^Tz+s^TM_{\mathrm{eff}}s,
\qquad
M_{\mathrm{eff}}=M_S+B_1^TT^{-T}T^{-1}B_1\succ0.
\]

Les contributions de plusieurs intérieurs s'additionnent. Les injections
peuvent être appliquées par facteurs locaux sans stocker une base globale.
La métrique a une signification physique sur l'espace compatible C₀s=0 ;
son prolongement aux valeurs incompatibles dépend de la reconstruction.

**Contre-exemple à la stabilité déduite du seul rang global.** Pour
G_ε=[ε,1,0], ε≠0 et M=I, la seule valeur singulière vaut √(1+ε²), toujours
séparée de zéro. Pourtant, éliminer la première coordonnée donne
x_I=−s₁/ε et M_eff=diag(1+ε⁻²,1). Le mauvais conditionnement vient
du choix de coordonnées.

Avec la racine carrée symétrique de M_S, posons :

\[
W=T^{-1}B_1M_S^{-1/2},\qquad\theta=\|W\|_2.
\]

Alors M_S^(−1/2) M_eff M_S^(−1/2)=I+WᵀW. Toutes ses valeurs propres
sont dans [1,1+θ²]. C'est une mesure du transfert mécanique, qui tient
compte à la fois du pivot, des couplages et de la masse.

**Construction exacte avec pivots retardés.** Dans une SVD
de F, on n'élimine que les directions σ_j≥τ, avec τ>0. Les autres coordonnées
rejoignent l'interface avec leurs équations σ_j y_j+b_j s=0 intactes.
La correction de masse des directions éliminées satisfait :

\[
0\preceq\Delta M_S
=B_a^T\Sigma_a^{-2}B_a
\preceq\frac{\|B_a\|_2^2}{\tau^2}I.
\]

Ce report est exact ; il élargit potentiellement l'interface. Une valeur
numériquement ambiguë peut ainsi rester une inconnue au lieu d'être déclarée
nulle. Cela diffère du seuil de suppression du prototype actuel.

Après composition, la masse peut comporter un bloc hors diagonale J :
M=[H J;Jᵀ K]. Le changement h=x_I+H⁻¹Js donne la masse
diag(H,K−JᵀH⁻¹J) et remplace B par B−AH⁻¹J. La construction précédente
s'applique de nouveau. **Elle reste locale seulement si J touche
l'interface retenue.** Il faut donc organiser les éliminations avec les
couplages mécaniques courants, au-delà du seul graphe de G.

Pour une factorisation par fronts denses de tailles f_v, comprenant
variables et lignes de contraintes du front, une majoration élémentaire
est O(Σ_v f_v³) en opérations arithmétiques à précision fixée et
O(Σ_v f_v²) coefficients pour leurs facteurs.
La conclusion O(Nw³) exige effectivement O(N) fronts de taille au plus w,
après les reports de pivots et la propagation des couplages. Ni cette
hypothèse ni une meilleure borne ne sont établies pour tous les mécanismes.

## 3. Une singularité du jacobien peut cacher une obstruction d'ordre deux

Considérons trois points du plan reliés par des barres, avec
c_ij(p)=½(‖p_i−p_j‖²−ℓ_ij²). Pour un triangle non collinéaire, les trois
lignes du jacobien sont indépendantes : une auto-contrainte exigerait,
à chaque sommet, l'équilibre de deux directions indépendantes ; tous ses
coefficients sont donc nuls. Le rang vaut trois et la nullité trois.

Aux positions p₁=(0,0), p₂=(1,0), p₃=(2,0), les lignes vérifient :

\[
G_{13}=2(G_{12}+G_{23}),\qquad
\lambda=(-2,-2,1)\in\ker G^T.
\]

Le rang vaut deux et la nullité quatre. Le graphe n'a pas changé.
Mais la mobilité infinitésimale supplémentaire n'est pas intégrable en
mouvement fini. Tout chemin contraint de classe C² doit satisfaire :

\[
Gv=0,\qquad Ga+D^2c[v,v]=0,
\qquad\lambda^TD^2c[v,v]=0.
\]

Or, pour cette auto-contrainte :

\[
\lambda^TD^2c[v,v]
=-2\|v_1-v_2\|^2-2\|v_2-v_3\|^2+\|v_1-v_3\|^2
=-\|v_1-2v_2+v_3\|^2.
\]

La vitesse v₁=v₃=0, v₂=(0,1) est dans ker G mais donne −4. Elle ne
peut être la vitesse d'un chemin contraint C². **Une dépendance linéaire
à une pose ne permet donc pas de retirer la contrainte non linéaire.**
Les contraintes originales et leurs dérivées doivent rester accessibles.

La géométrie gouverne également le conditionnement. En remplaçant p₃
par (2,ε), le rang redevient trois mais :

\[
G(\varepsilon)^T\lambda
=\varepsilon(0,-1,0,2,0,-1)^T,
\qquad
\sigma_{\min}(G(\varepsilon))
\leq\sqrt{2/3}\,|\varepsilon|.
\]

Le graphe garde les mêmes séparateurs tandis que la stabilité se dégrade
arbitrairement. L'ordre combinatoire et le contrôle géométrique sont deux
conditions distinctes.

## 4. Du résidu KKT à une erreur dans la métrique physique

Considérons min ½xᵀAx−fᵀx sous Gx=0, A symétrique. Choisissons S tel
que SᵀMS=I, avec M définie positive et S carrée inversible, puis posons
x=Su, A₀=SᵀAS, h=Sᵀf et B=DGS, D diagonale
inversible. Alors ‖x‖_M=‖u‖₂. Soient K=ker B, P son projecteur orthogonal
et H=(PA₀P)|_K. Supposons, pour 0<rang B<n :

\[
H\succeq\alpha I_K,\quad\alpha>0,\qquad
0<\beta\leq\sigma_{\min}^{+}(B),\qquad
\kappa\geq\|PA_0(I-P)\|_2.
\]

A₀ peut être indéfinie hors de K. La solution u est unique ; ses
multiplicateurs sont définis modulo ker Bᵀ.

Pour un candidat û et des multiplicateurs η̂ quelconques, définir :

\[
r=h-A_0\widehat u-B^T\widehat\eta,\quad g=B\widehat u,
\quad c=(I-P)\widehat u=B^\dagger g,
\quad e=u-P\widehat u.
\]

**Proposition a posteriori.** Avec s=Pr+PA₀c, on a exactement :

\[
He=s,\qquad u-\widehat u=e-c,\qquad e\perp c,
\qquad e^TA_0e=s^TH^{-1}s.
\]

Il en résulte la borne physique :

\[
\boxed{
\|x-\widehat x\|_M^2
\leq
\frac{(\|Pr\|_2+\kappa\|g\|_2/\beta)^2}{\alpha^2}
+\frac{\|g\|_2^2}{\beta^2}.}
\]

**Preuve.** Projeter l'équilibre exact donne PA₀u=Ph. Comme Pû=û−c,
PA₀(u−Pû)=Pr+PA₀c. La coercivité de H et
‖c‖≤‖g‖/β donnent la majoration, puis l'orthogonalité e⊥c donne la somme
des carrés. Pour le candidat corrigé admissible SPû, l'écart d'énergie
vaut exactement ½eᵀA₀e. Cette interprétation énergétique ne s'applique
pas à un candidat inadmissible lorsque A₀ est indéfinie hors du noyau.

Les constantes sont nécessaires. Avec B=(0,1),
A₀=[α κ;κ 0], h=0, û=(−κt/α,t) et η̂=κ²t/α, on a r=0 mais
‖u−û‖²=t²+(κt/α)². La borne est atteinte. L'équilibre seul ne suffit
pas à contrôler le déplacement.

Si K={0}, u=0 et ‖û‖≤‖Bû‖/β, sans constante α. Si B=0, on retrouve
le cas défini positif usuel, sans β. Si des mouvements rigides donnent
α=0, il faut fixer une jauge, vérifier la compatibilité du chargement et
travailler dans le quotient par ces mouvements.

**Précaution sur le noyau approché.** Si P̃ est un projecteur orthogonal
de même rang que P et
‖BP̃‖≤δ, alors θ=‖P−P̃‖≤δ/β. Cela ne préserve pas automatiquement la
coercivité. Définir la compression symétrique
A_⊥=((I−P)A₀(I−P))|_{K⊥}, puis ν≥max(0,−λ_min(A_⊥)). Une condition
suffisante est :

\[
\alpha(1-\theta^2)-2\kappa\theta-\nu\theta^2>0.
\]

Elle résulte de la décomposition d'un vecteur de range P̃ en ses parties
dans K et K⊥. Le contre-exemple A₀=diag(1,−N), B=(0,1),
B̃=(−ε,1) donne sur ker B̃ une raideur (1−Nε²)/(1+ε²) : une
perturbation aussi petite que souhaité peut détruire la coercivité quand
N est grand. Un contrôle du noyau original ne suffit pas.

Ces bornes deviennent des certificats seulement avec des constantes et
des résidus correctement majorés. Des estimations empiriques de α, β ou κ
ne constituent pas des bornes garanties en arithmétique flottante.

## 5. Réduire la dimension : viser les solutions utiles

Pour tout espace strictement réduit V_k⊊K, les projecteurs vérifient
‖P−P_{V_k}‖=1 : il existe un vecteur unitaire de K orthogonal à V_k.
Exiger une petite erreur globale de projecteur interdirait donc toute
réduction dimensionnelle réelle.

La solution de Galerkin u_k∈V_k satisfait en revanche :

\[
\|u-u_k\|_H=\min_{v\in V_k}\|u-v\|_H.
\]

**Preuve.** Le résidu admissible s_k=P(h−A₀u_k) est orthogonal à V_k,
donc l'erreur l'est pour le produit scalaire de H. Le théorème de
Pythagore donne l'optimalité.

Supposons qu'un opérateur positif T sur K vérifie
aH⁻¹≼T≼bH⁻¹, avec 0<a≤b. Alors :

\[
\frac{s_k^TTs_k}{b}
\leq\|u-u_k\|_H^2
\leq\frac{s_k^TTs_k}{a}.
\]

En enrichissant V_k par la direction Ts_k et en recalculant Galerkin :

\[
\boxed{\|u-u_{k+1}\|_H^2
\leq(1-a/b)\|u-u_k\|_H^2.}
\]

**Preuve.** La diminution minimale du carré de l'erreur le long de Ts_k
vaut (s_kᵀTs_k)²/(s_kᵀTHTs_k). Or THT≼bT ; cette diminution est donc
au moins s_kᵀTs_k/b, puis au moins a/b fois l'erreur précédente.
Le nouveau Galerkin fait au moins aussi bien.
Si s_k=0, la solution est déjà exacte et aucun enrichissement n'est nécessaire.

Une borne plus précise, avec les mêmes hypothèses, est
‖u−u_{k+1}‖_H²≤((b−a)/(b+a))²‖u−u_k‖_H². Pour la vérifier, poser
C=H^(1/2)TH^(1/2), w=H^(1/2)(u−u_k) ; son spectre implique
C²≼(a+b)C−abI. Le gain relatif de la recherche linéaire vaut
(wᵀCw)²/[(wᵀw)(wᵀC²w)] et est au moins 4ab/(a+b)². En posant
m=(wᵀCw)/(wᵀw), cette dernière inégalité se réduit à
[(a+b)m−2ab]²≥0.

Cela donne un principe d'adaptation à erreur énergétique imposée. Pour
en tirer un coût compétitif, il reste à construire T par calculs locaux
avec un rapport b/a maîtrisé et à borner le coût de chaque enrichissement.
Les constantes indépendantes de la taille ou des contrastes ne sont pas
acquises simplement en choisissant un préconditionneur.

## 6. Une interface scalaire peut exiger beaucoup de modes dynamiques

Prenons un système linéaire sans amortissement, M=diag(M_I,M_S)≻0,
K=Kᵀ et K_II≻0. Les charges intérieures sont nulles. En régime harmonique,
soient φ_j les vecteurs propres orthonormaux de
A=M_I^(−1/2)K_II M_I^(−1/2), avec valeurs λ_j>0 croissantes, et
b_j=K_SI M_I^(−1/2)φ_j. L'élimination exacte des coordonnées intérieures
donne, hors des pôles :

\[
S(\omega)=K_{SS}-\omega^2M_S
-\sum_j\frac{b_jb_j^T}{\lambda_j-\omega^2}.
\]

**Preuve.** Résoudre la première ligne du système harmonique par
(K_II−ω²M_I)⁻¹, substituer dans la seconde, puis employer la décomposition
spectrale de A. L'expression est une raideur dynamique d'interface.

Pour une interface scalaire, p valeurs λ_j distinctes et p couplages
b_j≠0 donnent p pôles distincts, sans annulation possible de leur résidu.
Un modèle du second ordre avec r coordonnées intérieures possède au plus
r tels pôles en la variable ω², à coefficients constants et sans loi
de bord dépendant de la fréquence. S'il reproduit S exactement sur un
intervalle ouvert sans pôle, l'identité des fonctions rationnelles impose
les mêmes pôles ; **r≥p**. Une interface de dimension un ne garantit donc
pas une réduction dynamique exacte de dimension constante.

Pour une approximation contrôlée, conserver les r premiers modes
dynamiquement et remplacer le résolvant des autres par sa valeur statique
donne d'abord :

\[
S_r(\omega)=K_{SS}-\omega^2M_S
-\sum_{j\leq r}\frac{b_jb_j^T}{\lambda_j-\omega^2}
-\sum_{j>r}\frac{b_jb_j^T}{\lambda_j}.
\]

Cette première approximation ne transporte pas encore la masse du
prolongement statique des modes omis. Elle doit être distinguée de la
réduction mécanique par congruence de la section 2.

Si Ω²<λ_{r+1}, définir T_tail=Σ_{j>r} b_jb_jᵀ/λ_j. Pour |ω|≤Ω,
hors des pôles conservés :

\[
\boxed{
0\preceq S_r(\omega)-S(\omega)
\preceq\frac{\Omega^2}{\lambda_{r+1}-\Omega^2}\,T_{\mathrm{tail}}.}
\]

**Preuve.** Chaque coefficient omis a pour différence
ω²/[λ_j(λ_j−ω²)], positive et au plus
Ω²/[λ_j(λ_{r+1}−Ω²)]. Multiplier par b_jb_jᵀ positif et sommer donne
les inégalités de Loewner.

La queue dépend des couplages, et pas seulement des fréquences propres.
Si, à la fréquence considérée, S(ω)≽aI avec a>0, alors pour la même
charge d'interface f, s=S⁻¹f et s_r=S_r⁻¹f vérifient :

\[
\|s-s_r\|\leq\frac{\|S_r-S\|}{a}\,\|s_r\|.
\]

Cela suit de S(s−s_r)=(S_r−S)s_r. À proximité d'une résonance globale,
une borne uniforme a peut disparaître. Les charges intérieures non nulles
produisent aussi un second membre réduit fréquentiel ; les conditions
initiales produisent une réponse libre. Une réduction de S seule ne
contrôle pas ces contributions. L'amortissement, la non-linéarité et le
contact demandent d'autres hypothèses.

**La conservation de la masse améliore l'ordre de l'erreur.** Dans les
coordonnées massiques propres, utiliser le prolongement :

\[
y_j=q_j\quad(j\leq r),\qquad
y_j=-\frac{b_j^Ts}{\lambda_j}\quad(j>r).
\]

La restriction par congruence de la masse ajoute alors exactement
M_tail=Σ_{j>r} b_jb_jᵀ/λ_j² à M_S. Les coordonnées retenues sont
orthogonales aux modes omis, donc il n'y a pas de terme de masse croisé
entre q et ce prolongement. Le Schur du modèle ainsi réduit vaut :

\[
S_r^M(\omega)=S_r(\omega)-\omega^2M_{\mathrm{tail}}.
\]

En soustrayant les deux premiers termes du développement du résolvant,
on obtient exactement :

\[
S_r^M-S=\sum_{j>r}
\frac{\omega^4}{\lambda_j^2(\lambda_j-\omega^2)}b_jb_j^T,
\qquad
\boxed{0\preceq S_r^M-S
\preceq\frac{\Omega^4}{\lambda_{r+1}-\Omega^2}M_{\mathrm{tail}}.}
\]

À spectre et couplages fixés, la borne passe ainsi de l'ordre Ω² à
l'ordre Ω⁴ à basse fréquence. Ce gain vient du transport cohérent de
la métrique, sous les hypothèses linéaires précédentes. Il ne constitue
ni une invention de réduction modale ni un gain de temps déjà mesuré.
Avec des charges intérieures approximées statiquement, l'erreur de
second membre peut rester d'ordre Ω² ; le gain sur l'opérateur ne
devient donc pas automatiquement un gain du même ordre sur toute réponse.

L'espace ainsi construit est celui de Craig–Bampton avec modes
d'interface fixe et prolongements statiques complets : il suffit de
changer les coordonnées retenues par
q=q_CB−Λ_r⁻¹B_rᵀs pour passer d'un ansatz à l'autre.
Ici B_r=[b₁,…,b_r] et Λ_r=diag(λ₁,…,λ_r).
Cette antériorité est concrète : Exudyn dispose déjà d'une chaîne
[Hurty–Craig–Bampton](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/ModelOrderReductionAndComponentModeSynthesis.html).
L'apport recherché pour Vinkulum est le contrôle adaptatif des interfaces
et de l'erreur dans des compositions plus générales.

**Un contrôle limité aux ports.** La queue peut s'écrire sans calculer
tous les modes omis :

\[
T_{\mathrm{tail}}=K_{SI}K_{II}^{-1}K_{IS}
-\sum_{j\leq r}\frac{b_jb_j^T}{\lambda_j}.
\]

Cela demande une résolution statique par colonne d'interface ; ce coût
n'est pas nul et l'exactitude de la soustraction doit être contrôlée.
De même, M_tail s'obtient à partir de
K_SI K_II⁻¹ M_I K_II⁻¹ K_IS, en retranchant les contributions retenues.
Soit D=Ω² T_tail/(λ_{r+1}−Ω²). Vérifier S_r(ω)−D≽aI>0 suffit
à garantir S(ω)≽aI, sans former le Schur dynamique exact.

Plus généralement, même si S_r est indéfinie, la condition
σ_min(S_r)>‖D‖ donne :

\[
\sigma_{\min}(S)\geq\sigma_{\min}(S_r)-\|D\|>0,
\qquad
\|s-s_r\|\leq
\frac{\|D\|}{\sigma_{\min}(S_r)-\|D\|}\,\|s_r\|.
\]

La première inégalité suit de ‖Sv‖≥‖S_rv‖−‖(S_r−S)v‖ pour tout
vecteur unitaire. Cette vérification agit sur la dimension des ports,
avec des bornes validées sur la queue et sa première valeur propre omise.
Un contrôle à quelques fréquences ne garantit pas une bande entière :
il faut des majorations uniformes sur les intervalles considérés.
Les mêmes preuves s'appliquent au modèle conservant la masse, en
remplaçant S_r par S_r^M et D par Ω⁴ M_tail/(λ_{r+1}−Ω²).

## 7. Conséquences concrètes pour l'architecture

Ces preuves séparent les décisions qui restent à implémenter :

1. Transmettre les compatibilités de coker A et conserver les mouvements
   intérieurs libres, avec transport explicite de la masse.
2. Retarder les directions mal conditionnées au lieu de les déclarer
   nulles ; mesurer l'interface élargie et l'amplification du transfert.
3. Conserver les contraintes non linéaires originales et tester les
   compatibilités d'ordre deux lorsque le rang change.
4. Piloter une réduction supplémentaire par erreur de solution, bande
   fréquentielle, charges et couplages d'interface. Une petite erreur
   de matrice ou un petit résidu brut ne suffit pas.

La prochaine expérience utile est une condensation par sous-structures
qui expose conjointement ces quatre informations. L'assemblage d'outils
mathématiques existants ne démontre aucune avance mondiale ; la cible
vérifiable est un coût inférieur à erreur physique commune, avec les cas
singuliers et les régimes défavorables conservés.
