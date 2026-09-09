# Du résidu des ports au champ mécanique reconstruit

Carnet du 7 septembre 2026. Ce document prolonge les
[preuves de relèvement](PORTS_RELEVEMENT_PREUVES.md) et les
[preuves Krylov](PORTS_KRYLOV_PREUVES.md). Il examine les contrats présents
dans `ci/ports_releves.py` et `ci/chaine_sous_structures.py`, puis les
bornes mises en œuvre expérimentalement dans `ci/champ_interieur.py`.
Les enveloppes locales couvrent la bande annoncée ; le contrôle global
et l'arrêt relatif du prototype portent sur les fréquences demandées.

Les propositions sont redérivées en arithmétique exacte, sans revendication
de nouveauté. Une évaluation flottante de leurs termes reste conditionnelle
tant que les constantes et les erreurs de calcul ne sont pas encadrées.
Aucune mesure de performance ni modification de code n'est fournie ici.

## 1. Domaine, champ exact et erreur à port imposé

Pour une sous-structure, le modèle de référence est K = DᵀD, avec M
symétrique. Les coordonnées sont partitionnées en intérieur I et port S.
On suppose

\[
K_{II}\succeq\lambda_*M_{II}\succ0,\qquad
0\le z=\omega^2\le\Omega^2<\lambda_*.
\]

Posons A(z) = K_II−zM_II, α(z) = λ_*−z et
α_* = λ_*−Ω² > 0. Il n'y a pas de charge appliquée aux coordonnées
intérieures. Les éventuelles charges sont appliquées aux coordonnées
conservées par l'assemblage. Une charge intérieure demanderait une
solution particulière et son résidu supplémentaire.

Une matrice inversible W convertit les coordonnées normalisées du port
s en déplacements physiques u_S = Ws. Le champ exact U(z)s vérifie

\[
U_S=W,\qquad
A(z)U_I+(K_{IS}-zM_{IS})W=0.
\]

Soit Û(z) un champ candidat de même trace, Û_S = W. Son résidu dans
le modèle D d'entrée est

\[
R(z)=D_I^TDÛ-z(MÛ)_I,
\qquad D_I=D[:,I].
\]

En posant E = Û−U, on a exactement

\[
\boxed{E_S=0,\qquad E_I=A(z)^{-1}R(z).}
\]

Les matrices E et R portent p colonnes, une par coordonnée normalisée
de port. Pour une valeur s imposée, l'erreur et le résidu sont Es et Rs.
Cette identité ne dépend ni du choix QR ou LU, ni d'une orthogonalité
parfaite de la base candidate. Elle exige un résidu dans le même modèle
physique que le champ exact auquel on compare.

## 2. Bornes massiques, dynamiques et de déformation

Pour une matrice Z de forces intérieures, définir

\[
\|Z\|_{M_{II}^{-1},2}=\|M_{II}^{-1/2}Z\|_2,
\qquad G_R=R^TM_{II}^{-1}R.
\]

**Proposition 1.** Les trois erreurs de Gram satisfont

\[
\boxed{E_I^TM_{II}E_I\preceq\frac{G_R}{\alpha(z)^2},}
\]

\[
\boxed{E_I^TA(z)E_I\preceq\frac{G_R}{\alpha(z)},\qquad
(DE)^TDE\preceq\frac{\lambda_*G_R}{\alpha(z)^2}.}
\]

**Preuve.** Poser H = M_II⁻¹/²K_IIM_II⁻¹/² ≽ λ_*I. Dans les
coordonnées massiques, E_I devient
(H−zI)⁻¹M_II⁻¹/²R. La norme de (H−zI)⁻¹ est au plus 1/α(z),
ce qui donne la première inégalité. La seconde suit de
A(z)⁻¹ ≼ M_II⁻¹/α(z). Enfin K_II = A(z)+zM_II donne

\[
E_I^TK_{II}E_I
\preceq\left[\frac1\alpha+\frac z{\alpha^2}\right]G_R
=\frac{\lambda_*}{\alpha^2}G_R.
\]

Comme E_S = 0, DE = D_IE_I et EᵀME = E_IᵀM_IIE_I. Les termes
M_IS ne sont donc pas oubliés : leur contribution à cette erreur locale
est exactement nulle parce que sa trace est nulle.

Si une enveloppe uniforme établie donne
sup_z ‖R(z)‖_(M_II⁻¹,2) ≤ δ, alors pour toute la bande

\[
\boxed{
\|E(z)s\|_M\le c_M\|s\|_2,\qquad
\|DE(z)s\|_2\le c_D\|s\|_2,}
\]

\[
c_M=\frac\delta{\alpha_*},\qquad
c_D=\frac{\sqrt{\lambda_*}\,\delta}{\alpha_*}.
\]

Le Gram G_R permet une borne directionnelle plus fine que δ²I quand
le port s est connu. Il suffit alors d'utiliser sᵀG_Rs. Cette évaluation
peut traiter une seule combinaison de charges au lieu de tous les ports.

### Lien avec l'audit existant

L'expansion dans D d'entrée du carnet de relèvement produit un terme
`majorant_residu` égal à une majoration b_res = δ²/α_* de l'erreur
énergétique. Les conséquences pour le champ sont donc

\[
\boxed{c_M=\sqrt{b_{\mathrm{res}}/\alpha_*},\qquad
c_D=\sqrt{\lambda_*b_{\mathrm{res}}/\alpha_*}.}
\]

Il faut employer la constante de coercivité du modèle D d'entrée.
`majorant_total` ajoute l'écart du fonctionnel de Schur renvoyé : ce terme
affecte la réponse de l'assemblage, mais n'est pas un résidu supplémentaire
du champ à port déjà imposé. On peut l'utiliser comme majoration plus
grossière, sans confondre ces deux origines d'erreur.

De même, `norme_champ_dynamique` est actuellement une borne des
coefficients réduits Y. Elle n'est ni une norme physique du champ complet,
ni une borne de son erreur ; un changement d'échelle de la base peut
modifier ces coefficients sans changer le champ.

## 3. Le Schur est quadratique en l'erreur ; le champ ne l'est pas

Le Schur énergétique du candidat et celui du champ exact satisfont

\[
\boxed{S_{c,D}-S_D=R^TA(z)^{-1}R=E_I^TA(z)E_I.}
\]

Le résidu intervient au carré dans l'énergie. Il intervient linéairement
dans E_I = A(z)⁻¹R. Une erreur de Schur de l'ordre ε ne fournit donc
qu'une borne de champ d'ordre √ε, multipliée par les constantes de
stabilité appropriées. Elle ne donne pas une erreur de déplacement ε.

Un contre-exemple scalaire le montre sans estimation. Prenons
K_II = M_II = 1, z = 1−α, α > 0, et un couplage de port
b = √(εα). À déplacement de port unitaire, le candidat intérieur nul
donne une erreur de Schur ε et une erreur de champ √(ε/α).
En choisissant α = ε², le Schur converge tandis que l'erreur de champ
diverge. La raideur complète [[1,b],[b,1]] reste positive pour ε petit.
Le phénomène vient de l'approche d'une résonance intérieure, pas d'un
défaut combinatoire du graphe.

L'évaluation de la précision doit donc conserver au moins trois objets :
erreur de Schur, erreur massique du champ et erreur de déformation.
La précision d'un seul de ces objets ne doit pas être renommée comme
précision des deux autres.

## 4. Observables : déplacement, déformation et efforts

Pour un observable linéaire ℓᵀu_I, la proposition 1 donne

\[
|\ell^TE_Is|
\le\|M_{II}^{-1/2}\ell\|_2\,\|E_Is\|_{M_{II}}.
\]

Pour plusieurs observables regroupés dans O_I, remplacer la première
norme par ‖O_IM_II⁻¹/²‖₂. À port imposé, une composante O_Su_S
n'a pas d'erreur puisque les traces coïncident. Après une résolution
sous force, l'erreur de port doit aussi être incluse.

Une observable de déformation O = TD peut utiliser directement
‖OEs‖₂ ≤ ‖T‖₂‖DEs‖₂. Les contraintes et efforts élémentaires
demandent leurs facteurs constitutifs et unités propres. Une norme
d'énergie intégrée ne fournit pas automatiquement une borne uniforme
sur la contrainte maximale ponctuelle. Sur une suite de maillages, la
norme duale d'une observation ponctuelle peut croître : aucune équivalence
de normes indépendante du maillage n'est supposée.

**Correction par adjoint local.** Soit v̂ une approximation de la solution
de A(z)v = ℓ et τ = ℓ−A(z)v̂ son résidu. Alors

\[
\boxed{\ell^TE_Is
=\widehat v^TRs+\tau^TA(z)^{-1}Rs,}
\]

\[
\left|\ell^TE_Is-\widehat v^TRs\right|
\le\frac{\|\tau\|_{M_{II}^{-1}}\,\|Rs\|_{M_{II}^{-1}}}{\alpha(z)}.
\]

La preuve consiste à remplacer ℓ par A(z)v̂+τ. Un adjoint précis
pour une observable peut donc donner une estimation plus fine qu'une
majoration du champ complet. Il ne suffit pas de calculer le terme signé
v̂ᵀRs sans borner le reste.

## 5. Masse assemblée et interfaces partagées

Notons s l'indice d'une sous-structure dans cette section. Les coordonnées
globales conservées sont y. L'application E_s donne les coordonnées
normalisées de son port ; ses déplacements physiques valent W_sE_sy.
Les intérieurs des sous-structures sont disjoints. Les traces partagées
doivent décrire le même déplacement physique après ces transformations.

La norme massique complète s'assemble selon

\[
\|u\|_{M_G}^2
=\sum_s u_s^TM_su_s+y^TM_{\mathrm{ext}}y.
\]

Pour utiliser cette expression comme somme de carrés, on suppose
M_s ≽ 0 et M_ext ≽ 0, avec M_G définie positive sur les coordonnées
physiques conservées. Si M_G n'est que semi-définie, on obtient une
semi-norme qui ne détecte pas toutes les composantes du champ.
`AssemblagePorts` vérifie la taille et la symétrie des matrices externes.
`ControleChamp` ajoute un contrôle numérique de positivité de M_ext et
K_ext ; il ne constitue pas une preuve en arithmétique encadrée. La
positivité des matrices locales complètes M_s reste une hypothèse du
modèle : la seule positivité de leurs blocs M_II ne la démontre pas.

Une masse nodale commune ne doit pas être ajoutée intégralement dans
chaque sous-structure qui touche ce nœud. Il faut distribuer ses
contributions, ou en conserver une seule dans M_ext. Pour des masses
consistantes avec M_IS non nul, la règle correcte est l'assemblage des
matrices complètes, et non un partage des seules diagonales.

### Vérification de la chaîne segmentée actuelle

Dans `chaine_sous_structures.py`, chaque segment reçoit m/2 à chacune
de ses extrémités et m aux nœuds intérieurs. À une interface partagée,
les deux demi-masses donnent bien m. À la racine fixée, la demi-masse
ne contribue pas au mouvement. Au dernier nœud, le supplément externe
m/2 complète la demi-masse du segment. Dans les coordonnées globales
normalisées, ce supplément vaut (m/2)·echelle² sur la dernière diagonale.

Les ressorts appartiennent chacun à un seul segment. La reconstruction
écrit deux fois une interface commune, mais les deux valeurs doivent
coïncider : cela ne justifie ni deux degrés de liberté physiques, ni
deux masses complètes. Un éventuel défaut de coïncidence calculé en
machine doit être mesuré ; il n'est pas couvert par l'hypothèse E_S = 0.

La conversion des forces compte également. Le code utilise
u_port = echelle·y et une force normalisée terminale égale à 1.
Par travail virtuel, la force physique terminale vaut alors 1/echelle.
Une comparaison avec l'oracle de chaîne à force physique unitaire exige
de multiplier la force normalisée par echelle, ou d'adapter la force
demandée à l'oracle.

### Pourquoi garder M_IS dans les Grams de champ

Le terme M_IS disparaît d'une erreur locale à trace nulle. Il reste
présent dans le champ candidat et dans la propagation d'une erreur de
port. Il faut donc calculer Û_sᵀM_sÛ_s avec la masse complète.
Par exemple, pour M = [[1,c],[c,1]] et le vecteur (1,−1), l'énergie
massique vaut 2(1−c), alors que supprimer les termes croisés donne 2.
Le rapport devient arbitrairement grand lorsque c approche 1. Une
interpolation de masse qui ignore M_IS n'a donc pas de justification
uniforme à partir des seules diagonales positives.

## 6. Sous force donnée : erreur de port puis erreur de champ

Soit S_G le Schur assemblé exact dans D d'entrée et Ŝ_G le Schur
renvoyé. L'assemblage des enveloppes locales donne une borne
‖Ŝ_G−S_G‖₂ ≤ B. Ce défaut peut avoir les deux signes lorsque les
écarts de fonctionnel sont inclus.

Pour une force f appliquée aux coordonnées conservées, soit ŷ une
solution approchée et r_y = f−Ŝ_Gŷ son résidu. Si une borne établie
σ_min(Ŝ_G) ≥ σ̂ > B est disponible, alors

\[
\boxed{
\|y-\widehat y\|_2\le b_y
=\frac{\|r_y\|_2+B\|\widehat y\|_2}{\widehat\sigma-B}.}
\]

**Preuve.** La perturbation des valeurs singulières donne
σ_min(S_G) ≥ σ̂−B. Puis
S_G(y−ŷ) = r_y+(Ŝ_G−S_G)ŷ fournit l'inégalité.
La formule de `reponse` correspond au cas r_y = 0 en arithmétique
exacte ; un audit flottant complet doit garder aussi ce résidu.

Pour chaque sous-structure, soit c_M,s la constante locale uniforme de
la section 2. Définissons

\[
A_M=\sum_s c_{M,s}^2E_s^TE_s,\qquad
\eta_{M,\mathrm{op}}=\sqrt{\lambda_{\max}(A_M)},\qquad
\eta_M(\widehat y)=\sqrt{\widehat y^TA_M\widehat y}.
\]

Ces termes bornent les erreurs locales à trace nulle, respectivement
pour tous les ports globaux et pour le port déjà calculé. Des Grams
résiduels directionnels peuvent les remplacer pour réduire le pessimisme.

La norme de propagation du champ candidat est calculable avec de petits
Grams :

\[
\widehat C_M(z)^2
=\lambda_{\max}\left[
M_{\mathrm{ext}}+\sum_sE_s^TÛ_s(z)^TM_sÛ_s(z)E_s\right].
\]

**Proposition 2.** Le champ physique reconstruit û à partir de ŷ et
le champ exact u sous la même force satisfont

\[
\boxed{
\|u-\widehat u\|_{M_G}
\le\eta_M(\widehat y)
+[\widehat C_M(z)+\eta_{M,\mathrm{op}}]b_y.}
\]

**Preuve.** Noter U_G et Û_G les applications globales exactes et
approchées du port vers le champ physique. Écrire

\[
u-\widehat u
=U_G(y-\widehat y)+(U_G-Û_G)\widehat y.
\]

Le second terme a une trace globale nulle. Sa norme massique est donc
la somme des normes intérieures locales, majorée par η_M(ŷ). Les
masses externes n'y contribuent pas. De plus
‖U_G‖ ≤ ‖Û_G‖+‖U_G−Û_G‖ ≤ Ĉ_M+η_M,op.
Appliquer l'inégalité triangulaire donne la proposition. Les deux termes
de la décomposition ne sont pas supposés orthogonaux en masse.

Pour mesurer uniquement les intérieurs, remplacer Ĉ_M² par

\[
\widehat C_{M,I}^2
=\lambda_{\max}\sum_sE_s^TÛ_{I,s}^TM_{II,s}Û_{I,s}E_s.
\]

La même preuve donne une borne dans la norme
(Σ_s e_I,sᵀM_II,se_I,s)^(1/2). Cette norme intérieure ne comprend
ni les ports ni les masses externes. Elle doit être nommée comme telle,
et ne nécessite que la positivité des M_II,s.

### Déformation assemblée

En remplaçant c_M,s par c_D,s, former A_D, η_D,op et η_D(ŷ).
Pour une raideur externe positive semi-définie K_ext, poser

\[
\widehat C_D^2
=\lambda_{\max}\left[K_{\mathrm{ext}}
+\sum_s E_s^T(D_sÛ_s)^T(D_sÛ_s)E_s\right].
\]

Alors

\[
\boxed{\|D_G(u-\widehat u)\|_2
\le\eta_D(\widehat y)+(\widehat C_D+\eta_{D,\mathrm{op}})b_y.}
\]

D_G regroupe les déformations élémentaires et une racine de K_ext.
Si K_ext est indéfinie, ce dernier terme n'est pas une norme de
déformation ; on peut encore contrôler les déformations des éléments
présents dans les D_s, sans nommer l'ensemble énergie positive.

## 7. Quand la borne sous force est-elle uniforme sur la bande ?

La coercivité des intérieurs donne des constantes c_M,s et c_D,s
uniformes. Elle ne donne pas une borne uniforme de l'inverse global.
Il faut en plus établir

\[
\inf_{0\le z\le\Omega^2}\sigma_{\min}(\widehat S_G(z))
\ge\widehat\beta>B.
\]

Si sup ‖f(z)‖₂ ≤ F et sup ‖r_y(z)‖₂ ≤ R_y, alors

\[
Y_{\max}=\frac{F+R_y}{\widehat\beta},\qquad
B_y=\frac{R_y+BY_{\max}}{\widehat\beta-B}
\]

majorent respectivement ‖ŷ‖ et ‖y−ŷ‖. Une enveloppe C̄_M de Ĉ_M
donne donc la borne uniforme explicite

\[
\boxed{\sup\|u-\widehat u\|_{M_G}
\le\eta_{M,\mathrm{op}}Y_{\max}
+(\overline C_M+\eta_{M,\mathrm{op}})B_y.}
\]

Une enveloppe locale de la propagation candidate peut se calculer sans
former une reconstruction dense globale. Si Û_s = X₀,s−J_IV_sY_s
et ‖Y_s‖₂ ≤ y_s,max, alors

\[
\|M_s^{1/2}Û_s\|_2
\le\|M_s^{1/2}X_{0,s}\|_2
+\|M_{II,s}^{1/2}V_s\|_2\,y_{s,\max}.
\]

Avec le membre de droite noté γ_s, on peut choisir
C̄_M² = λ_max(M_ext+Σ_s γ_s²E_sᵀE_s). Les deux normes locales
se calculent avec des Grams ; M_IS reste dans le premier. La même
construction avec D_s donne une enveloppe de déformation.

Une grille de fréquences ne prouve pas le minimum global requis. Une
voie vérifiable consiste à borner ‖Ŝ_G'(z)‖₂ par L sur chaque
intervalle. Pour son centre z_c et sa demi-largeur h,

\[
\sigma_{\min}(\widehat S_G(z))
\ge\sigma_{\min}(\widehat S_G(z_c))-Lh.
\]

Cette inégalité vient de la Lipschitz-continuité des valeurs singulières.
Les dérivées des petites résolvantes et leurs dénominateurs séparés de
zéro permettent de construire L ; les évaluations et arrondis doivent
encore être encadrés pour un certificat machine.

Si la bande contient une résonance globale du modèle conservatif, cette
marge peut être nulle alors que tous les intérieurs sont coercifs.
Une borne uniforme de réponse sous force générale peut alors ne pas
exister. Découper la bande en domaines non résonants est une modification
du domaine annoncé ; ajouter un amortissement est une modification du
modèle et demande une nouvelle preuve. Aucun des deux ne doit être
supposé silencieusement.

## 8. Critère relatif sans connaître la solution exacte

Pour n'importe laquelle des normes précédentes, supposons disposer de
e_abs ≥ ‖u−û‖ et posons a = ‖û‖. Si e_abs < a, l'inégalité
triangulaire donne

\[
\boxed{\frac{\|u-\widehat u\|}{\|u\|}
\le\frac{e_{\mathrm{abs}}}{a-e_{\mathrm{abs}}}.}
\]

Une précision relative τ est donc assurée si
e_abs ≤ τa/(1+τ). Ce critère n'utilise pas le champ exact. Il faut
cependant évaluer a dans la même norme physique : la norme euclidienne
des coefficients réduits ne convient pas. Une erreur observable scalaire
ε_o donne de même ε_o/(|ô|−ε_o), lorsque |ô| > ε_o.
Près d'un zéro de l'observable, une borne absolue peut rester exploitable
alors que le critère relatif ne l'est plus.

Pour une affirmation relative uniforme, il faut une minoration uniforme
de a, pas seulement son évaluation à quelques fréquences. La masse
permet parfois une minoration sans solution intérieure exacte. Pour
chaque sous-structure, son complément de Schur massique est

\[
M_{\mathrm{trace},s}
=W_s^T(M_{SS,s}-M_{SI,s}M_{II,s}^{-1}M_{IS,s})W_s\succeq0.
\]

Minimiser l'énergie massique à trace fixée donne pour tout champ candidat

\[
\|\widehat u\|_{M_G}^2
\ge\widehat y^T\left[M_{\mathrm{ext}}
+\sum_sE_s^TM_{\mathrm{trace},s}E_s\right]\widehat y.
\]

Si la matrice entre crochets majore γ²I avec γ > 0, alors
‖û‖_M ≥ γ‖ŷ‖. Pour une force fixe non nulle et une résolution
réduite exacte, ‖ŷ‖ ≥ ‖f‖/‖Ŝ_G‖ fournit une minoration
uniforme dès qu'une majoration uniforme de ‖Ŝ_G‖ est établie.
Il n'existe pas de minoration analogue générale pour la seule norme
intérieure : des ports peuvent bouger sans entraîner certains intérieurs.

## 9. Enrichissement avec un correcteur statique contrôlé

La réponse ne doit pas être enrichie uniquement parce qu'un résidu de
Schur est petit ou grand. Le choix dépend de la quantité demandée :
masse du champ, déformation, observable ou stabilité du port.

**Proposition 3 : correction statique.** Pour une fréquence fixée, un
candidat à trace correcte peut être corrigé par

\[
\widehat u_I^+=\widehat u_I-K_{II}^{-1}R.
\]

Si e = û_I−u_I est l'erreur, alors

\[
\boxed{e^+=zK_{II}^{-1}M_{II}e.}
\]

L'opérateur K_II⁻¹M_II est autoadjoint dans les métriques K_II et
M_II, avec valeurs propres au plus 1/λ_*. Par conséquent

\[
\boxed{\|e^+\|_{M_{II}}\le\rho\|e\|_{M_{II}},\qquad
\|e^+\|_{K_{II}}\le\rho\|e\|_{K_{II}}.}
\]

La même contraction vaut dans la métrique A(z), puisque dans les
coordonnées propres généralisées les trois métriques sont diagonales
et l'erreur est multipliée par z/λ_j. Cette preuve donne une contraction
uniforme sur la bande, avec une seule résolution statique exacte par
correction. Elle ne donne pas son coût ni la qualité d'une LU approchée.

Partant d'un relèvement X₀ dont le résidu initial vaut R₀−zC, les
corrections après k itérations appartiennent au sous-espace fixe

\[
\operatorname{span}\left\{
(K_{II}^{-1}M_{II})^jK_{II}^{-1}R_0,
(K_{II}^{-1}M_{II})^jK_{II}^{-1}C:
0\le j<k\right\}.
\]

Sa dimension est au plus 2kp, souvent moindre en présence de dépendances
exactes. Un Galerkin exact dans l'espace affine du relèvement enrichi
par ce sous-espace minimise l'erreur A(z). Il fait donc au moins aussi
bien dans cette norme que le k-ième correcteur statique. Si δ₀ majore
le résidu initial en norme duale de masse sur la bande, cela donne

\[
\boxed{\sup\|e_{V}\|_{M_{II}}
\le\frac{\rho^k\delta_0}{\alpha_*},\qquad
\sup\|e_{V}\|_{K_{II}}
\le\frac{\sqrt{\lambda_*}\,\rho^k\delta_0}{\alpha_*}.}
\]

Les normes et δ₀ peuvent ici être prises pour un port imposé ou en norme
d'opérateur sur tous les ports. L'inclusion exacte du sous-espace est une
hypothèse : un seuil de déflation ou des solves inexacts peuvent la perdre.
La convergence est uniforme en fréquence sous les constantes affichées.
Elle n'est pas automatiquement uniforme en nombre de degrés de liberté :
ρ, δ₀ et les métriques peuvent varier avec le maillage ou la partition.
Il ne s'agit pas non plus d'une comparaison de temps avec HCB ou avec
un solveur multicorps complet.

Pour une correction approchée ĉ, définir ξ = R−K_IIĉ. L'identité
devient e⁺ = zK_II⁻¹M_IIe+K_II⁻¹ξ, d'où

\[
\|e^+\|_{M_{II}}
\le\rho\|e\|_{M_{II}}+\frac{\|\xi\|_{M_{II}^{-1}}}{\lambda_*}.
\]

Un défaut uniformément borné par ξ_max crée donc un plancher possible
ξ_max/[λ_*(1−ρ)]. Les arrondis de l'application de la correction
s'ajoutent encore. Le contrôle de ξ dans D d'entrée permet d'adapter
la précision des solves ; accepter un nombre fixe de corrections ne
fournit pas cette garantie.

### Portée de l'adaptation actuellement implémentée

`adapter_champ` évalue les bornes massiques et de déformation aux seules
fréquences fournies. En cas d'échec, il enrichit tous les intérieurs par
rondes, en réutilisant leur factorisation statique et leur base. Il ne
met pas en œuvre un marquage sélectif des sous-structures, et ne démontre
pas que la dimension obtenue est optimale. La contraction de la
proposition 3 décrit un correcteur statique exact de référence ; elle
n'est pas une garantie de contraction de chaque ronde flottante du code.

Les constantes locales utilisées dans ces contrôles sont uniformes sur
la bande. Les valeurs singulières globales, normes de reconstruction et
critères relatifs sont évalués ponctuellement. Une acceptation sur trois
fréquences n'implique donc pas une acceptation sur une grille plus fine,
ni sur toute la bande. Près d'une résonance globale, le terme de défaut
du modèle de Schur amplifié par l'inverse peut dominer malgré un meilleur
champ à port fixé. L'enrichissement seul ne garantit pas de supprimer ce
plancher, et une marge insuffisante doit rester un refus du contrôle.

## 10. Pont conditionnel vers une tangente sous précontrainte

Une extension conservatrice demande la tangente symétrique dans les
coordonnées, ou avec la connexion, appropriées. Supposons qu'elle se
décompose en une partie matérielle et une partie géométrique :

\[
K_{\mathrm{tan}}=D^TD+G,\qquad G=G^T,
\qquad D_I^TD_I\succeq\lambda_D M_{II},\qquad
G_{II}\succeq-\gamma M_{II},\quad\gamma\ge0.
\]

**Proposition 4.** Si λ_eff = λ_D−γ > Ω², alors

\[
A_{\mathrm{tan}}(z)=K_{\mathrm{tan},II}-zM_{II}
\succeq(\lambda_D-\gamma-z)M_{II}\succ0.
\]

La preuve est l'addition des deux minorations. Avec le résidu complet

\[
R=D_I^TDÛ+(GÛ)_I-z(MÛ)_I,
\]

l'identité E_I = A_tan(z)⁻¹R et la borne massique restent valides en
remplaçant λ_* par λ_eff. Les couplages G_IS et M_IS doivent être
inclus dans ce résidu. Pour la déformation matérielle D, la borne est

\[
\boxed{\|DEs\|_2
\le\frac{\sqrt{\lambda_D}}{\lambda_D-\gamma-z}
\|Rs\|_{M_{II}^{-1}}.}
\]

En effet, D_IᵀD_I = A_tan(z)+zM_II−G_II
≼ A_tan(z)+(z+γ)M_II. Les deux bornes résiduelles précédentes
donnent le facteur [α_eff+(z+γ)]/α_eff² = λ_D/α_eff²
pour l'erreur de Gram. Le Schur du fonctionnel de cette même tangente
reste quadratique en l'erreur intérieure. La norme de déformation D
mesure ici la partie matérielle ; elle ne doit pas être renommée énergie
totale positive lorsque la partie géométrique peut être indéfinie.

Cette proposition fournit une condition suffisante à établir pour un
modèle précontraint. Elle ne fournit pas γ, et ne prouve pas qu'une
minoration calculable avec le coût souhaité existe dans tous les cas.
Une tangente de transport non symétrique exprimée directement à gauche
ne relève pas de cette preuve symétrique. Il faut d'abord établir sa
relation avec la Hessienne conservatrice pertinente. De même, corriger
avec (D_IᵀD_I)⁻¹ sous précontrainte donne l'opérateur d'erreur
(D_IᵀD_I)⁻¹(zM_II−G_II) : la seule minoration de G_II ne
démontre pas sa contraction. La proposition 3 ne s'y transpose donc pas
automatiquement.

## 11. Expériences qui peuvent réfuter une précision de champ annoncée

Les témoins suivants doivent distinguer les résultats mathématiques et
les diagnostics flottants, avec les mêmes unités et forces physiques.

- Un intérieur scalaire presque résonant, avec faible couplage au port,
  confronte une petite erreur quadratique de Schur à une grande erreur
  linéaire de champ. Les facteurs α et √ε doivent apparaître.
- Une même assemblée avec masse externe positive accrue peut approcher
  une résonance globale sans changer les coercivités intérieures. Une
  borne de champ sous force doit alors se dégrader ou refuser la marge.
- Une masse complète à M_IS non nul vérifie les Grams de reconstruction
  et la conversion des normes ; sa version diagonalisée artificiellement
  est un autre modèle, à ne pas utiliser comme oracle de substitution.
- Une partition puis un raffinement du nombre de segments vérifient que
  les demi-masses d'interface, le supplément terminal et la conversion de
  force reconstruisent la même chaîne physique. Le champ complet et les
  déformations se comparent à l'oracle, pas seulement le dernier port.
- Un maillage affiné distingue erreur massique intégrée et observable
  ponctuelle. Une constance de la première ne prouve pas celle de la
  seconde lorsque sa norme duale croît.
- Des corrections statiques exactes puis volontairement inexactes
  confrontent la contraction ρ et le terme ξ. Une déflation numérique
  doit être testée avec les directions effectivement supprimées.

Un enrichissement peut cibler le port le plus contributif à A_M ou A_D,
une correction résiduelle chargée par ŷ, ou l'adjoint d'une observable.
Un enrichissement chargé pour une seule force ne couvre pas automatiquement
les autres forces. Si la difficulté vient de la marge globale, améliorer
le seul intérieur le moins cher ne résout pas le problème de stabilité.
Si elle vient d'une marge intérieure faible, retenir explicitement la
direction souple ou modifier la partition peut être nécessaire ; le
nouveau modèle doit alors être réaudité avec ses propres constantes.

Le résultat à démontrer reste une erreur de champ dans une métrique
annoncée, pour une bande et une classe de charges annoncées. Une borne
de Schur seule, une grille de fréquences seule ou un indicateur de norme
de coefficients ne démontre pas cette conclusion.
