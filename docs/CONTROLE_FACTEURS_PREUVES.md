# Contrôle par facteurs : normes d'opérateur et images de réparation

Cette note établit les inégalités des transferts proposés dans le prototype
de contrôle. Les tests Fraction associés vérifient ces inégalités sur de
petites matrices et exhibent leurs hypothèses nécessaires. Ils ne certifient
ni les réponses binary64 ni un gain de temps. Le certificat dirigé d'inertie
du complément conserve sa portée distincte.

## 1. Périmètre qui conserve les champs

On note D l'opérateur de déformation, E l'injection des coordonnées
intérieures, T le relèvement des ports et des directions conservées, W la
base intérieure et Y les coefficients calculés. Le contrôleur forme
X = fl(T − E fl(WY)), Z = fl(DX), MX, le Schur, les solutions par charge et
leurs normes physiques exactement comme auparavant. Seules des **bornes
supérieures de normes d'opérateur** ou d'images de réparation sont remplacées.
Le résidu du petit système Y, le défaut de contrainte, le défaut du
fonctionnel bilinéaire, l'enveloppe anisotrope et la marge du Schur restent
nécessaires. Une compression exacte en arithmétique réelle ne permet pas de
supprimer leurs contreparties calculées.

En particulier, un Gram réduit peut être utile pour majorer la plus grande
valeur propre sans être utilisable pour la norme d'une petite combinaison
signée. Pour Z = [[1,1],[0,ε]], q = (1,−1), ε = 2⁻³⁰, le Gram binary64 peut
être la matrice de uns. Il donne qᵀGq = 0 alors que ‖Zq‖ = ε. Les normes
des champs par charge sont donc conservées dans l'espace physique. Aucune
hypothèse supplémentaire sur une racine de la masse globale, éventuellement
semi-définie positive et couplée, n'est introduite par ces deux transferts.

## 2. Norme d'extension depuis le Gram déjà calculé

Dans cette section Z désigne la matrice **stockée** des images physiques.
Son Gram exact est G = ZᵀZ et le Gram calculé est Ĝ. Sous le modèle habituel
des produits scalaires sans débordement ni sous-flux non comptabilisé,
avec η = γ_N = Nu/(1−Nu) < 1,

    |G − Ĝ| ≤ η |Z|ᵀ|Z|.

La diagonale donne G_ii ≤ Ĝ_ii/(1−η). Choisissons d_i ≥ √(Ĝ_ii/(1−η)).
Cauchy–Schwarz implique (|Z|ᵀ|Z|)_ij ≤ d_i d_j. La matrice symétrique
non négative

    H_ij = min(|Ĝ_ij|, |Ĝ_ji|) + η d_i d_j

domine donc |G| entrée par entrée. Pour tout vecteur w strictement positif,

    ‖Z‖₂² ≤ β(w) := max_i (Hw)_i / w_i.

En effet xᵀGx ≤ |x|ᵀH|x|, puis le rayon spectral de H est majoré par
la norme infinie de diag(w)⁻¹ H diag(w). L'irréductibilité de H n'est pas
nécessaire. Le choix w = 1 donne Gershgorin. Quelques itérations peuvent
proposer un meilleur w, mais seule l'évaluation finale de tous les rapports
est un majorant ; le quotient de Rayleigh d'un vecteur approché n'en est
pas un. Une composante nulle ou négative de w invalide cette preuve.

Une version scalaire plus simple est

    ‖Z‖₂² ≤ ‖sym_exact(Ĝ)‖∞ + η tr(Ĝ)/(1−η).

L'erreur d'une symétrisation calculée doit alors être ajoutée. La version H
évite cette opération et peut être plus précise. La correction du Gram est
indispensable même pour la seule norme d'opérateur : dans l'exemple de la
section 1, les sommes de lignes du Gram arrondi valent 2, mais le quotient
de Rayleigh exact en (1,1) vaut 2 + ε²/2.

Une diagonale nulle ne démontre pas que la colonne physique est nulle :
le carré de 2⁻⁶⁰⁰ s'arrondit à zéro. Si le modèle relatif est invalidé par
le sous-flux, un repli sur une norme de Frobenius physique robuste peut
fournir un majorant d'opérateur. Le Gram seul ne permet plus cette déduction.

## 3. Image de réparation de petit rang, avec résidu direct

La réparation théorique H = CJ provient de J = (BᵀC)⁻¹BᵀW. Le contrôleur
stocke cependant fl(H), puis F_D = fl(D_i fl(H)) et F_M = fl(R_M fl(H)).
Chacune de ces deux matrices stockées F est la cible à majorer. Leur rang
peut dépasser le rang théorique de C.

On prépare une matrice mince Q et un petit facteur L, par exemple depuis
un QR de D_i C ou R_M C puis L = fl(RJ). Leur provenance ne change pas
l'identité exacte entre objets stockés

    F = QL + Δ.

On audite **directement** le résidu Δ, avec e_j ≥ ‖Δ[:,j]‖₂ et
κ ≥ ‖Q‖₂. Pour le vecteur stocké y et v̂ = fl(Ly), si
ε_L(y) ≥ ‖Ly − v̂‖₂, alors

    ‖Fy‖₂ ≤ κ(‖v̂‖₂ + ε_L(y)) + Σ_j e_j |y_j|.

Pour dominer aussi l'ancien produit calculé fl(Fy), il suffit d'ajouter
γ_k Σ_j f_j |y_j|, où f_j ≥ ‖F[:,j]‖₂, sous le même modèle relatif.
Cette preuve couvre les défauts de formation de H, des images et du QR
sans supposer que fl(H) a un rang exact donné. Les deux chemins doivent
utiliser le même y stocké, ici fl(Yq) ; viser Yq exact demanderait aussi
le budget de formation de ce produit.

Contre-exemple de rang : ε = 2⁻²⁷, C = (1,1+ε)ᵀ et J = (1,1−ε).
Le produit exact CJ a rang un. Son coefficient inférieur droit 1−ε²
s'arrondit à 1 ; le déterminant du produit stocké devient ε². Pour
y = (1−ε,−1), Jy = 0 mais fl(CJ)y = (0,−ε²)ᵀ. Le résidu direct est
indispensable pour cette charge.

## 4. Auditer les résidus calculés et l'orthogonalité

Toutes les matrices ci-dessous sont des objets stockés. Si
Δ̂ = fl(F − fl(QL)), avec a colonnes contractées dans QL, on a

    |(F−QL) − Δ̂| ≤ γ_a |Q||L| + u/(1−u) |Δ̂|.

Le dernier terme dépend du **résultat de la soustraction**. Il ne peut
être supprimé parce que le résidu observé semble petit. Un résidu nul
observé ne rend pas le premier terme nul : pour Q = 1+2⁻²⁷,
L = 1−2⁻²⁷ et F = 1, le produit arrondi vaut 1 mais le résidu exact
vaut 2⁻⁵⁴.

Si Ω̂ = fl(fl(QᵀQ)−I), un audit suffisant est

    ζ ≥ ‖Ω̂‖₂ + ‖γ_N |Q|ᵀ|Q| + u/(1−u)|Ω̂|‖F,
    κ = √(1+ζ).

La soustraction de I peut être exacte par Sterbenz lorsque ses hypothèses
sont vérifiées. En revanche Q = (1,2⁻³⁰)ᵀ possède un Gram calculé égal
à 1 et un Gram exact strictement supérieur à 1. Ignorer l'erreur du
produit d'orthogonalité sous-estime κ.

Les normes employées dans ces audits doivent elles-mêmes être majorées.
Pour une colonne (1,2⁻²⁷,…,2⁻²⁷) de 128 lignes, une sommation séquentielle
des carrés peut donner 1 alors que le carré exact vaut 1+127·2⁻⁵⁴.
Un QR audité avec des normes sous-estimées n'est donc pas un certificat.
Le rescaling par le maximum absolu évite les carrés extrêmes ; il faut
encore couvrir les erreurs des divisions, carrés, sommation, racine et
remise à l'échelle. Un modèle mixte |fl(a)−a| ≤ u|a|+τ/2, τ étant le
plus petit sous-normal, apporte le budget absolu absent du modèle γ seul.
Ainsi, pour une colonne non nulle, m = max |v_i|, a_i = fl(v_i/m)
et q = fl(Σ fl(a_i²)), un budget conservateur est

    ‖v‖₂² ≤ m² [q+(3n+2)τ] / [1−γ_(n+3)].

On exige γ_(n+3) < 1. Les trois erreurs relatives par composante sont
celles de la division utilisée deux fois dans le carré et du carré ;
les additions apportent au plus n erreurs supplémentaires. Les termes
absolus du modèle mixte, sur des composantes normalisées de module au
plus 1, sont couverts par (3n+2)τ. Ce budget doit être évalué vers le
haut, ainsi que sa racine et la remise à l'échelle. La colonne nulle
est traitée exactement et le cas non représentable doit être refusé.
Plus explicitement, si z_i = fl(a_i²), le modèle mixte et |a_i| ≤ 1
donnent (v_i/m)² ≤ [z_i+(5/2)τ]/(1−u)³. Une sommation de n termes
non négatifs donne Σz_i ≤ [q+nτ/2]/(1−u)ⁿ. Leur combinaison est
majorée par [q+3nτ]/(1−u)^(n+3), puis par l'expression ci-dessus,
puisque (1−u)^(n+3) ≥ 1−(n+3)u ≥ 1−γ_(n+3).
Les inflations binary64 et nextafter du prototype restent des **audits
évalués**, tant que chaque opération de la preuve n'est pas encadrée.

## 5. Alternative QR de l'extension complète

Le QR préparé de F = [Â, B̂], Â = fl(DT), B̂ = fl(D_i W), donne
F = QR+Δ avec erreurs par colonnes δ_T, δ_W et κ ≥ ‖Q‖₂. Si

    L = R_T − R_W Y,  L̂ = fl(R_T − fl(R_W Y)),

alors ε_L ≥ ‖γ_k |R_W||Y| + u/(1−u)|L̂|‖F suffit pour la petite
formation. On obtient

    ‖fl(DX)‖₂ ≤ κ(‖L̂‖₂+ε_L)
                 + ‖δ_T+δ_W|Y|‖₂ + ε_formation(Y).

Aucun rang plein de F n'est nécessaire. Une troncature doit placer la
partie supprimée dans Δ. Cette alternative évite les grandes SVD mais
ajoute une préparation ; les preuves ne décident pas de son intérêt temporel.

Voici un budget de formation préparé par colonnes. Soient γ_D et γ_i
des constantes couvrant les longueurs des produits scalaires de D et D_i,
γ_k celle du produit WY et c_k = γ_k + u(1+γ_k). Définissons

    A₊ = |D||T|,  A_i₊ = |D_i||T_i|,  B₊ = |D_i||W|,
    H₀ = 2γ_D A₊ + u(1+γ_D) A_i₊,
    H₁ = [γ_D(1+c_k)+c_k+γ_i] B₊.

Pour X_i = fl(T_i−fl(WY)), ports copiés exactement,

    |fl(DX) − (Â−B̂Y)| ≤ H₀+H₁|Y|.

La preuve utilise successivement
|X_i−(T_i−WY)| ≤ u|T_i|+c_k|W||Y|, l'erreur γ_D|D||X| du
produit final, puis les erreurs de préparation de Â et B̂. Si h₀ et h₁
majorent les normes euclidiennes de leurs colonnes, on peut choisir
ε_formation(Y) = ‖h₀+h₁|Y|‖₂. Le coût en ligne est indépendant du
nombre de lignes physiques. Employer seulement une borne y_max de la
solution exacte Y est incorrect pour la solution calculée ; l'expression
en |Y| stocké évite cette hypothèse, sinon il faut y_max + erreur_Y.

## 6. Limite d'un transfert futur du Schur

Remplacer seulement des majorants laisse le Schur stocké identique. Un
futur Schur calculé depuis les facteurs nécessite une erreur supplémentaire.
Si Z = QL+Δ et ‖QᵀQ−I‖₂ ≤ ζ, alors

    ‖ZᵀZ−LᵀL‖₂ ≤ ζ‖L‖₂²
                    +2√(1+ζ)‖L‖₂‖Δ‖₂+‖Δ‖₂².

Il faut encore ajouter les erreurs des petits produits, de la soustraction,
de la partie massique pondérée par ω² et de sa racine éventuelle. La marge
globale devient σ_min(S_stocké) moins **toutes** ces erreurs et les défauts
déjà présents. Près d'un pôle global, même une petite erreur supplémentaire
peut supprimer l'acceptation. Une identité théorique de Gram ne justifie
donc pas le transfert de la résolution ou des normes des petites charges.
