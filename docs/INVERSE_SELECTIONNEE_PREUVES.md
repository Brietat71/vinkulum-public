# Minoration spectrale par inverse sélectionnée et écart énergétique encadré

Carnet du 7 septembre 2026. Les implémentations sont
`python/vinkulum/_ports/inverse_selectionnee.py` et
`python/vinkulum/_ports/certificat_spectral.py` ; les modules `ci/` homonymes
conservent les imports des expériences antérieures. Le premier
calcule une trace avec la récurrence classique de Takahashi. Le second
encadre cette récurrence et l'écart au facteur d'énergie fourni pour
obtenir une minoration de coercivité sur les valeurs binary64 d'entrée.
Ce sont deux contrats distincts : le résultat flottant du premier ne
devient pas un certificat par l'ajout d'une marge arbitraire.
L'[API publique de réduction 0.10.0](REDUCTION_PORTS.md) utilise le second
quand l'appelant ne fournit pas de constante spectrale.

## 1. Domaine et raison d'utiliser une trace

Soit D_I une matrice réelle rectangulaire d'énergie intérieure, M une
masse symétrique définie positive, E une diagonale strictement positive
et R une matrice triangulaire supérieure de diagonale positive. Posons

\[
K_D=D_I^TD_I,\qquad K_Q=E^{-1}R^TRE^{-1},\qquad
S=(R^TR)^{-1}.
\]

Dans une factorisation QR exacte de D_IE, on aurait K_D = K_Q. Le
prototype ne suppose pas cette égalité pour le facteur calculé en machine.
L'ordre de R est celui des coordonnées intérieures ; toute permutation
de coordonnées doit aussi être appliquée à M, E et D_I.

La trace de flexibilité massique est

\[
\tau=\operatorname{tr}(M K_Q^{-1})
=\operatorname{tr}((EME)S)>0.
\]

Si λ_j sont les valeurs propres généralisées de (K_Q,M), alors

\[
\tau=\sum_j\lambda_j^{-1},\qquad
\boxed{\lambda_{\mathrm{trace}}=\tau^{-1}
\le\lambda_{\min}(K_Q,M).}
\]

**Preuve.** Les valeurs propres de M¹/²K_Q⁻¹M¹/² sont positives
et valent λ_j⁻¹. Leur maximum est inférieur à leur somme. L'inégalité
est une minoration de coercivité, pas une approximation précise du
premier mode. Elle peut perdre un facteur aussi grand que le nombre
de coordonnées lorsque les valeurs propres sont toutes comparables.
Le coût évité est celui d'une inverse complète ou d'une recherche de
mode ; la largeur de bande exploitable peut en contrepartie diminuer.

La trace demande uniquement S_ii pour une masse diagonale. Pour une
masse symétrique générale, elle demande

\[
\tau=\sum_i e_i^2M_{ii}S_{ii}
+2\sum_{i<j\,;\,M_{ij}\ne0}e_ie_jM_{ij}S_{ij}.
\]

Omettre une entrée hors du motif de R change cette trace. Un bloc
d'inertie 3×3 non diagonal ne peut donc pas être traité comme trois
masses indépendantes sous prétexte que seules les diagonales de
l'inverse sont immédiatement disponibles.

## 2. Récurrence exacte et fermeture des dépendances

Noter N_i = {k > i : R_ik ≠ 0}. Puisque RS = R⁻ᵀ,

\[
\boxed{S_{ij}=-R_{ii}^{-1}\sum_{k\in N_i}R_{ik}S_{kj},
\quad j>i,}
\]

\[
\boxed{S_{ii}=R_{ii}^{-2}
-R_{ii}^{-1}\sum_{k\in N_i}R_{ik}S_{ki}.}
\]

La première formule vient du zéro de R⁻ᵀ au-dessus de la diagonale ;
la seconde vient de sa diagonale R_ii⁻¹. On calcule i en ordre
décroissant, les entrées hors diagonale de la ligne i avant S_ii.
Pour k,j > i, l'entrée symétrique S_kj appartient à une ligne dont
l'indice est strictement supérieur à i : elle est déjà calculée si
le motif est fermé.

Le motif initial contient toutes les diagonales, toutes les paires
(i,k) avec k dans N_i, et toutes les paires demandées par M. Pour
chaque paire (i,j), j > i, il faut ajouter les paires non orientées
(k,j), k dans N_i. Après réorientation par le plus petit indice,
leur ligne est strictement supérieure à i. Un passage croissant sur
i propage donc toutes les dépendances ; les ajouts ne nécessitent
pas de revenir à une ligne déjà traitée. Les diagonales ajoutées
étaient déjà présentes.

**Proposition 1.** Sur le motif ainsi fermé, la récurrence retourne
les entrées exactes de S en arithmétique exacte.

**Preuve.** La dernière diagonale vaut R_nn⁻². Supposons toutes
les lignes strictement supérieures à i correctes. La fermeture fournit
chaque terme de la première formule, donc les S_ij sélectionnés sont
corrects. Les S_ik nécessaires à la seconde sont inclus dans le motif
initial ; S_ii est donc correct également. L'induction descendante
achève la preuve. Aucune entrée inconnue n'est remplacée par zéro.

Pour le motif structurel rempli d'un facteur de Cholesky, les voisins
supérieurs de chaque sommet forment une clique ; le motif usuel de
l'inverse sélectionnée est déjà fermé. La fermeture explicite reste
utile pour un R quelconque, pour des zéros numériques dans le facteur,
ou pour des entrées de masse demandées hors de ce motif.

## 3. Budgets et portée algorithmique

Soit J_i l'ensemble des colonnes hors diagonale sélectionnées dans la
ligne i. Le nombre de coefficients stockés vaut

\[
n+\sum_i|J_i|,
\]

et les visites symboliques et sommes numériques sont contrôlées par

\[
O\left(n+\sum_i |N_i|(1+|J_i|)\right).
\]

Le code compte les visites de dépendances, les opérations scalaires
de la récurrence et les coefficients sélectionnés. Ces compteurs ne
sont ni des mesures de mémoire en octets ni un décompte matériel exact
des opérations. La validation des entrées lit aussi leurs coefficients.
Une bande et des blocs physiques de taille bornée peuvent ainsi donner
un coût linéaire en nombre de corps. Pour un graphe arbitraire, le
remplissage peut devenir quadratique en stockage et cubique en travail.
La sélection ne démontre donc pas un coût linéaire général.

Le prototype refuse explicitement un dépassement de budget ; il ne
remplace pas ce refus par une inverse dense ou n résolutions triangulaires.
Une entrée non sélectionnée demandée après le calcul est refusée.
Des masses hors motif sont traitées par fermeture, dans ces mêmes budgets.
Leur positivité n'est prise en charge que pour des blocs connexes de
taille bornée, six par défaut. Des blocs 3×3 d'inertie et des translations
diagonales entrent dans ce domaine, y compris après une permutation.
Une masse consistante connectant une grande partie du maillage requiert
un autre procédé de validation ou un refus : elle n'est pas couverte
silencieusement par l'API actuelle.

## 4. Pourquoi une trace QR flottante ne suffit pas pour D d'entrée

La minoration 1/τ vise K_Q. Le facteur calculé par QR n'est pas
exactement un facteur de D_I au sens des nombres réels représentés par
ses entrées. Un eigensolve convergé ou une multiplication du résultat
par 0,999 ne démontre pas la relation d'ordre nécessaire sur K_D.

Posons ΔK = K_D−K_Q. Si un majorant établi satisfait

\[
\|K_Q^{-1/2}\Delta K K_Q^{-1/2}\|_2\le\eta<1,
\]

alors

\[
K_D\succeq(1-\eta)K_Q.
\]

Si τ̄ ≥ τ est également établi, on obtient

\[
\boxed{\lambda_{\min}(K_D,M)
\ge\frac{1-\eta}{\overline\tau}>0.}
\]

Cette assertion prouve aussi le rang colonne plein de D_I ; celui-ci
n'a pas à être préalablement affirmé à partir d'un seuil numérique.
Si η ≥ 1, la méthode ne conclut pas. Cela ne signifie pas que D_I
est nécessairement singulière, seulement que cette estimation de son
écart relatif est insuffisante.

## 5. Majorant d'écart utilisant seulement les diagonales sélectionnées

Soit d_i = (K_Q⁻¹)_ii = e_i²S_ii. Pour les vecteurs canoniques v_i,

\[
K_Q^{-1/2}\Delta K K_Q^{-1/2}
=\sum_{ij}\Delta K_{ij}(K_Q^{-1/2}v_i)(K_Q^{-1/2}v_j)^T.
\]

La norme spectrale d'un produit extérieur uvᵀ vaut ‖u‖₂‖v‖₂.
L'inégalité triangulaire fournit donc

\[
\boxed{\eta\le\sum_{ij}|\Delta K_{ij}|\sqrt{d_i d_j}.}
\]

Le membre de droite est un choix admissible de majorant, souvent
pessimiste : il détruit les annulations entre contributions. Pour un
stockage du triangle supérieur, les contributions i < j doivent être
comptées deux fois. Toutes les d_i sont déjà calculées par la sélection.

Le calcul de ΔK additionne les produits extérieurs de chaque ligne de
D_I puis soustrait ceux de RE⁻¹. Il ne passe pas par le Gram calculé
en float. Les nombres de produits et de paires distinctes sont bornés
explicitement. Le travail de cette étape dépend de la somme des carrés
des nombres de coefficients par ligne. Il peut donc être élevé même
si le nombre total de coefficients du facteur paraît modéré.

Ce majorant emploie les diagonales de K_Q⁻¹ ; il n'exige pas que
l'ensemble des entrées de ΔK fasse aussi partie de l'inverse sélectionnée.
C'est précisément ce qui évite une seconde fermeture sur le motif du
Gram d'écart.

## 6. Encadrement dirigé : ce que le certificat établit

`certifier_spectral` interprète chaque valeur binary64 d'entrée comme
un nombre réel exact, obtenu par `Decimal.from_float`. Il recommence
la récurrence sélectionnée avec des intervalles, en utilisant des
contextes privés pour les arrondis vers −∞ et +∞. Le contexte Decimal
global de l'appelant n'est ni utilisé pour les opérations calculées,
ni modifié.

Les entrées exactes de S appartiennent aux intervalles calculés par
induction dans la proposition 1 : chaque somme, produit et quotient
est arrondi vers l'extérieur. La trace, ses termes de masse croisés et
les diagonales physiques e_i²S_ii sont ensuite encadrés de la même
façon. La construction par produits extérieurs encadre chaque entrée
de ΔK dans le modèle D d'entrée.

Le calcul de η utilise la valeur absolue maximale de chaque intervalle
ΔK_ij et les majorants des d_i. Il majore les produits, racines et
sommes. Une subtilité d'implémentation est traitée explicitement :
`Decimal.sqrt` est correctement arrondi au plus proche, même lorsque
le contexte demande un arrondi vers +∞. Le code prend son successeur
Decimal, qui est un majorant. Il ne suppose pas que changer l'arrondi
du contexte change la règle de `sqrt`.

Pour la masse, le contrôle `eigvalsh` flottant du premier prototype
n'est pas accepté comme preuve. Chaque petit bloc est refactorisé
par LDL en intervalles ; les bornes inférieures de tous les pivots
doivent être strictement positives. Cela démontre sa positivité exacte.
Les entrées doivent être exactement symétriques. Les données sparse
non canoniques sont refusées avant toute somme implicite de doublons
pour D, R et M ; sinon la conversion pourrait changer le modèle par
un arrondi qui ne figure pas dans le certificat.
Le contrôle est répété après conversion en CSR : certains formats,
dont LIL, n'exposent pas ce renseignement avant conversion. L'API prenant
un objet d'inverse sélectionnée reconstruit sa sélection depuis les R,
M et E actuels. Elle ne fait pas confiance à des diagonales, voisins ou
résultats mémorisés devenus incohérents après une mutation de l'objet.
Les indices compressés ou COO sont inspectés directement : un drapeau
`has_canonical_format` peut lui aussi être périmé après une mutation.

La soustraction 1−η̄ et la division par τ̄ sont arrondies vers le
bas. La conversion finale en float est également dirigée vers le bas
par comparaison au Decimal obtenu et, si nécessaire, par son voisin
float inférieur. L'API renvoie le minorant float, ses valeurs Decimal,
les majorants η̄ et τ̄, les budgets consommés et la portée du résultat.

Le statut `certification_machine=True` concerne seulement cette
minoration spectrale des valeurs d'entrée. Il ne certifie pas :

- la fidélité de D à une poutre continue, à une géométrie mesurée ou à
  des paramètres constitutifs incertains ;
- l'égalité de DᵀD exact avec un assemblage natif de K ayant suivi
  un autre chemin d'arrondi ;
- les autres opérations flottantes de réduction, d'assemblage, de
  reconstruction ou de contrôle relatif du champ ;
- une tangente sous précontrainte K_tan = DᵀD+G. Il faut alors une
  minoration supplémentaire de G dans la métrique massique ;
- une marge globale uniforme près d'une résonance de l'assemblage.

Une précision Decimal fixe peut ne pas suffire sur une entrée très
mal conditionnée. Les intervalles peuvent croître, un pivot peut ne
pas être séparé de zéro ou η̄ peut dépasser un. Le comportement
correct est alors un refus, ou un nouveau calcul explicitement budgété
avec une précision supérieure. La complétude du procédé n'est pas
revendiquée. Le prétraitement flottant peut également refuser avant
qu'un encadrement à plus grande précision ait été tenté.

## 7. Témoins persistants

`ci/test_inverse_selectionnee.py` compare de petites matrices SPD
rectangulaires à une inverse dense indépendante, avec échelles très
différentes, masses par blocs, permutations, facteur non chordal et
couplage de masse hors motif. Un témoin de grande bande interdit la
densification et les résolutions denses ; sa diagonale inverse possède
une récurrence scalaire indépendante. Les budgets et entrées hors du
contrat doivent produire des refus explicites.
Un témoin de console à 512 éléments utilise D réellement exporté par
`Noyau.facteurs_materiels_poutres`, avec racine supprimée et six ports
au bout. Il vérifie la chaîne complète QR, sélection et certificat,
en interdisant la densification des grandes matrices. La raideur native
assemblée par un autre chemin d'arrondi n'est pas son modèle certifié.

Le certificat est confronté à des calculs exacts `Fraction` sur une
petite entrée perturbée : la trace rationnelle exacte doit appartenir
à son intervalle ; les pivots LDL rationnels de K_D−λ_certM doivent
être positifs. Un contexte Decimal ambiant de faible précision et
d'arrondi différent doit laisser inchangés les résultats du certificat.
Des carrés rationnels exacts vérifient les majorants de racines. Un
facteur incompatible avec D doit être refusé lorsque η̄ ne conclut pas.
Ces tests vérifient l'implémentation ; la justification de la minoration
reste celle des propositions précédentes.

## 8. Antériorité et sources primaires

L'inverse sélectionnée n'est pas une nouveauté de Vinkulum. La
récurrence, ses versions par blocs et la fermeture sur un graphe rempli
s'inscrivent dans les travaux de Takahashi et collaborateurs (1973),
Erisman–Tinney (1975), puis les développements modernes de SelInv.
Lin et collaborateurs donnent la dérivation par compléments de Schur
et la propriété de fermeture dans
[SelInv — An Algorithm for Selected Inversion of a Sparse Symmetric
Matrix (2011)](https://web.stanford.edu/~lexing/diagex.pdf), §2.

Le découpage entre sélection demandée, fermeture symbolique et calcul
numérique par blocs apparaît encore dans
[Abdul Fattah, Ltaief, Rue et Keyes, GPU-Accelerated Parallel Selected
Inversion for Structured Matrices Using sTiles (2025)](https://arxiv.org/pdf/2504.19171),
§III. Leur étude de matrices structurées ne démontre pas un faible coût
universel pour des demandes hors motif ou pour un graphe mécanique
arbitraire, et ne fournit pas le certificat d'arrondi du présent code.

La distinction entre une factorisation approchée et une conclusion
spectrale vérifiée est cohérente avec
[Rump et Ogita, Verified Error Bounds for Matrix Decompositions
(2024)](https://www.tuhh.de/ti3/paper/rump/RuOg24a.pdf), §5.
Cette source emploie notamment des inverses approchées ; elle ne prouve
pas le coût creux du prototype ici. Le majorant par diagonales de la
section 5 et son utilisation dans le certificat sont démontrés dans
ce carnet, sans revendication de priorité bibliographique.
