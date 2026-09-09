# Contrôle réduit : facteurs d’image, annulations et garanties numériques

Recherche primaire ciblée, 8 septembre 2026 ; code lu à HEAD 724236e dans /tmp/vinkulum-energie-native-travail. Aucun solveur ni banc exécuté. Les six références ci-dessous suffisent pour identifier l’antériorité, le transfert algorithmique et les limites. Le principe de compression stable des résidus était déjà explicite en 2014 ; l’apport possible de Vinkulum est son intégration aux contraintes, aux normes physiques et aux décisions de refus, pas la revendication d’un principe inédit.

## Diagnostic et choix immédiat

Dans ci/controle_complement.py, reponses, lignes 129–203, chaque fréquence reconstruit X=T−E_i WY, puis DX, MX, le Schur S=(DX)ᵀDX−ω²XᵀMX. Les normes de champs sont formées après reconstruction physique ; les images de réparation H sont encore appliquées en grande dimension. La préparation de docs/KRYLOV_CONTRAINT_COMPRESSION.md traite un autre coût : certifier l’enveloppe des résidus de préparation.

L’accélération mathématiquement pertinente consiste à factoriser une image commune une fois, puis à appliquer ses petits facteurs aux coefficients avant de calculer une norme. Une image H=FJ de rang contraint petit offre un premier transfert particulièrement ciblé : QR de F, action RJv, majoration du défaut H−QRJ. Cela évite de répéter H v en grande dimension sans tronquer un rang numérique arbitrairement.

La branche nouvelle ControleFacteurs(extension="gram"), lue après le diagnostic, prend une voie plus limitée : elle conserve le Schur et les normes des champs historiques, majore seulement la plus grande valeur propre du Gram de DX et comprime les images H. Les critiques des petites formes quadratiques ne réfutent pas cet usage du Gram. Ses corrections sont calculées en binary64 ; le module déclare correctement certification_machine=False. Aucun résultat de vitesse n’est déduit des articles.

## Six sources primaires vérifiées

**1. Buhr, Engwer, Ohlberger, Rave, A numerically stable a posteriori error estimator for reduced basis approximations of elliptic equations.** Prépublication arXiv déposée le 30 juillet 2014 ; contribution WCCM 2014 répertoriée par l’équipe de Münster. [Notice primaire](https://arxiv.org/abs/1407.8005), [PDF, 7 pages](https://arxiv.org/pdf/1407.8005), [liste institutionnelle](https://www.uni-muenster.de/AMM/ohlberger/team/stephan_rave.shtml).

Pages 2–3 : problème hilbertien coercif, décomposition affine et borne résiduelle du théorème 4.1. Pages 3–4, équations (7)–(10) : la formule quadratique préassemblée du résidu peut rencontrer un plancher d’arrondi de l’ordre de √u ; une base orthonormale de l’image résiduelle permet de former d’abord sa combinaison linéaire, puis sa norme. L’algorithme 1 réorthogonalise ; le §4.3 donne le coût supplémentaire de préparation et un coût en ligne indépendant de la dimension physique. C’est le raccord direct. Ce n’est pas une preuve de certification par intervalles de toute l’implémentation ; le coût des solveurs de Riesz annoncé suppose des préconditionneurs adaptés.

**2. Casenave, Ern, Lelièvre, Accurate and online-efficient evaluation of the a posteriori error bound in the reduced basis method.** ESAIM:M2AN 48(1), 207–229, publication en ligne le 10 janvier 2014. [DOI](https://doi.org/10.1051/m2an/2013097), [PDF auteur de la version publiée](https://cermics.enpc.fr/~ern/PDFs/14_CEL_M2AN.pdf).

Définitions 2.3–2.4, p. 210 : indépendance en ligne vis-à-vis de N et dépendance affine. Proposition 3.2, p. 213 : analyse de pires erreurs d’arrondi expliquant les seuils proportionnels à u et √u, selon la formule. Proposition 3.7, p. 214 : les erreurs des solutions préparatoires imposent leur propre seuil ; réécrire la norme ne les efface pas. Proposition 4.5, p. 218 : évaluation par interpolation empirique bien définie et efficace en ligne. Proposition 4.9, p. 220 : stabilisation itérative équivalente en arithmétique exacte. Pour Vinkulum, l’interpolation est moins immédiate qu’une image affine explicite de quelques dizaines de colonnes : son défaut devrait être borné avant de devenir un majorant uniforme.

**3. Rump et Ogita, Verified Error Bounds for Matrix Decompositions.** SIAM Journal on Matrix Analysis and Applications 45(4), 2155–2183, publié en ligne le 11 novembre 2024. [Éditeur/DOI](https://epubs.siam.org/doi/10.1137/24M165096X), [manuscrit auteur, 30 pages](https://www.tuhh.de/ti3/paper/rump/RuOg24a.pdf).

Le §5, p. 20 du manuscrit (page PDF 19), étudie une matrice haute de plein rang colonne. Avec une approximation de l’inverse de R, il ramène la vérification à une matrice proche d’une isométrie ; les produits associés évitent de former directement AᵀA. Le lemme 5.1 contrôle une isométrie à partir d’un défaut d’orthogonalité et d’un écart matriciel. Cette voie vérifie des facteurs, mais l’hypothèse de plein rang interdit son application littérale à toute image commune dépendante. Pour notre besoin de normes, borner directement A−QR et QᵀQ−I peut suffire, sans inverser R ni sélectionner un rang. Une certification vérifiée supplémentaire doit être chiffrée dans la préparation.

**4. Lange et Rump, Accurate floating-point matrix residuals.** Manuscrit de 50 pages annoncé soumis pour publication en 2026 ; aucune acceptation vérifiée. [Annonce officielle](https://www.tuhh.de/ti3/publications.shtml?author=rump), [PDF auteur](https://www.tuhh.de/ti3/paper/rump/MatrixResidualsFinal.pdf).

Lemme 4 et corollaire 5, p. 7 : exactitude de produits de parties flottantes sous conditions de grille, d’échelle et de taille des sommes, avec absence de débordement. Le §8.4, pp. 43–45, exploite les produits réellement non nuls : l’équation (41) contraint le découpage à partir de leur nombre maximal. Transfert possible : évaluer précisément les défauts de facteurs et les produits critiques de petits résidus, avec une précision accrue là où elle réduit effectivement la borne. Cela ne supprime pas à soi seul la dépendance en N. Les résultats rapides sur produits creux dépendent de la répartition des valeurs entre parties ; ils ne constituent pas une promesse universelle de coût ni un certificat obtenu simplement en augmentant la précision.

**5. Humphry et Yano, Efficient Hyperreduction for Large-Scale Problems: Exploiting Reducible Constraint Manifolds in Empirical Quadrature Procedure.** International Journal for Numerical Methods in Engineering 126(23), e70204 ; première publication le 7 décembre 2025. [Article éditeur](https://onlinelibrary.wiley.com/doi/10.1002/nme.70204), [liste institutionnelle de l’auteur](https://arrow.utias.utoronto.ca/~myano/publications).

§4.2, équations (31)–(33) : actualisation QR lors de l’ajout de colonnes. §5, équation (34) : pour la solution de moindres carrés Rρ=Qᵀb, le résidu s’écrit (I−QQᵀ)b ; cette expression est utilisée sélectivement lorsque le petit résidu devient difficile à évaluer directement. Le transfert est une organisation des résidus par projection et un traitement particulier des fins de convergence. Pour un Y seulement approché, supprimer son résidu de résolution serait incorrect. L’article ne fournit pas une garantie machine générale d’exactitude relative de ce projecteur ; sa méthode d’hyperréduction introduit par ailleurs une approximation inutile pour notre première image affine explicite.

**6. Smetana, Zahm, Patera, Randomized Residual-Based Error Estimators for Parametrized Equations.** SIAM Journal on Scientific Computing 41(2), A900–A926, publié le 28 mars 2019 ; prépublication 2018. [DOI éditeur](https://epubs.siam.org/doi/abs/10.1137/18M120364X), [prépublication et PDF lus](https://arxiv.org/abs/1807.10489).

Proposition 2.1 et corollaire 2.2, pp. 5–6 de la prépublication : projections gaussiennes et probabilité de contrôle uniforme sur un ensemble fini de vecteurs. Proposition 3.1 et corollaire 3.2, pp. 10–11 : l’approximation des problèmes duaux doit être contrôlée ; une marge positive sépare le dénominateur de zéro. Une covariance définie positive permet de choisir la norme. C’est une possibilité de présélection pour des images beaucoup plus vastes, pas un remplacement du contrat déterministe actuel. Des charges adaptées après tirage ou toutes leurs combinaisons nécessitent un contrôle du sous-espace, au-delà d’une union finie sur des vecteurs fixés.

## Dérivation propre : voie plus ambitieuse, à mesurer séparément

Posons A=[T,E_iW], C=[I;−Y] ; alors X=AC. Pour les essais actuels, A possède 31/39/47 colonnes. Si M=LᵀL, préparer

    DA=Q_D R_D,    LA=Q_M R_M,    G_D=R_D C,    G_M=R_M C.

En arithmétique exacte, les normes d’extension valent ||G_D||₂ et ||G_M||₂, le Schur vaut G_DᵀG_D−ω²G_MᵀG_M, et les normes de la charge j valent ||G_D q_j||₂ et ||G_M q_j||₂. On applique les facteurs AVANT la norme, sans former qᵀCᵀRᵀRCq. Une QR séparée des petites images H protège leur échelle. Une QR mince sans troncature fonctionne aussi quand R est singulier : aucun inverse n’est nécessaire.

Le travail en ligne des métriques dépend alors des petites dimensions ; retourner les champs physiques exige toujours leur reconstruction en N. La préparation QR, les audits et la factorisation massique sont à compter. Un facteur de M_ii seul ne représente pas une masse complète comportant des couplages port–intérieur. Une masse complète semi-définie positive peut nécessiter un facteur rectangulaire ; une norme relative sur sa direction nulle reste indéfinie.

Pour une image exacte A₀ et des facteurs stockés, supposons ||A₀−QR||₂≤ε et ||QᵀQ−I||₂≤η<1. Alors

    max(0, √(1−η)||Rc||₂−ε||c||₂) ≤ ||A₀c||₂
      ≤ √(1+η)||Rc||₂+ε||c||₂.

Le défaut ε doit inclure la formation initiale de l’image. Les produits en ligne et la reconstruction physique ajoutent leurs propres erreurs. Cette estimation reste additive : aucune réécriture ne garantit une précision relative uniforme pour une combinaison presque nulle.

Pour le Schur, notons G=RC, e=ε||C||₂ et β=√(1+η). La différence entre le Gram physique et GᵀG est majorée par

    δ = η||G||₂² + 2β||G||₂ e + e².

Il faut sommer δ_D+ω²δ_M et l’erreur de formation du petit Schur ; un défaut de facteur massique ajoute ||M−LᵀL||₂||X||₂². Ce supplément réduit la marge sur σ_min(S) et majore l’action résiduelle par δ_S||q||₂. Pour une borne relative, on utilise une MINORATION de la norme du champ candidat dans b/(n_min−b), avec refus si le dénominateur n’est pas positif.

## Contre-épreuves et recommandation bornée

Une contre-épreuve exacte minuscule suffit à rejeter la stratégie des petits Grams de réponses : A=((1,1),(0,2⁻³⁰)), c=(1,−1). La norme exacte vaut 2⁻³⁰ ; le Gram binary64 peut perdre le terme 2⁻⁶⁰ ajouté à 1 et annoncer zéro. Ce défaut ne concerne pas la majoration de la plus grande valeur propre employée par la nouvelle branche.

Autres contrôles ciblés : image de rang déficient sans inversion de R ; facteur Q non isométrique ; réparation H beaucoup plus petite que T ; masse couplée ; champ quasi nul ; proximité d’un pôle du Schur conservé ; combinaison de charges annulantes. Aucune accélération ne doit provenir de la suppression d’une de ces difficultés.

Priorité : mesurer la branche Gram pour la grande norme et les facteurs de H, coût de préparation inclus, avec les champs et décisions historiques. Ensuite seulement, comparer la compression complète des métriques, en incluant les défauts ci-dessus. Ni l’interpolation empirique, ni les projections aléatoires, ni la précision accrue ne garantissent à elles seules le gain algorithmique recherché. Aucun des articles ne démontre que Vinkulum dépasserait un solveur concurrent sur cette campagne.
