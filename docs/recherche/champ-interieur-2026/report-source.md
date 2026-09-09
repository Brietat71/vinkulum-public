# Vinkulum : des interfaces au champ physique

## Recherche et décisions — 7 septembre 2026

**Question : comment traduire davantage de raisonnement mathématique en une précision physique réellement contrôlée ?**

Le travail établit un critère d'erreur sur les déplacements, les déformations et les observables, puis l'utilise pour enrichir les intérieurs sans reconstruire leurs facteurs statiques. Il amorce aussi l'accès aux facteurs d'énergie directement dans le noyau Rust.

### Le résultat qui change le critère d'arrêt

L'erreur d'un transfert énergétique aux interfaces est quadratique en résidu. L'erreur du champ intérieur est linéaire. Une très petite erreur du premier ne suffit donc pas à accepter la seconde. Cette différence découle d'une identité exacte de condensation, et non d'une observation empirique isolée.

Sur une chaîne discrète de 2 048 inconnues, découpée en huit sous-structures, le contrôle physique porte la dimension intérieure réduite de 32 à 64. L'erreur relative de déformation passe de 3,70 × 10⁻⁸ à 2,89 × 10⁻¹³, en maximum sur les trois fréquences demandées. Le majorant final vaut 2,81 × 10⁻¹¹ et accepte la cible 10⁻¹⁰.

Deux expériences proches d'une résonance sont conservées comme refus. L'une a une solution assez précise mais une borne trop pessimiste ; l'autre dépasse effectivement la tolérance. Enrichir un intérieur ne corrige pas toute erreur d'assemblage amplifiée par une résonance globale.

### Portée exacte

Ces résultats concernent des modèles linéaires conservatifs et des constantes spectrales fournies. Les bornes locales couvrent une bande ; l'acceptation relative globale concerne seulement les fréquences demandées. Les opérations sont évaluées en doubles, sans certification des arrondis. Aucune nouvelle supériorité chronométrée sur Exudyn, MBDyn ou Simpack n'est établie.

---

## 1. Du résidu au champ : preuve locale

Posons K = DᵀD et A(ω) = K − ω²M. L'intérieur I et l'interface S partitionnent les coordonnées physiques. On suppose M_II positive et K_II ≥ λ₀ M_II, avec 0 ≤ ω ≤ Ω et Ω² < λ₀.

Pour un port imposé, le candidat û et la solution u ont la même trace. Si R = [Aû]_I et e = û − u, alors e_S = 0 et A_II e_I = R. Les blocs M_IS peuvent être non nuls : ils figurent dans R et disparaissent seulement de la différence à trace nulle.

En normalisant par M_II, l'opérateur devient H − ω²I avec H = M_II⁻¹ᐟ² K_II M_II⁻¹ᐟ² et spectre contenu dans [λ₀, +∞). Posons η = ‖R‖_(M_II⁻¹) et α = λ₀ − ω².

**Déplacements : ‖e‖_M ≤ η / α.**

La norme de (H − ω²I)⁻¹ est au plus 1/α. Cela donne directement la première inégalité, sans hypothèse d'orthogonalité de Galerkin du candidat.

**Déformations : ‖De‖ ≤ √λ₀ · η / α.**

La fonction t/(t − ω²)² décroît pour t > ω². Sa valeur maximale sur le spectre est donc λ₀/α². L'inégalité résulte du calcul fonctionnel spectral appliqué à la norme énergétique de l'erreur.

### Pourquoi le Schur peut sembler excellent

Pour la matrice candidate de relèvement X et le Schur exact S, l'identité d'énergie donne XᵀAX − S = RᵀA_II⁻¹R. Le terme de droite est quadratique en R. Dans le cas scalaire A_II = 1, une erreur intérieure ε produit une erreur énergétique ε². Il n'existe donc pas de majoration universelle linéaire du champ par cette seule erreur énergétique.

Dans le code, le Schur rapporté peut différer de XᵀAX en flottants. L'audit précédent conserve déjà cet écart séparément. La nouvelle borne du champ exploite le résidu dans D d'entrée ; elle ne confond pas ces deux erreurs.

---

## 2. Assembler les garanties, contrôler une observable

Chaque sous-structure possède un majorant local c_s du champ à trace imposée. Son application E_s transforme les ports globaux y en coordonnées locales normalisées. Les intérieurs sont disjoints ; leurs erreurs à trace nulle se somment quadratiquement. Les masses d'interface et les masses externes sont comptées dans la norme du relèvement. Les masses locales complètes et les énergies externes sont supposées positives semi-définies.

Soit B un majorant en norme de l'erreur du Schur assemblé Ŝ, σ sa plus petite valeur singulière et r_y = f − Ŝŷ le résidu de sa résolution. Si σ > B :

**b_y = (‖r_y‖ + B‖ŷ‖)/(σ − B)** majore l'erreur des ports.

Définissons η_loc² = Σ c_s²‖E_s ŷ‖², η_op² = λ_max(Σ c_s²E_sᵀE_s) et C la norme du relèvement candidat dans la métrique physique assemblée. Une décomposition algébrique du champ exact et candidat donne :

**‖u − û‖ ≤ η_loc + (C + η_op)b_y.**

La formule s'applique séparément à la norme massique et à la norme de déformation, avec leurs constantes respectives. Si σ ≤ B, le prototype refuse de conclure. La coercivité des intérieurs ne remplace pas cette stabilité globale.

### Un critère relatif calculable

Si B_u majore l'erreur absolue et N = ‖û‖ > B_u, alors ‖u − û‖/‖u‖ ≤ B_u/(N − B_u). Le dénominateur tient compte de l'écart entre la norme inconnue et celle du candidat. Si N ≤ B_u, aucun succès relatif n'est annoncé.

### Observable linéaire et problème dual

Pour J(u_I) = ℓᵀu_I, choisir un dual approché ẑ et son résidu r_d = ℓ − A_II ẑ. La correction du candidat est J_corr = J(û_I) − ẑᵀR. L'identité exacte du reste donne :

**|J(u_I) − J_corr| ≤ ‖r_d‖_(M_II⁻¹) · η / α.**

Un dual exact annule le reste en arithmétique exacte, sans exiger l'enrichissement du champ primal. Cela ouvre un contrôle ciblé des grandeurs utiles ; le coût de construction du dual doit rester compté. L'antériorité est celle des méthodes primal-dual, explicitée par [Rannacher, 2003](https://arxiv.org/pdf/math/0305006).

---

## 3. Cinq expériences, y compris les refus

Chaîne discrète : n masses et n ressorts unitaires, racine fixée, force physique terminale 1/√n. Les sous-structures sont construites indépendamment. L'oracle est la solution analytique du système discret complet ; il n'intervient pas dans l'arrêt.

Bande locale : Ω = 0,8 · 2 sin(π/(2n)). Tolérance relative : 10⁻¹⁰ dans les deux normes physiques. Les trois premiers essais demandent les fréquences 0, 0,37Ω et Ω. Les chiffres sont des maxima sur ces points.

| n / segments | Directions avant → après | Erreur déformation avant | Erreur après | Majorant après |
|---|---|---|---|---|
| 64 / 4 | 24 → 32 | 5,23 × 10⁻¹⁰ | 1,26 × 10⁻¹³ | 1,68 × 10⁻¹² |
| 2048 / 8 | 32 → 64 | 3,70 × 10⁻⁸ | 2,89 × 10⁻¹³ | 2,81 × 10⁻¹¹ |
| 2048 / 32 | 128 → 192 | 3,58 × 10⁻¹¹ | 2,81 × 10⁻¹³ | 7,42 × 10⁻¹¹ |

Les trois critères passent. Ajouter respectivement 4, 8 et 32 ports globaux aux dimensions intérieures indiquées. Les facteurs sont réutilisés ; les directions précédentes sont conservées. Chaque ronde enrichit tous les intérieurs, sans prétention de sélection optimale.

Les deux essais suivants demandent 17 fréquences uniformes sur [0, Ω], dont une proche du premier mode global. Après trois rondes, leurs critères restent refusés :

| n / segments | Directions finales | Erreur déformation | Majorant | Lecture |
|---|---|---|---|---|
| 2048 / 8 | 80 | 2,55 × 10⁻¹² | 3,42 × 10⁻⁹ | Borne conservatrice |
| 2048 / 32 | 320 | 3,55 × 10⁻¹⁰ | 9,25 × 10⁻⁸ | Erreur trop grande |

L'élargissement de l'ensemble de fréquences révèle donc une limite qui serait invisible dans un bilan fondé uniquement sur les trois points précédents. Aucune acceptation entre les points testés n'est déduite de la grille.

Onze tests indépendants couvrent aussi les masses couplées, l'assemblage dense complet de petits systèmes, les identités primal-dual, les budgets et les refus. Sources, données et manifeste sont conservés dans docs/bancs/champ-interieur-2026. L'expérience ne mesure pas les performances temporelles.

---

## 4. Passer au noyau sans perdre l'énergie physique

Les poutres du noyau disposent déjà de déformations γ, κ et de rigidités effectives positives. Le nouvel accès Rust calcule leurs six déformations pondérées z et le jacobien matériel D par dérivation première. Il ne reconstruit pas D en factorisant une matrice K déjà arrondie.

**U = ½‖z‖² et f = −Dᵀz.**

Dans des coordonnées différentiables, la dérivée seconde de l'énergie comprend DᵀD et les termes géométriques Σ z_j Hess(z_j). Au repos sans contraintes internes, le second terme s'annule. Sous précontrainte, il doit être conservé. Les conventions de rotation et de transport des efforts doivent également être respectées.

Cette distinction est testable et évite de transformer arbitrairement une tangente indéfinie en énergie positive. L'accès ajouté reste interne à Rust : les ports Python et les trajectoires réduites n'en sont pas encore consommateurs. Le superélément corotationnel existant est conservé ; son import de K ne fournit pas un facteur élémentaire d'origine.

### Une condition suffisante pour une extension ultérieure

Pour une tangente conservative symétrique K_tan = DᵀD + G, si D_IᵀD_I ≥ λ_D M_II et G_II ≥ −γ M_II, alors K_tan,II − ω²M_II ≥ (λ_D − γ − ω²)M_II. Cette implication algébrique fournit le contrat de stabilité à établir. Elle ne donne pas automatiquement γ, et le résidu doit inclure G.

Le cas non linéaire demande aussi un voisinage, une constante de stabilité et un contrôle de variation de la dérivée. [Smetana–Taddei, 2023](https://doi.org/10.1137/22M148402X) illustre précisément ces exigences dans un cadre elliptique ; elles ne sont pas acquises pour un multicorps avec contact et événements.

### Décision scientifique

La prochaine difficulté identifiée est le défaut d'assemblage amplifié par les résonances, puis la stabilité des tangentes précontraintes. La priorité consiste à mieux contrôler ces termes. Les preuves actuelles justifient un enrichissement physique et son refus éventuel ; elles ne justifient pas une affirmation de rupture de performance globale.

---

## Sources et accès

1. Rolf Rannacher, ICM 2002, prépublication 2003. Adaptive Finite Element Methods for Partial Differential Equations. Propositions 1–2 : représentation primal-dual avec reste cubique, nul dans le cas linéaire. [Texte primaire consulté](https://arxiv.org/pdf/math/0305006).

2. Kathrin Smetana et Anthony T. Patera, 2016. Optimal Local Approximation Spaces for Component-Based Static Condensation Procedures. Proposition 3.5 : largeur de Kolmogorov égale à √λ_(n+1) pour l'opérateur de transfert compact et les métriques prescrites. Le spectral greedy vise un ensemble fini d'apprentissage. [Article primaire consulté](https://ris.utwente.nl/ws/portalfiles/portal/168414349/15m1009603.pdf).

3. Andreas Buhr et Kathrin Smetana, 2018. Randomized Local Model Order Reduction. Propositions 3.7–3.8 : majorants probabilistes d'erreur de transfert, avec constantes de métrique et budget d'échec adaptatif. Cette méthode aléatoire n'est pas implémentée dans le prototype présenté. [Article primaire consulté](https://arxiv.org/pdf/1706.09179).

4. Andreas Buhr, 2018. Exponential Convergence of Online Enrichment in Localized Reduced Basis Methods. Théorèmes 1 et 4 : contraction et optimalité dans des espaces locaux recouvrants, sous coercivité et partition de l'unité stable. La contraction ne s'étend pas directement à notre dynamique globalement indéfinie. [Texte primaire consulté](https://arxiv.org/pdf/1710.02104).

5. Kathrin Smetana et Tommaso Taddei, 2023. Localized Model Reduction for Nonlinear Elliptic Partial Differential Equations, SISC 45(3), A1300–A1331. Le préprint établit une borne locale-globale et un contrôle de Brezzi–Rappaz–Raviart sous constantes quantitatives. [Publication](https://doi.org/10.1137/22M148402X), [préprint consulté](https://www.math.u-bordeaux.fr/~ttaddei/data/KSTT_arxiv.pdf).

Lecture ciblée arrêtée le 7 septembre 2026 après couverture des cinq familles d'énoncés et vérification indépendante des passages décisifs. Aucun code concurrent consulté. Les identités adaptées sont détaillées dans docs/CHAMP_INTERIEUR_PREUVES.md ; les résultats et commandes dans docs/CHAMP_INTERIEUR_PROTOTYPE.md ; le contrat Rust dans docs/FACTEURS_ENERGIE_NOYAU.md.

Le contrôle primal-dual, les espaces optimaux de transfert et l'enrichissement résiduel sont des antériorités établies. L'apport étudié ici est leur composition avec les facteurs énergétiques et l'assemblage de Vinkulum. Aucune nouveauté mondiale n'est revendiquée.
