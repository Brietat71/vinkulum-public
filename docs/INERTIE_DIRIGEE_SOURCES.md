# Inertie creuse vérifiée : pistes ciblées pour Vinkulum

Recherche arrêtée le 8 septembre 2026. Cinq références primaires ci-dessous ; aucun code concurrent consulté, aucun solveur relancé. Les pages indiquées sont imprimées, donc indice PDF = page − 1, sauf pagination S439–S442. Les textes des théorèmes récents ont été vérifiés dans les PDF/HTML ; leurs captures visuelles n'ont pas été disponibles. Le scan de 1995 ne permettait pas de relire son théorème : l'attribution précise vient de sa reprise explicite en 2026. Aucun gain de performance de ces méthodes n'est démontré sur Vinkulum.

## Références et hypothèses utiles

Cette lecture bibliographique est distincte de la
[campagne du prototype par congruences dirigées](INERTIE_CONTRAINTE_PROTOTYPE.md).
Cette dernière fournit des mesures propres à Vinkulum ; elle ne mesure
pas les implémentations des articles ni le repli unilatéral proposé ci-dessous.

1. **S. M. Rump, _Verified Computation of the Solution of Large Sparse Linear Systems_, ZAMM 75, S439–S442, 1995.** [Notice institutionnelle](https://www.tuhh.de/ti3/publications.shtml?author=rump), [PDF original](https://www.tuhh.de/ti3/paper/rump/Ru95a.pdf). Antériorité exacte du procédé d'inertie par décalages et résidus, explicitement identifié comme théorème 1.1 de cette référence dans Part II (2026), p. 2. La nouveauté recherchée ne peut donc pas être simplement « employer LDL et l'inertie pour vérifier un système creux ».

2. **S. M. Rump, _Verified error bounds for sparse systems Part II: Inertia-based bounds, least squares and nonlinear systems_, 2026, 48 pages, annoncé à paraître dans SIMAX.** [PDF](https://www.tuhh.de/ti3/paper/rump/sparselss_II_final.pdf), [statut officiel](https://www.tuhh.de/ti3/publications.shtml?author=rump). Aucun DOI définitif vérifié. Théorème 1.1, p. 2 : deux facteurs des décalages opposés, inerties concordantes et majorants de résidus inférieurs au décalage excluent zéro. Imposer explicitement des facteurs de congruence inversibles. Section 4, pp. 13–14 : signes des blocs 2×2 par déterminant et trace ; les transformations sans erreur traitent la cancellation de deux produits. Cela certifie le bloc flottant donné, pas automatiquement le système original. Section 5, pp. 15–17 : décalages dirigés, factorisation en deux facteurs pour réduire chaque résidu à un produit scalaire, précision renforcée et remplissage supplémentaire. Hypothèses arithmétiques pp. 3–5 : modes dirigés réellement respectés, maîtrise des dépassements/sous-flux. Un nombre de chiffres annoncé par une bibliothèque ne constitue pas ce contrat. SHA-256 du PDF local : `2fc442b164b4c5e6fd68d2440c785ca95d370de75909520173824923b70132e6`.

3. **T. Terao, Y. Watanabe, K. Ozaki, _Verified error bounds for the singular values of structured matrices with applications to computer-assisted proofs for differential equations_, prépublication arXiv v1 déposée le 14 février 2025.** [Notice et version](https://arxiv.org/abs/2502.09984), [PDF](https://arxiv.org/pdf/2502.09984v1#page=13). Le manuscrit porte la date du 17 février ; publication évaluée non vérifiée. Pour leur matrice de masse \(W=R^HR\succ0\), le théorème 3, p. 13, transforme la borne sur \(R^{-H}AR^{-1}\) en comptage des positives de \(G_\theta=[[\theta W,A^H],[A,\theta W]]\). Le théorème 4, p. 14, utilise **un seul décalage**,
   \(\|G_\theta+\tau I-\widehat L\widehat D\widehat L^H\|_2\le\delta<\tau\), avec exactement \(n\) positives de \(\widehat D\) et une borne structurelle opposée sur le comptage. L'inverse de \(R\) n'est pas formée. Le théorème 2, p. 12, propose aussi deux tests SPD sur des matrices de Gram ; les auteurs signalent remplissage et conditionnement aggravé. La matrice augmentée double la dimension : aucune raison de l'ajouter à notre KKT déjà symétrique.

4. **M. Lange, S. M. Rump, _Accurate floating-point matrix residuals_, 2026, 50 pages, soumis pour publication.** [Statut institutionnel](https://www.tuhh.de/ti3/publications.shtml?author=rump), [PDF](https://www.tuhh.de/ti3/paper/rump/MatrixResidualsFinal.pdf#page=44). Le lemme 4, p. 7, fournit une condition d'exactitude pour le produit scalaire des parties significatives : coefficients multiples de puissances de la base, produit des échelles représentable, taille des sommes intermédiaires bornée ; aucun dépassement. La section 8.4, pp. 43–45, adapte les décompositions aux non-zéros stockés. Avec \(n_{ij}\) nombre de produits non nuls, équations (40–42), le choix de l'échelle commune satisfait \(k\max_{ij}n_{ij}/(4\tau^2)\le\beta^p\), permettant d'additionner exactement les contributions de même importance. C'est une piste pour **resserrer un résidu**, sans augmenter la précision de toute la factorisation. Les commentaires de coût presque indépendant du nombre de parties supposent une répartition favorable des chiffres/non-zéros ; ce n'est pas une complexité universelle. Le bénéfice BLAS dense diminue dans le creux. Ne pas transposer leurs chronométrages à Python/Decimal.

5. **R. Parker, M. Garcia, R. Bent, _Exploiting block triangular submatrices in KKT systems_, prépublication arXiv v1, 20 février 2026.** [Notice](https://arxiv.org/abs/2602.17968), [texte primaire, sections 3.1–3.4](https://arxiv.org/html/2602.17968v1). Le théorème 3 donne l'inertie \((m,m,0)\) de \([[H,J^T],[J,0]]\) lorsque \(J\) est **carrée inversible**. Les auteurs l'identifient comme cas particulier d'un résultat de Gould (1985). Leur accélération exploite en plus une Jacobienne triangulaire par blocs ; elle ne constitue pas une certification des arrondis. Le complément de Schur remplit les colonnes de liaison. Notre covecteur spectral, long et dense, ne fournit pas cette structure à lui seul. Ne pas confondre permutations différentes des lignes/colonnes, utiles pour résoudre, et congruence symétrique, nécessaire pour transporter l'inertie.

## Conclusion mathématique propre à notre KKT

Notons \(D_e\) l'opérateur d'énergie, pour éviter de le confondre avec les pivots LDL :

\[
 C_\gamma=\begin{pmatrix}D_e^TD_e-\gamma M&B\\B^T&0\end{pmatrix},
 \qquad B\in\mathbb R^{n\times s},\quad\operatorname{rang}B=s.
\]

Le produit \(D_e^TD_e\) désigne le produit **exact des données d'entrée**. Pour \(Z\) engendrant \(\ker B^T\), la congruence donne

\[
\operatorname{inertie}(C_\gamma)
=\operatorname{inertie}(Z^T(D_e^TD_e-\gamma M)Z)+(s,s,0).
\]

Cette identité, déjà documentée dans Vinkulum, n'impose aucune dominance diagonale. Pour en déduire une valeur propre généralisée, la masse doit être positive sur le sous-espace considéré. La masse native diagonale strictement positive satisfait cette condition.

**Le LDL par intervalles est une première option correcte**, si chaque Schur encadre le Schur exact issu du même système et chaque pivot est certifié inversible. Pour un bloc symétrique \([[a,b],[b,c]]\), un intervalle strictement négatif pour \(ac-b^2\) donne une positive et une négative ; un déterminant strictement positif avec signe certifié de la trace donne deux signes identiques. Si le déterminant contient zéro : raffiner/changer de pivot/refuser. Le centre de l'intervalle ne tranche jamais. Supprimer un remplissage jugé petit ou régulariser un pivot changerait la matrice à certifier.

**Repli proposé : un seul décalage négatif, preuve indépendante.** Prenons \(\tau>0\), \(H=\widehat L\widehat D\widehat L^T\), avec \(\widehat L\) inversible, au moins \(n\) valeurs propres positives de \(H\), et

\[
 \|C_\gamma-\tau I-H\|_2\le\delta<\tau.
\]

Par Weyl, les \(n\) plus grandes valeurs propres de \(C_\gamma\) sont strictement positives. La structure et le rang de \(B\) assurent déjà au moins \(s\) négatives. La dimension \(n+s\) impose donc exactement \((n,s,0)\). Nous adaptons ici l'idée unilatérale du théorème 4 de Terao et al. ; ce n'est pas son énoncé original. Pour \(s=1\), une entrée exacte non nulle de \(B\) prouve son rang. Pour plusieurs colonnes, le rang doit aussi être certifié. Sans cette hypothèse, l'argument unilatéral n'est pas justifié.

Le résidu doit inclure formation de Gram, décalage et produits des facteurs. Soit E le résidu exact, et W un majorant entrée par entrée de sa valeur absolue. Alors \(\|E\|_2\le\sqrt{\|W\|_1\|W\|_\infty}\). Choisir δ égal à ce dernier majorant arrondi vers le haut. Si W est symétrique, la somme de ligne maximale suffit. Tous ces calculs doivent arrondir vers l’extérieur. On n’utilise pas le lemme 2.3 général de Part I réfuté dans docs/VERIFICATION_CREUSE_2026.md ; cela ne réfute pas sa branche SPD/Perron.

## Si les intervalles gonflent

1. Équilibrer par congruence diagonale de puissances de deux, en contrôlant la représentabilité. Cela conserve l'inertie ; les unités du bloc de contraintes peuvent sinon rendre les résidus normiques inutilisables.
2. Séparer les refus : matrice effectivement proche du seuil, croissance des pivots, largeur introduite par Gram, largeur propagée des Schur. Augmenter seulement la précision d'une dernière division ne répare pas des entrées déjà trop larges : reconstruire les fronts dépendants depuis des données certifiées ou recommencer à précision supérieure, et compter tout le coût.
3. Si la factorisation ponctuelle paraît stable mais la récurrence d'intervalles échoue, vérifier son résidu avec des produits scalaires compensés/expansions et le décalage unilatéral ci-dessus. La séparation stricte reste obligatoire. Une petite erreur rétrograde relative, seule, ne certifie aucune inertie.

Contre-épreuves minimales : \([[0,1],[1,0]]\), déterminant 2×2 presque nul par cancellation, permutation physique des coordonnées, longue chaîne de Schur, seuil exactement égal à une valeur propre contrainte, et contrainte de rang déficient. Mesurer temps complet, mémoire/remplissage et refus, à seuil et entrées identiques. Conserver le certificat dirigé existant comme témoin ; aucune reprise, diminution de seuil ou régularisation silencieuse.
