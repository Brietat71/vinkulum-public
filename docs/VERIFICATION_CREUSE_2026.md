# Vérification creuse : examen mathématique des manuscrits de Rump (2026)

La vérification par résidus de factorisation ouvre une piste pour réduire le
coût des certificats de Vinkulum. Elle ne démontre actuellement aucun gain
sur notre certificat par trace inverse sélectionnée et arithmétique dirigée.
Une borne générale du premier manuscrit, prise telle qu'imprimée, admet
le contre-exemple exact ci-dessous. Cette observation ne réfute pas à elle
seule les algorithmes complets, leur implémentation ou la branche SPD.

## Sources et périmètre

La [page officielle TUHH](https://www.tuhh.de/ti3/publications.shtml?author=rump),
consultée le 8 septembre 2026, annonce deux manuscrits de S. M. Rump à paraître
dans SIMAX en 2026. Aucun DOI ni date de publication définitive n'a été vérifié.

| Manuscrit officiel | Pages | SHA-256 du PDF examiné |
| --- | ---: | --- |
| [Verified error bounds for sparse systems Part I: The splitting of a matrix into two factors](https://www.tuhh.de/ti3/paper/rump/sparselss_I_final.pdf) | 49 | `aaf94ce47e5023f6e22bdcd3f80c984669b69ceb1bd709e62fc276bb54df8a34` |
| [Verified error bounds for sparse systems Part II: Inertia-based bounds, least squares and nonlinear systems](https://www.tuhh.de/ti3/paper/rump/sparselss_II_final.pdf) | 48 | `2fc442b164b4c5e6fd68d2440c785ca95d370de75909520173824923b70132e6` |

La pagination citée est celle imprimée ; la page 9 correspond à l'indice PDF
8. Les pages 6, 9, 13 et 21 de Part I ont aussi été examinées visuellement.
Ce travail porte sur des énoncés mathématiques et nos calculs indépendants ;
aucune implémentation concurrente n'a été inspectée.

## Énoncé strictement réfuté

Le [lemme 2.3, équation (2.6), p. 9 de Part I](https://www.tuhh.de/ti3/paper/rump/sparselss_I_final.pdf#page=9)
énonce, pour des matrices binary64 compatibles et un produit arrondi au plus proche,

\[
\|\operatorname{fl}(AB)-AB\|_2
\le u\sum_k\min(\mu_k,\nu_k)\rho_k\sigma_k.
\]

Ici \(\mu_i,\rho_i\) sont le nombre de non-zéros et la norme euclidienne de
la ligne \(i\) de \(A\) ; \(\nu_j,\sigma_j\) concernent la colonne \(j\)
de \(B\). La norme matricielle est spectrale, explicitement utilisée p. 6.

**Contre-exemple et preuve indépendants.** Posons

\[
u=2^{-53},\quad q=2^{-27},\quad d=2^{-10},\qquad
A=\begin{pmatrix}1&q\\0&d\end{pmatrix},\quad
B=\begin{pmatrix}d&1\\0&q\end{pmatrix}.
\]

Tous les coefficients sont exactement représentables. Les déterminants
\(d\) et \(dq\) sont non nuls. Le produit exact est

\[
AB=\begin{pmatrix}d&1+q^2\\0&dq\end{pmatrix}.
\]

Un calcul séquentiel binary64, dans chacun des deux ordres de sommation,
arrondit \(1+2^{-54}\) à \(1\). Toutes les autres entrées sont exactes.
L'erreur ne possède qu'une entrée non nulle ; sa norme spectrale vaut
donc exactement \(q^2=2^{-54}\).

Or

\[
\mu=(2,1),\quad\nu=(1,2),\quad
\rho=(\sqrt{1+q^2},d),\quad\sigma=(d,\sqrt{1+q^2}).
\]

Le membre droit annoncé est \(2ud\sqrt{1+q^2}\). Sans évaluer de racine
flottante, \(1+q^2<4\) démontre

\[
2ud\sqrt{1+q^2}<4ud=2^{-61}<2^{-54}.
\]

L'inégalité imprimée est donc fausse pour ces facteurs inversibles.
La preuve utilise l'égalité incorrecte
\(\|xy^T\|_2=y^Tx\) ; l'identité valable est
\(\|xy^T\|_2=\|x\|_2\|y\|_2\).

## Témoin général conservatif, dérivé indépendamment

Considérons des produits puis sommes séquentiels, sans FMA, en arrondi
au plus proche. Les données et résultats exacts intermédiaires non nuls
restent dans la plage normale finie ; chaque opération satisfait ainsi
le modèle relatif \(\operatorname{fl}(z)=z(1+\delta)\), \(|\delta|\le u\).
Les zéros exacts sont permis.

Pour l'entrée \((i,j)\), soit \(k\) le nombre de produits non nuls.
Chaque terme subit au plus \(2k\) facteurs d'arrondi. Si \(2ku<1\),

\[
|E_{ij}|\le W_{ij}
:=\gamma_{2k}\sum_\ell |a_{i\ell}b_{\ell j}|,
\qquad \gamma_m=\frac{mu}{1-mu}.
\]

En effet, pour \(m\) facteurs, leurs écarts au produit unité sont majorés
par \(\gamma_m\) :
\((1-u)^m\ge1-mu\) et \((1+u)^m\le(1-mu)^{-1}\).
La majoration entrée par entrée suit par inégalité triangulaire.
Pour \(k=0\), l'erreur et la borne sont nulles. Enfin,

\[
\|E\|_2^2\le
\min\left(\sum_{ij}W_{ij}^2,\ \|W\|_1\|W\|_\infty\right).
\]

Cette borne n'utilise pas les termes diagonaux d'un majorant non symétrique
pour majorer sa norme. Le test calcule tout son membre droit en `Fraction`,
avec le dénominateur de \(\gamma_m\), et refuse les opérations hors hypothèses.
Il vérifie aussi des annulations et des ordres de sommation différents.
Ce témoin pédagogique n'est pas une implémentation de certification rapide.

## Ce qui reste exploitable et ce qui doit être mesuré

Dans Part I, le lemme 2.10 et le corollaire 2.12, pp. 12–14, combinent
un majorant de résidu Cholesky, Perron–Frobenius et des sommes préfixes.
Le théorème 6.1, p. 21, certifie une borne spectrale par Cholesky décalé.
Ces résultats ne supposent aucune dominance diagonale. La contre-épreuve
précédente ne les réfute pas.

Le [Part II, pp. 15–17, équations (5.1–5.2)](https://www.tuhh.de/ti3/paper/rump/sparselss_II_final.pdf#page=16)
utilise deux décalages opposés, des inerties concordantes et des résidus
majorés pour exclure zéro du spectre. Il signale le coût du remplissage
supplémentaire des facteurs décalés.

**Raccord propre à Vinkulum.** Pour \(M\succ0\) et une contrainte
\(B\in\mathbb R^{n\times s}\) de rang \(s\), posons

\[
C_\gamma=\begin{pmatrix}D^TD-\gamma M&B\\B^T&0\end{pmatrix}.
\]

Si \(Z\) engendre \(\ker B^T\), la congruence donne, dans l'ordre
(positif, négatif, nul),

\[
\operatorname{Inertia}(C_\gamma)
=\operatorname{Inertia}\bigl(Z^T(D^TD-\gamma M)Z\bigr)+(s,s,0).
\]

Certifier l'inertie \((n,s,0)\) démontre donc que la valeur propre
généralisée minimale sur le complément dépasse \(\gamma\), sans
construire \(Z\), une inverse dense ou une différence de grandes traces.

Les résidus devront porter sur \(D^TD\) exact, pas seulement son assemblage
flottant. Le facteur actuel issu de QR n'est pas automatiquement admissible
dans un théorème sur un facteur produit par Cholesky. Les contraintes
occupent \(O(ns)\) coefficients ; pivotage, remplissage et factorisations
supplémentaires peuvent absorber le gain. Les sommes préfixes et produits
creux offrent du parallélisme, sans supprimer les dépendances des facteurs.

Conserver le certificat dirigé actuel et mesurer séparément ce candidat :
temps complet, mémoire, remplissage et seuil certifié, à données identiques.
Les [preuves et dix tests rationnels d'inertie](INERTIE_COMPLEMENT_PREUVES.md)
éprouvent notamment \(K=\operatorname{diag}(10^{-30},1,2,3)\),
\(M=I\), \(B=e_1\) : acceptation de \(\gamma=1/2\), refus de
\(\gamma=3/2\), singularité à \(\gamma=1\). Ils n'établissent pas encore
un certificat machine creux à grande échelle.

Les contre-épreuves du produit matriciel s'exécutent depuis la racine :

```sh
python -m unittest discover -s ci -p test_bornes_produits_creux.py -v
```
