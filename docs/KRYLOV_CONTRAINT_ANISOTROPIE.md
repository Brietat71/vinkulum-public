# Revue du contrôle et variante anisotrope

Note du 8 septembre 2026. Les formules visent l'arithmétique exacte ;
leurs évaluations binary64 conservent le statut non certifié des arrondis.
La [preuve de l'enveloppe](KRYLOV_CONTRAINT_PREUVES.md) donne les matrices
résiduelles et les conventions partagées avec le contrôleur.

## 1. Contrat du contrôle

La norme massique duale du complément utilise une racine supérieure
intérieure R telle que M_ii=RᵀR. Les réflecteurs de R⁻ᵀB retirent les
réactions avant toute norme. Le modèle du Schur et les normes des champs
demandés emploient M original et D original.

La réparation vise le right-inverse théorique C*=C0(BᵀC0)⁻¹ et le champ
contraint associé. La normalisation finale des directions peut amplifier
BᵀV : le prototype répare ce défaut après normalisation puis le conserve
explicitement dans son audit. Il ajoute aussi le défaut du petit solve
et l'écart de fonctionnel énergétique.

Les calculs de rang, de normes et de résidus restent flottants. Ils ne
prouvent pas les identités machine exactes. L'enrichissement invalide le
contrôleur préparé ; modifier arbitrairement ses tableaux en place n'est
pas pris en charge.

La racine intérieure est limitée aux composantes connexes SPD de taille
au plus six. Certaines masses consistantes plus largement couplées restent
hors de ce prototype. La masse complète peut être PSD, notamment avec des
ports sans masse : aucune racine SPD globale n'est exigée. Le cas où aucun
vecteur de Krylov n'est nécessaire est également traité.

## 2. Enveloppes par colonne

Noter R(mu) le résidu du champ réparé idéal, et H_M un facteur tel que
S_M=H_M^T H_M soit la métrique duale contrainte. Le nombre de colonnes
m=p+s comprend tous les ports et toutes les coordonnées retenues.

Pour chaque colonne j, construire
\[
\delta_j\ge\sup_{\mu\in[0,\rho]}\|H_MR(\mu)e_j\|_2.
\]
La même préparation Bernstein suffit. Si B_{q,k} sont les coefficients
matriciels Bernstein déjà construits, F_d=H_M F et t≥||Theta||,
\[
\boxed{\delta_j^{(q)}=
\max_k\|B_{q,k}e_j\|_2+
\frac{\rho^q\|F_d\Theta^{q-2}\|_2\,\|D e_j\|_2}{1-\rho t}.}
\]
On peut prendre le minimum sur les profondeurs q pour chaque j séparément.
Chaque borne reste uniforme : les profondeurs choisies ne doivent pas
nécessairement coïncider entre colonnes.

Conserver aussi l'enveloppe d'opérateur actuelle delta_op. Les deux
familles donnent
\[
\delta_*=\min(\delta_{\rm op},\|\delta\|_2),\qquad
\boxed{\eta(q)=\min\big(\delta_*\|q\|_2,\ \sum_j\delta_j|q_j|\big).}
\]
Pour tout mu et tout q, ||H_M R(mu)q||≤eta(q). Ce résultat est simultané
pour tous les q ; il reste donc valable lorsque q est lui-même la réponse
calculée à cette fréquence.

La seconde majoration suit de la somme des colonnes ; la première de la
norme d'opérateur. La relation ||R||≤||delta||₂ découle de Cauchy–Schwarz.
Les calculs supplémentaires sont des normes de colonnes des coefficients
déjà préparés et de D. Aucune factorisation supplémentaire n'est nécessaire.
Le bénéfice vient des directions peu affectées par le résidu. Réduire
seulement l'amplitude d'une charge multiplie toutes ces bornes absolues
par la même amplitude et n'améliore pas, à lui seul, leur qualité relative.

## 3. Action du défaut de Schur : conserver un facteur global

À la fréquence z, poser alpha=lambda_c−z>0 et
\[
A=\bar S-S=R^T R_zR\succeq0.
\]
La coercivité donne
\[
A\preceq R^T S_MR/\alpha,\quad
\|A\|\le\delta_*^2/\alpha,\quad
q^T A q\le\eta(q)^2/\alpha.
\]
Comme A est PSD,
\[
\|Aq\|^2=q^T A^2q\le\|A\|q^T A q.
\]
Donc
\[
\boxed{\|(\bar S-S)q\|_2\le\frac{\delta_*\,\eta(q)}{\alpha}.}
\]
Cette formule améliore delta_*²||q||/alpha tout en conservant le facteur
global nécessaire à l'action d'une matrice.

Une dérivation composante par composante donne aussi
\[
|(Aq)_i|\le\delta_i\,\eta(q)/\alpha,
\qquad
\|Aq\|\le\|\delta\|_2\,\eta(q)/\alpha.
\]
L'emploi de delta_* au lieu de ||delta||₂ peut être plus fin.

Attention : |A_ij|≤delta_i delta_j/alpha est une majoration ABSOLUE des
coefficients. Elle ne prouve pas A≼delta delta^T/alpha en ordre de Loewner.

## 4. Numérateur anisotrope avec marge globale

Le Schur réellement renvoyé Shat comprend encore les défauts de
réparation du champ non contraint, du petit solve et de l'évaluation de
l'énergie. Écrire
\[
\widehat S-S=A+E_f,\qquad \|E_f\|\le\epsilon_f.
\]
Une variante immédiatement compatible avec le contrôle actuel est
\[
b_{\rm global}=\delta_*^2/\alpha+\epsilon_f,\qquad
\mathfrak m=\sigma_{\min}(\widehat S)-b_{\rm global}.
\]
Si cette marge est strictement positive, pour qhat calculé et
r=f−Shat qhat,
\[
\boxed{
b_q=
\frac{\|r\|+
       \delta_*\eta(\widehat q)/\alpha+
       \epsilon_f\|\widehat q\|}
     {\mathfrak m}.}
\]
La marge ne dépend pas favorablement du seul qhat : elle doit protéger
l'inverse dans toutes les directions. Le numérateur peut exploiter le
vecteur réellement demandé. Garder alpha_* au lieu d'alpha dans tout
le contrôle reste valide mais plus pessimiste.

On peut conserver exactement la marge globale du code existant et
remplacer seulement son terme delta²||q||/alpha par la formule anisotrope.
Cela évite tout changement dans le raisonnement de stabilité.

## 5. Réparation et petit solve à la direction demandée

Pour le champ brut X et le champ réparé Xbar=X+E H Y, une borne locale est
\[
e_M(q)=\eta(q)/\alpha+\|H Yq\|_M,\qquad
e_D(q)=\sqrt{\lambda_c}\eta(q)/\alpha+\|D_i H Yq\|.
\]
Si l'on souhaite seulement des constantes uniformes par colonne, poser
\[
y_j=
\frac{\|a_0e_j\|+\rho\|a_1e_j\|}{1-\rho t},\quad
\upsilon(q)=\min(Y_{\max}\|q\|,\sum_j y_j|q_j|).
\]
Alors ||H Yq||_M≤||H||_M upsilon(q), et de même dans D.
Le facteur H a rang au plus s : garder C et B^T W séparés permet de
calculer son action sans dilater un opérateur de rang faible.

À une fréquence donnée, avec
tau=(a0+mu a1)−(I−mu Theta)Yhat,
\[
\epsilon_Y(q)=\|\tau q\|/(1-\mu t)
\]
majore ||(Yhat−Y)q||. Cela remplace le produit
||tau||||q||/(1−mu t), sans coût de grande dimension. Ajouter
||M_{ii}^{1/2}W|| epsilon_Y(q) et ||D_iW|| epsilon_Y(q) aux erreurs
de champ locales. Pour la réparation, on peut utiliser
\[
\|HYq\|_M\le\|HYhat\,q\|_M+\|H\|_M\epsilon_Y(q).
\]
Toutes ces normes doivent être évaluées après l'action sur q, ou au moyen
de facteurs de norme, sans soustraction de grands Grams.

La propagation de l'erreur de coordonnées reste globale :
\[
b_M=e_M(\widehat q)+(\|\widehat X\|_{M,2}+c_{M,\rm op})b_q,
\]
et pareil pour D. On ne peut pas remplacer c_M,op dans cette extension
par une constante favorable à la seule direction qhat : la direction de
q_exact−qhat est inconnue.

## 6. Raffiner aussi le défaut du fonctionnel, si nécessaire

Pour deux champs U et U+J, noter
u_D=||DU||, j_D=||DJ||,
u_D(q)=||DUq|| et j_D(q)=||DJq||, avec analogues massiques.
L'expansion
\[
(U+J)^TA(U+J)-U^TAU=U^TAJ+J^TAU+J^TAJ
\]
donne l'action majorée
\[
\boxed{\beta_f(q)=
(u_D+j_D)j_D(q)+j_Du_D(q)
+z[(u_M+j_M)j_M(q)+j_Mu_M(q)].}
\]
Des majorants des quatre normes directionnelles conviennent.
La borne d'opérateur globale reste
2u_Dj_D+j_D²+z(2u_Mj_M+j_M²).

Cette formule s'applique séparément à J=EHY pour la réparation et à
J=E W(Y−Yhat) pour le petit solve, en choisissant U de façon cohérente.
On peut alors remplacer epsilon_f||qhat|| dans le numérateur par la somme
des beta_f(qhat), tout en conservant la même marge globale. Ce raffinement
n'est utile que si les défauts de fonctionnel dominent encore le résidu ;
il n'est pas nécessaire à la validité de la première variante.

## 7. Contre-exemples et limite de l'approche par colonnes

Prendre alpha=1, H_M R=(epsilon,1), avec 0<epsilon<1. Alors
\[
A=R^TR=
\begin{bmatrix}\epsilon^2&\epsilon\\\epsilon&1\end{bmatrix},
\quad q=e_1.
\]
La norme résiduelle directionnelle est eta(q)=epsilon, tandis que
\[
\|Aq\|=\epsilon\sqrt{1+\epsilon^2}.
\]
La borne isotrope vaut 1+epsilon² ; la borne anisotrope
delta_* eta(q)=epsilon sqrt(1+epsilon²) est exacte.
Pour epsilon=10^{-6}, le gain est voisin d'un million dans ce majorant.
Mais eta(q)²=epsilon² sous-estime l'action : le résidu faible de la
première colonne se couple à la seconde colonne forte. L'énergie
q^T A q=epsilon² ne peut pas être renommée norme de Aq.

Autre exemple : R=(1,−1), delta=(1,1).
Pour q=(1,−1), q^T R^T Rq=4 mais q^T delta delta^Tq=0 :
la prétendue majoration de Loewner par delta delta^T est fausse.
Pour q=(1,1), le résidu réel est nul alors que sum delta_j|q_j|=2 :
les seules normes de colonnes perdent les compensations entre charges.

Si cette dernière limite devient importante, une variante préserve les
matrices Bernstein :
\[
\eta_q(\widehat q)=
\max_k\|B_{q,k}\widehat q\|+
\frac{\rho^q\|F_d\Theta^{q-2}\|\,\|D\widehat q\|}{1-\rho t}.
\]
Elle reste une enveloppe uniforme pour tout q demandé, et permet le même
raisonnement d'action de Schur avec delta_* eta_q(q)/alpha.
Pour éviter de stocker tous les grands B_{q,k}, un QR mince peut conserver
un facteur triangulaire R_k avec ||B_{q,k}q||=||R_kq|| en exact.
Il faut appliquer ce facteur à q, sans reconstruire un Gram puis
évaluer q^T Gram q. Les arrondis des QR et petites combinaisons restent
à qualifier ; aucun certificat machine n'est acquis par cette réduction.

Le [prototype](../ci/controle_complement.py) applique l’enveloppe par
colonnes au numérateur et garde une marge globale sur les coordonnées
conservées. Les arrondis restent évalués sans certificat machine.
