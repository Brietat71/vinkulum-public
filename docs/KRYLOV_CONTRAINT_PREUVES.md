# Enveloppe uniforme du Krylov sur le complément contraint

Note du 8 septembre 2026, issue de la relecture des contrats 0.11.0 et de
[RETENTION_INTERIEURE_PREUVES.md](RETENTION_INTERIEURE_PREUVES.md). Les identités sont en arithmétique exacte.
Des évaluations binary64 de ces majorants restent non certifiées tant que
leurs arrondis ne sont pas encadrés. Aucun gain de temps n'est déduit ici.

## 1. Données et inconnues conservées

Noter K=D_i^T D_i, M=M_ii, N=ker(B^T), avec K,M SPD. Le B effectivement
stocké définit exactement N. On dispose d'un minorant lambda_c tel que
v^T K v≥lambda_c v^T Mv pour v∈N et Omega²<lambda_c.

Le relèvement complet T=[L,E Phi] porte m=p+s coordonnées conservées :
ports normalisés et s coordonnées intérieures retenues. E injecte les
coordonnées intérieures. Avec la raideur complète mathématique D^T D et
la masse complète M_full, poser
\[
G(z)=G_0-zG_1,\quad
G_0=D_i^T(DT),\quad G_1=(M_{\rm full}T)_i.
\]
Toutes les colonnes de T interviennent, ainsi que tous les couplages K/M.
En particulier, B=fl(M Phi) ne permet pas d'annuler un couplage massique
exact. Un changement de normalisation de Phi n'autorise pas à changer B
sans mettre à jour le contrat de complément et son certificat.

## 2. Petits systèmes et signes du résidu

Partir d'une base V⊂N et la normaliser dans l'énergie du modèle D :
V^T K V=I. Pour une base brute V0 et V0^T K V0=LL^T, il s'agit de
V=V0 L^{-T}. En exact,
\[
\mu=z/\lambda_c,\quad \rho=\Omega^2/\lambda_c<1,\quad
\Theta=\lambda_c V^T M V,\quad
a_0=V^T G_0,\quad a_1=-\lambda_c V^T G_1.
\]
La coercivité restreinte implique 0≺Theta≼I. Les coefficients et le champ sont
\[
Y(\mu)=(I-\mu\Theta)^{-1}(a_0+\mu a_1),\qquad X=T-EVY.
\]
On diagonalise Theta ou factorise les petits systèmes ; aucune
factorisation KKT de taille n+s n'est nécessaire par fréquence.

Définir les matrices constantes du modèle physique d'entrée
\[
P=KV=D_i^T(D_iV),\quad Q=\lambda_c MV,
\]
\[
E_0=G_0-Pa_0,\quad E_1=-\lambda_c G_1-Pa_1,\quad
F=Q-P\Theta,
\]
\[
C_1=E_1+Fa_0,\qquad D_*=\Theta a_0+a_1.
\]
Le résidu intérieur du champ est exactement
\[
\boxed{R(\mu)=G_0-\mu\lambda_cG_1-(P-\mu Q)Y
=E_0+\mu E_1+\mu FY}
\]
puis, puisque Y=a0+mu(I−mu Theta)^{-1}D_*,
\[
\boxed{R(\mu)=E_0+\mu C_1+\mu^2F(I-\mu\Theta)^{-1}D_*.}
\]

La seconde identité reste vraie pour T,V,Theta,a0,a1 arbitraires fixés,
si l'on utilise bien les définitions P,Q,E0,E1 ci-dessus et la même
résolvante pour Y. Elle n'exige donc pas un Galerkin parfait, une
orthogonalité énergétique parfaite, ni une factorisation Q exacte.
En revanche, sans ces hypothèses, on ne peut pas conclure ||Theta||≤1.
Il faut alors disposer séparément d'un majorant t≥||Theta||₂ avec rho*t<1.
Une valeur singulière simplement calculée n'est pas un encadrement machine.

## 3. Norme duale restreinte sans différence de Grams

Pour une matrice de forces R, la norme pertinente est
\[
\|R\|_{S_M,2}=\|H_M R\|_2,\qquad
S_M=M^{-1}-M^{-1}B(B^TM^{-1}B)^{-1}B^TM^{-1}.
\]
Une force de réaction B eta est nulle dans cette norme. Avec M=L_M L_M^T,
écrire Q_B=orth(L_M^{-1}B). Alors
\[
H_M=(I-Q_BQ_B^T)L_M^{-1},\qquad S_M=H_M^TH_M.
\]
Un QR de Householder de L_M^{-1}B permet d'appliquer la transformation
orthogonale et de ne conserver que les n−s lignes restantes. Les
réflecteurs occupent O(ns) coefficients, sans base dense globale de N.

Il faut transformer les colonnes résiduelles avant leur Gram : ne pas
calculer R^T M^{-1}R moins le Gram de la réaction. Pour M=I, B=e1,
r=(1,2^{-30}), cette différence vaut exactement 2^{-60}, mais les deux
grandes quantités s'annulent en binary64. Les coordonnées restantes
conservent directement 2^{-30}.

Cela ne prouve pas une bonne erreur relative pour toute projection presque
nulle : les produits et transformations ont leurs arrondis, et une réaction
énorme peut encore masquer une petite composante. Des réflecteurs stables
évitent la différence de Grams ; une certification de la petite composante
demande toujours une borne d'erreur absolue des opérations ou plus de
précision. Même le résidu physique initial doit être évalué avec soin.

## 4. Enveloppe uniforme par résolvante et polynôme

Pour q≥2, poser A0=E0, A1=C1 et Aj=F Theta^{j−2} D_* pour 2≤j≤q−1.
Alors, exactement,
\[
R(\mu)=P_{q-1}(\mu)+
\mu^qF\Theta^{q-2}(I-\mu\Theta)^{-1}D_*,
\quad P_{q-1}(\mu)=\sum_{j=0}^{q-1}\mu^jA_j.
\]
La somme des Aj pour j≥2 est vide à q=2. Avec rho*t<1,
\[
\boxed{\mathcal T_q=
\frac{\rho^q\|H_MF\Theta^{q-2}\|_2\,\|D_*\|_2}{1-\rho t}}
\]
majore la queue uniformément. Toutes ses matrices sont préparées une fois.
Il ne faut pas utiliser ||F Theta^{q−2}D_*|| à la place du produit ci-dessus :
l'ordre des facteurs et la résolvante restante doivent être respectés.

Une première borne est la somme triangulaire
\[
\delta_q^{\rm triangle}
=\sum_{j=0}^{q-1}\rho^j\|H_MA_j\|_2+\mathcal T_q.
\]
Une borne pouvant conserver davantage de compensations transforme le
polynôme en Bernstein sur [0,rho]. Avec d=q−1 et 0≤k≤d, définir
\[
B_k=\sum_{j=0}^{k}
\frac{\binom{k}{j}}{\binom{d}{j}}\rho^j H_MA_j.
\]
Pour u=mu/rho∈[0,1],
\[
H_MP_d(\rho u)=
\sum_{k=0}^{d}\binom{d}{k}u^k(1-u)^{d-k}B_k.
\]
Les poids sont non négatifs et de somme un. La convexité de la norme donne
\[
\boxed{\delta_q^{\rm Bernstein}
=\max_{0\le k\le q-1}\|B_k\|_2+\mathcal T_q,\qquad
\sup_{\mu\in[0,\rho]}\|R(\mu)\|_{S_M,2}\le\delta_q.}
\]
Prendre le minimum de majorants établis, sur q ou entre méthodes, est
permis. Former d'abord chaque B_k, puis sa norme, est essentiel aux
compensations. Par exemple P(mu)=(mu−1/2)² sur [0,1] donne la borne
triangulaire 9/4 et les contrôles Bernstein (1/4,−1/4,1/4), donc 1/4.

La subdivision par de Casteljau peut resserrer la partie polynomiale
sans refaire de factorisation. La queue reste une borne uniforme séparée ;
ni quelques points fréquentiels ni une grille dense ne prouvent la borne.

Un autre certificat algébrique est possible dans de très petites
dimensions : multiplier par un dénominateur positif connu puis vérifier
des coefficients de Bernstein PSD. Le test rationnel fourni réalise
cette vérification sur toute la bande lorsque Theta est scalaire.
L'adjugée d'un grand petit système et ses polynômes de haut degré ne sont
pas recommandées ici : leur coût et leur stabilité restent à examiner.

## 5. Base imparfaitement contrainte et réparation

Une petite valeur B^T V calculée n'est pas une contrainte exacte.
Choisir un right-inverse théorique C tel que B^T C=I et poser
\[
\Delta=B^T V,\quad H=C\Delta,\quad \bar V=V-H.
\]
Si l'on dispose seulement de C0 flottant, on peut définir mathématiquement
C=C0(B^T C0)^{-1}, sous rang établi. Une résolution flottante de ce petit
système n'est pas une représentation exactement admissible par elle-même.

Garder les coefficients Y existants et définir
\[
X=T-EVY,\qquad \bar X=T-E\bar VY=X+EHY.
\]
Le champ réparé possède les mêmes ports et covecteurs retenus que T.
Il est inutile de prétendre que les coefficients étaient Galerkin sur
bar V : le résidu de ce champ suffit. Pour l'enveloppe, remplacer simplement
P,Q par K bar V et lambda_c M bar V dans les définitions de la section 2.
Les mêmes Theta,a0,a1,Y restent utilisables sous rho||Theta||<1.

On obtient ainsi un majorant delta du résidu de bar X. En notant
\[
Y_{\max}=\frac{\|a_0\|_2+\rho\|a_1\|_2}{1-\rho t},
\qquad \alpha_*=\lambda_c-\Omega^2,
\]
les constantes d'erreur du champ réellement renvoyé X peuvent être prises
comme
\[
\boxed{c_M=\delta/\alpha_*+\|H\|_{M,2}Y_{\max},\quad
c_D=\sqrt{\lambda_c}\,\delta/\alpha_*+\|D_iH\|_2Y_{\max}.}
\]
Pour une réponse précise, les normes du correctif H Y q peuvent remplacer
ces majorants d'opérateur plus grossiers. Elles doivent porter sur les
champs calculés, conformément au contrat de petites combinaisons de charges.

La composante résiduelle réparée donne l'erreur de Schur énergétique
\[
0\preceq\bar X^T A\bar X-S_{\rm exact}
\preceq\delta^2 I/\alpha_*.
\]
Cela ne vaut pas directement pour X^T A X. Si le Schur renvoyé conserve
cette dernière énergie, ajouter l'écart de fonctionnel. Avec les matrices
complètes, une borne uniforme explicite est
\[
\epsilon_{\rm contrainte}\le
2(\|T^T K_{\rm full}E H\|+
\Omega^2\|T^T M_{\rm full}E H\|)Y_{\max}
\]
\[
+\big(\|H^T K H-V^T K H-H^T K V\|
+\Omega^2\|H^T M H-V^T M H-H^T M V\|\big)Y_{\max}^2.
\]
Les produits se calculent avec DT, D_iV, D_iH et la masse, sans former
K_full. Ils suivent de l'expansion de bar X=X+EHY. Une borne plus simple
par produits de normes est permise mais peut être beaucoup plus large.
Si le programme calcule directement l'énergie de bar X, ce terme n'est
pas nécessaire ; ses erreurs d'évaluation, elles, le restent.

## 6. Différence K_Q/K_D et erreurs du petit solve

Le facteur statique Q sert à construire les directions, pas à redéfinir
la physique. Pour K_Q=E_s^{-1}R_Q^T R_Q E_s^{-1}, l'application contrainte
exacte peut s'écrire
\[
S_Q=E_s R_Q^{-1}(I-Q_KQ_K^T)R_Q^{-T}E_s,\quad
Q_K=orth(R_Q^{-T}E_s B).
\]
Elle permet les semences S_Q G0, S_Q(lambda_c G1), puis les applications
S_Q(lambda_c M V), sans nouvelle factorisation. Le remplissage du facteur,
l'orthogonalisation et la largeur de V restent à mesurer.

Même si K_D≽(1−eta)K_Q fournit le certificat de lambda_c, les matrices
G0,P et les énergies doivent viser D d'entrée. Si l'on emploie un autre
schéma de petits coefficients, l'identité résiduelle générale reste vraie
avec le modèle D dans P,Q,E0,E1. Un résidu nul dans Q peut être non nul
dans D ; le test diag(K_Q)=(1,2,3), diag(K_D)=(1,3,3) l'éprouve exactement.

L'enveloppe décrite vise les coefficients idéaux de la petite résolvante.
Si Yhat est une résolution approchée et
tau=(a0+mu a1)−(I−mu Theta)Yhat, alors
\[
Yhat-Y=-(I-\mu\Theta)^{-1}\tau,
\]
\[
Rhat-R=(P-\mu Q)(I-\mu\Theta)^{-1}\tau.
\]
D'où le terme additionnel
\[
\frac{\|H_M P\|+\rho\|H_M Q\|}{1-\rho t}\|\tau\|.
\]
Il faut aussi tenir compte des erreurs de formation des résidus, produits
et normalisations pour annoncer une certification machine. Les observer
en doubles et rendre certification_machine=False décrit une évaluation
conditionnelle, pas une preuve complète des arrondis.

## 7. Schur rapide, assemblage et marge

Les projections constantes de T,EV ou E bar V permettent de calculer
l'énergie candidate avec de petites matrices :
\[
\widehat S=A_{00}-A_{0V}Y-Y^TA_{0V}^T+Y^TA_{VV}Y.
\]
La simplification A00−G_V^T Y est réservée au Galerkin exact associé.
Les défauts du fonctionnel réellement renvoyé s'ajoutent à
b=delta²/alpha_* et, si nécessaire, epsilon_contrainte.

Après assemblage, q regroupe tous les ports et toutes les coordonnées
retenues. Si ||Shat−S||≤b et sigma_min(Shat)≥sigma>b,
\[
\|q-qhat\|\le
(\|f-Shat\,qhat\|+b\|qhat\|)/(\sigma-b).
\]
Les extensions de champ utilisent ce q entier et les normes physiques
complètes. Une marge des seuls ports peut rester positive alors que le
bloc retenu possède un mode singulier ; la contre-épreuve diagonale le
montre. Le certificat du complément ne supprime aucune résonance globale.

Les douze tests Fraction contrôlent les signes, la queue, les identités
Bernstein, un certificat polynomial sur toute la bande, la cancellation
de Grams, la réparation et son énergie, le défaut D/Q et la marge du bloc
entier, ainsi que l'action anisotrope et ses contre-exemples. Ils complètent
la preuve sans établir un solveur machine certifié.


Le [contrôleur expérimental](../ci/controle_complement.py) met en œuvre
ces expressions, avec la [variante anisotrope](KRYLOV_CONTRAINT_ANISOTROPIE.md).
Les [contre-épreuves Fraction](../ci/test_krylov_contraint_preuves.py)
et les [tests physiques](../ci/test_krylov_contraint.py) sont indépendants
des petites matrices projetées utilisées par le candidat.
