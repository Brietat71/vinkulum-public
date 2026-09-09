# Retenir les résonances et condenser seulement leur complément

Note de composition, 8 septembre 2026. Dérivation indépendante à partir des
contrats de Vinkulum 0.11.0 (9062db8) et du carnet
[sur la trace complémentaire](TRACE_COMPLEMENT_SPECTRAL_PREUVES.md). Les identités ci-dessous sont en
arithmétique exacte. Un certificat dirigé du minorant spectral ne transforme
pas les résolutions flottantes et les bornes de champ en calculs machine
certifiés. Aucune nouveauté théorique ni domination d'un solveur multicorps
complet n'est revendiquée.

## 1. Le contrat qui change

Le champ physique se partitionne en intérieur i de taille n et interface s
de taille p. Poser
\[
\mathcal K=D^T D,\quad \mathcal A(z)=\mathcal K-z\mathcal M,\quad
K=\mathcal K_{ii},\quad M=\mathcal M_{ii},\quad z=\omega^2.
\]
On suppose K,M SPD, la masse complète PSD pour les normes physiques, et
B de taille n×r, de rang r, avec 0<r<n. Le sous-espace éliminé est
\[
N=\ker B^T,\qquad v^T K v\ge\lambda_c v^T Mv\quad(v\in N),
\qquad \Omega^2<\lambda_c.
\]
On ne suppose plus K−zM SPD sur tout l'intérieur. Les r coordonnées
complémentaires restent dans le système assemblé. Il est faux de remplacer
seulement lambda_min dans ControleChamp : à seuls ports fixés, l'erreur
peut suivre les directions retenues et n'appartient alors pas à N.

Choisir Phi telle que J=B^T Phi soit inversible. Avec un relèvement
arbitraire Psi, une normalisation de port inversible W et E l'injection
intérieure, définir
\[
L=\begin{bmatrix}\Psi\\W\end{bmatrix},\quad T=[L,E\Phi],\quad
q=\begin{bmatrix}y\\a\end{bmatrix},\quad x=Tq+Ew,\quad w\in N.
\]
C'est une décomposition bijective :
\[
y=W^{-1}x_s,\quad a=J^{-1}B^T(x_i-\Psi y),\quad
w=x_i-\Psi y-\Phi a.
\]
Ni Phi ni Psi ne sont supposés exacts au sens modal ou statique. Les
échelles des nouvelles coordonnées doivent être déclarées : la norme
euclidienne du bloc (y,a) dépend de ces unités.

Avec U=K^{-1}B et H=B^T U, le relèvement canonique C=UH^{-1} vérifie
B^T C=I. Choisir Phi=C permet a=B^T(x_i−Psi y), sans hypothèse B=M Phi.
Ce choix conserve les covecteurs sélectionnés par B, avec des représentants
qui ne sont pas nécessairement modaux. Il vérifie v^T KC=0 pour v∈N.
Si l'on remplace K_D par K_Q, cette orthogonalité ne permet pas de supprimer
le couplage dans le modèle D d'entrée. Un autre relèvement admissible est
C=Phi(B^T Phi)^{-1}; les deux diffèrent par une application dans N.

## 2. Inverse contrainte sans base dense Z

Définir S=K^{-1}−UH^{-1}U^T. Pour la preuve uniquement, une base Z de N donne
\[
S=Z(Z^T KZ)^{-1}Z^T,\qquad
R_z=Z[Z^T(K-zM)Z]^{-1}Z^T.
\]
L'inverse contrainte R_z existe sur la bande, même lorsque K−zM est
singulière sur l'intérieur complet. Sans construire Z,
\[
\boxed{R_z=(I-zSM)^{-1}S=S(I-zMS)^{-1}.}
\]
Preuve : les équations contraintes donnent w−zSMw=Sf. Réciproquement
cette équation impose w∈N et la même équation variationnelle. La
coercivité restreinte assure l'unicité. Sur N, SM est auto-adjoint positif
dans la métrique K, de norme au plus 1/lambda_c. Ainsi
\[
R_z=\sum_{j\ge0}z^j(SM)^jS\quad(z<\lambda_c).
\]
Cela motive un Krylov contraint; cela ne recommande pas de former la
matrice dense ambiante I−zSM. Un témoin indépendant sans Z est le KKT
\[
\begin{bmatrix}K-zM&B\\B^T&0\end{bmatrix}
\begin{bmatrix}w\\\mu\end{bmatrix}
=\begin{bmatrix}f\\0\end{bmatrix},\qquad w=R_z f.
\]
Le KKT est inversible sous les hypothèses ci-dessus. Sa factorisation
fréquentielle reste une résolution complète, pas un gain de réduction.

## 3. Le petit système conserve tous les couplages

Poser
\[
G(z)=E^T\mathcal A(z)T=G_0-zG_1,\quad
G_0=D_i^T(DT),\quad G_1=(\mathcal MT)_i,\quad
A_0=T^T\mathcal A T.
\]
Les colonnes comprennent les ports ET Phi. Le bloc condensé p+r est
\[
\boxed{S_c=A_0-G^T R_zG,\qquad X=T-E R_zG,\qquad x=Xq.}
\]
On a x_s=Wy, B^T(x_i−Psi y)=Ja et S_c=X^T A X.
Pour V⊂N de rang plein, poser
\[
H_V=V^T(K-zM)V,\quad Y=H_V^{-1}V^T G,\quad
\widehat X=T-EVY,\quad
\widehat S_c=A_0-G^T VH_V^{-1}V^T G.
\]
H_V est SPD et Shat_c=Xhat^T A Xhat. Garder Phi^T K V,
Phi^T M V, les couplages de masse du port et tous les termes de Psi.
B=M Phi exact annule Phi^T M V, mais généralement pas Phi^T K V.
Les produits D Phi et D V permettent les projections physiques sans
former D^T D. Si Y n'est pas le Galerkin exact, garder l'énergie candidate
complète; tout autre Schur renvoyé ajoute son écart de fonctionnel à la borne.

Pour une force complète f, avec f_i=E^T f,
\[
x=Xq+E R_z f_i,\qquad
S_cq=T^T f-G^T R_z f_i.
\]
Pour les seuls ports, le second membre est (W^T f_s,0). Des charges
intérieures nécessitent cette solution particulière et un audit de son
résidu; l'API actuelle ne les prend pas en charge.

## 4. Bornes à coordonnées conservées imposées

Soit Xhat un champ de mêmes ports ET mêmes covecteurs retenus que X.
F=Xhat−X=E F_i avec B^T F_i=0. Son résidu R=E^T A Xhat peut contenir
des réactions dans range(B), même si Xhat=X. L'identité correcte est
\[
\boxed{F_i=R_z R,}
\]
et non (K−zM)^{-1}R. La norme duale de N pour la masse s'écrit sans Z :
\[
S_M=M^{-1}-M^{-1}B(B^T M^{-1}B)^{-1}B^T M^{-1},\qquad
G_R=R^T S_M R.
\]
En effet S_M=Z(Z^T M Z)^{-1}Z^T, donc
r^T S_M r=sup_{v∈N,v≠0}(v^T r)^2/(v^T M v).
Elle annihile les réactions dans range(B). Comme 0≼S_M≼M^{-1},
le Gram complet R^T M^{-1}R reste une majoration valide mais pessimiste.
Il n'autorise pas un défaut des contraintes.

Avec alpha=lambda_c−z>0, les preuves spectrales sur la paire restreinte donnent
\[
\boxed{F^T\mathcal MF\preceq G_R/\alpha^2,\quad
F_i^T(K-zM)F_i\preceq G_R/\alpha,\quad
(DF)^TDF\preceq\lambda_c G_R/\alpha^2.}
\]
La dernière borne suit de K=(K−zM)+zM. Un minorant lambda_c suffit :
lambda/(lambda−z)^2 décroît pour lambda>z≥0. La stationnarité suivant N donne
\[
\boxed{\widehat X^T\mathcal A\widehat X-S_c
=F_i^T(K-zM)F_i=R^T R_z R\succeq0.}
\]
Une enveloppe uniforme ||R||_{S_M,2}≤delta donne
c_M=delta/(lambda_c−Omega²),
c_D=sqrt(lambda_c) delta/(lambda_c−Omega²) et
b_S=delta²/(lambda_c−Omega²).
Les développements Krylov doivent porter sur l'inverse contrainte et les
p+r colonnes du nouveau relèvement. Un échantillonnage ne prouve pas une
enveloppe uniforme.

Un défaut flottant B^T V petit n'est pas une égalité exacte.
Avec B^T C=I, la réparation mathématique Vbar=V−C(B^T V) appartient à N.
Pour Delta=B^T(Xhat_i−T_i), Xbar=Xhat−E C Delta a les bonnes coordonnées.
Les bornes s'appliquent à Xbar; borner Xhat demande d'ajouter les normes de
E C Delta et le défaut de fonctionnel. Ces opérations nécessitent elles
aussi un encadrement pour prétendre à une certification machine.

## 5. Assemblage et singularités

Chaque sous-structure fournit q_j=A_j q_global. q_j comprend ports et
coordonnées retenues locales, généralement privées : elles ajoutent des
colonnes à l'assemblage et ne sont pas collées comme des ports entre pièces.
Conserver la conformité des interfaces et la comptabilité des masses.

Les congruences et majorants de Loewner s'assemblent comme avant. La marge
globale porte toutefois sur TOUT le bloc conservé. Si
||Shat−S||₂≤b, sigma_min(Shat)≥sigma>b et r=f−Shat qhat, alors
\[
\|q-qhat\|_2\le(\|r\|_2+b\|qhat\|_2)/(\sigma-b).
\]
Les extensions du champ portent également sur tout q. Une seconde
élimination des coordonnées retenues nécessite une nouvelle justification;
elle peut échouer précisément aux anciennes résonances intérieures.

La congruence par [T,E Z] est inversible. L'inertie du système complet est
celle de S_c plus n−r valeurs positives. Le système complet est inversible
ssi S_c l'est. Une résonance globale n'est donc pas guérie. À marge nulle
ou non établie, refuser la borne de réponse. Un amortissement arbitraire ou
une pseudo-inverse changerait le problème.

## 6. Contre-exemples

1. Couplage de raideur omis : K=[[2,1],[1,3]], M=I, Phi=B=e1, z=1.
   lambda_c=3; Schur de la direction retenue=1−1/2=1/2.
   Supprimer le couplage donne 1 et divise sa réponse par deux.
2. Pôle intérieur traversable : dans l'ordre (i1,i2,port),
   Kfull=[[1,0,1],[0,4,0],[1,0,2]], Mfull=I, B=Phi=e1, z=1.
   L'intérieur diag(0,3) est singulier; le complément vaut 3.
   Le bloc conservé (port,a) est [[1,1],[1,0]], déterminant −1.
   Sous force de port unitaire, le champ unique est (1,0,0).
3. Pôle global conservé : Kfull=diag(1,4,2), Mfull=I, B=Phi=e1, z=1.
   Même complément coercif; le bloc conservé diag(1,0) est singulier.
4. Mauvaise paire : B=e1, Phi=e2 en dimension deux ont chacun rang un,
   mais B^T Phi=0. La décomposition est invalide.
5. Mauvais espace de borne : K=diag(epsilon,1), M=I, B=e1.
   lambda_c=1 mais e=e1∉N. À z=0, son résidu vaut epsilon e1.
   La borne mal appliquée annoncerait ||e||≤epsilon alors que ||e||=1.
6. Réaction : r=B mu peut être non nul pour un champ exact contraint,
   mais R_z r=S_M r=0. La norme massique duale complète surestime ce
   défaut; supprimer un résidu quelconque sous prétexte de réaction est faux.

## 7. Binary64 et chemin algorithmique

Contrat recommandé : les valeurs stockées de B définissent EXACTEMENT
les covecteurs cibles. Si B=fl(M Phi), ne pas invoquer Phi^T M Z=0.
Garder tous les couplages suffit aux identités pour ce B, sous J inversible.
Viser le produit mathématique M Phi exige son encadrement et la même
convention au certificat, aux projections, aux relèvements et aux résidus.

Si K_D≽(1−eta)K_Q et la trace tau_c,Q est encadrée sur le MÊME B,
lambda_c=(1−eta)/tau_c,Q minore la paire restreinte du modèle D.
Cela ne donne pas R_z,D=(I−z S_Q M)^{-1}S_Q, qui viserait le modèle Q.
Le défaut K_D−K_Q doit rester dans le résidu et les projections physiques.

Former S f=K^{-1}f−UH^{-1}U^T f peut annuler deux très grandes réponses.
Pour K=T_K^T T_K, une forme équivalente est
\[
S=T_K^{-1}(I-QQ^T)T_K^{-T},\quad Q=orth(T_K^{-T}B).
\]
Un QR de Householder de T_K^{-T} B projette en coordonnées d'énergie,
avant la seconde résolution, sans base dense du complément. Cela limite
certaines cancellations mais ne prouve aucune erreur relative sur une
composante presque nulle. Contrôler résidus et contraintes reste nécessaire.
Le Gram B^T K^{-1}B peut mettre au carré le conditionnement; la preuve de
rang et la résolution sont distinctes de la sélection modale initiale.

Stockage O(nr) pour K^{-1}B ou les réflecteurs, plus les facteurs creux.
Le remplissage dépend du graphe; aucune complexité linéaire universelle.

Ordre de validation : certificat dirigé et refus de rang; petites identités
rationnelles avec masse couplée et Phi non modal; témoin KKT sans Z aux
pôles intérieurs traversables et aux vrais pôles globaux; Krylov contraint
ensuite; confrontation 0–40 Hz aux mêmes oracles et charges, préparation
comprise. La seule extension de bande ne prouve aucun gain de temps.


Le [certificat dirigé et la qualification 0–40 Hz](COMPLEMENT_SPECTRAL_DIRIGE.md)
mettent maintenant en œuvre les deux premières étapes. La [preuve par
inertie](INERTIE_COMPLEMENT_PREUVES.md) ouvre une autre voie pour le seuil
spectral. Les [dix contre-épreuves de composition](../ci/test_retention_complement.py)
sont exécutées dans la CI ; le témoin reste hors API publique.
