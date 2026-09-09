# Inertie du complément contraint sans inverse du complément

Note du 8 septembre 2026. Résultat algébrique classique redérivé pour le
contrat de Vinkulum, éprouvé par dix tests Fraction dans
[ci/test_inertie_complement.py](../ci/test_inertie_complement.py). Cette note
conserve la dérivation initiale après le certificat de trace. Son
[transfert par arrondis dirigés](INERTIE_DIRIGEE_PREUVES.md) et la
[campagne de 84 essais](INERTIE_CONTRAINTE_PROTOTYPE.md) sont désormais
disponibles séparément ; les dix tests ci-dessous restent algébriques.

Soient A une matrice réelle symétrique n×n, B de taille n×s de rang colonne
plein, et Z une base de N=ker(B^T). A n'a pas à être définie positive et
Z^T A Z peut être singulière. L'inertie désigne, dans cet ordre, les nombres
de valeurs propres positives, négatives et nulles.

Choisir Y tel que B^T Y=I_s, par exemple Y=B(B^T B)^{-1} pour la preuve.
Le changement x=Zz+Yy est bijectif. Le formulaire quadratique de
\[
C=\begin{bmatrix}A&B\\B^T&0\end{bmatrix}
\]
devient
\[
z^T Z^T A Zz+2z^T Z^T A Yy+y^T Y^T A Yy+2y^T\mu.
\]
Avec le changement triangulaire inversible
\[
\nu=\mu+Y^T A Zz+\tfrac12Y^T A Yy,
\]
il vaut exactement
\[
z^T Z^T A Zz+2y^T\nu.
\]
Aucun inverse de Z^T A Z n'est employé. Le dernier terme porte la matrice
hyperbolique [[0,I],[I,0]], d'inertie (s,s,0). La loi de Sylvester donne donc
\[
\boxed{\operatorname{Inertia}(C)
=\operatorname{Inertia}(Z^T A Z)+(s,s,0).}
\]
Le test de congruence explicite vérifie cette égalité matricielle avec
contraintes obliques, masse couplée, matrices indéfinies et complément nul.

Pour A_gamma=K−gamma M et une masse définie positive sur N, on obtient
\[
\boxed{(K-\gamma M)|_N\succ0
\iff\operatorname{Inertia}(C_\gamma)=(n,s,0)
\iff\gamma<\lambda_{\min}(K|_N,M|_N).}
\]
Un test d'inertie positif fournirait ainsi un minorant spectral strict, sans
soustraire deux traces. À la frontière singulière, les zéros restent visibles
dans le KKT ; la formule ne perd pas son sens.

Pour K=diag(10^{-30},1,2,3), M=I et B=e1 :

| gamma | Inertie du complément | Inertie du KKT | Coercivité stricte |
|---|---|---|---|
| 1/2 | (3,0,0) | (4,1,0) | oui |
| 3/2 | (2,1,0) | (3,2,0) | non |
| 1 | (2,0,1) | (3,1,1) | non, frontière |

L'intérieur complet est déjà indéfini au premier point ; cela n'empêche
pas la coercivité du complément. La trace totale 10^30+11/6 et le terme
retiré 10^30 s'annulent en binary64, tandis que le test exact d'inertie
reste une décision indépendante de cette soustraction.
Une autre contre-épreuve accepte gamma=3/4, alors que le minorant par
trace vaut seulement 6/11. L'inertie peut donc confirmer une bande que ce
minorant ne suffit pas à accepter, même avec des traces évaluées exactement.

La condition de rang n'est pas facultative. Pour B=[e1,2e1], le KKT possède
un multiplicateur nul. Au premier point sa signature est (4,1,1).
Ajouter (s,s,0) avec s=2 à celle de ker(B^T) donnerait même une dimension
incorrecte. Plus généralement, si le rang réel est r<s, la correction est
(r,r,s−r), après séparation des multiplicateurs redondants. La formule
simple proposée doit refuser un B redondant plutôt que masquer ces zéros.

Le calcul indépendant utilisé par les tests procède par congruence exacte :
un pivot diagonal non nul apporte son signe puis son complément de Schur ;
si toute la diagonale restante est nulle, un coefficient hors diagonale
b non nul fournit le pivot [[0,b],[b,0]], de signature (1,1,0), dont
l'inverse exacte est [[0,1/b],[1/b,0]]. Une matrice restante nulle donne sa
nullité. Ce procédé termine sans racines, valeurs propres flottantes,
tolérance numérique ou hypothèse d'inversibilité globale. Les permutations
physiques, celles mêlant multiplicateurs et physique, et les changements
inversibles de base des contraintes préservent l'inertie vérifiée.

Cette preuve suggérait une alternative à tester lorsque la trace dirigée
perd trop de chiffres. Elle ne prouve pas qu'une LDL flottante ordinaire
certifie l'inertie. Les obligations de transfert sont :

- garder le même B exact stocké, le même M et le même modèle K dans le
  certificat et la réduction ;
- encadrer les pivots et les erreurs de factorisation congruente, y compris
  les blocs 2×2, et refuser toute signature non établie ;
- accroître la précision ou modifier le point gamma lorsque les intervalles
  ne permettent pas de séparer les signes, sans transformer un intervalle
  contenant zéro en pivot régulier ;
- mesurer remplissage, coût des décisions répétées et mémoire ; le calcul
  Fraction de petite dimension ne fournit aucune complexité creuse.

Si le certificat vise K_Q et qu'un encadrement séparé établit
K_D≽(1−eta)K_Q, eta<1, une coercivité certifiée de K_Q−gamma M sur le même
N implique celle de K_D−(1−eta)gamma M. Pour couvrir une bande Omega,
on peut choisir gamma≥Omega²/(1−eta), avec quotient encadré vers le haut,
puis établir la coercivité stricte à ce point ; l'égalité du point choisi
avec ce quotient est permise puisque le test d'inertie est strict. Un
point supérieur laisse une réserve explicite. Cette composition exige l'audit
énergétique D/Q ; former fl(D^T D) et certifier son inertie n'établirait pas
automatiquement celle du produit mathématique D^T D.

Aucune nouvelle API ou extension publique de bande n'est livrée par ces
tests. Le certificat machine et ses mesures appartiennent au
[prototype ultérieur](INERTIE_CONTRAINTE_PROTOTYPE.md).

Voir l’[audit des manuscrits de 2026](VERIFICATION_CREUSE_2026.md) pour les
méthodes de vérification creuse envisagées et leurs conditions de transfert.
