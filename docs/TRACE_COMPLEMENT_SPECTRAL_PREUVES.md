# Élargir la bande en conservant les directions résonantes

Ce carnet dérive une voie pour dépasser le refus de bande de la réduction
matérielle 0.10.0, encore présent dans la 0.11.0. Ses identités sont désormais
éprouvées par un [certificat dirigé expérimental](COMPLEMENT_SPECTRAL_DIRIGE.md)
et un témoin mécanique sur 0–40 Hz. Les [preuves de composition](RETENTION_INTERIEURE_PREUVES.md)
traitent les coordonnées retenues, les couplages et l'erreur du champ.
Cette extension reste **hors API publique** ; les confrontations des roues
0.10.0 et 0.11.0 conservent leurs résultats et leurs refus.

## 1. Le problème à résoudre

Le certificat actuel exige que l'opérateur intérieur éliminé soit coercif
sur toute la bande demandée. Le minorant par trace porte sur son premier
mode. Lorsque la bande dépasse une fréquence propre intérieure, améliorer
uniquement la précision de ce minorant ne peut pas rendre l'opérateur
coercif.

Une autre décomposition est possible : conserver quelques coordonnées
intérieures avec les ports, puis n'éliminer que leur complément. La
positivité à démontrer porte alors sur un sous-espace différent. Cela ne
supprime aucune résonance physique du système complet.

## 2. Inverse contrainte sans construire une base globale du complément

Soient \(K,M\) symétriques définies positives, \(B\in\mathbb R^{n\times s}\)
de rang colonne plein, \(0<s<n\), et \(Z\) une base quelconque de
\(\ker B^\mathsf T\). Posons

\[
H=B^\mathsf TK^{-1}B,\qquad
S=K^{-1}-K^{-1}BH^{-1}B^\mathsf TK^{-1}.
\]

Alors

\[
\boxed{S=Z(Z^\mathsf TKZ)^{-1}Z^\mathsf T.}
\]

**Preuve.** Pour toute charge \(f\), les équations stationnaires de
\(\frac12x^\mathsf TKx-f^\mathsf Tx\), sous \(B^\mathsf Tx=0\), donnent
\(Kx+B\mu=f\), puis \(H\mu=B^\mathsf TK^{-1}f\) et \(x=Sf\).
En écrivant \(x=Zy\), elles donnent aussi
\(y=(Z^\mathsf TKZ)^{-1}Z^\mathsf Tf\). La stricte convexité assure
l'unicité, donc les deux opérateurs sont égaux. En particulier \(S\)
est positive semi-définie de rang \(n-s\), bien qu'elle soit écrite
comme une différence.

La trace massique se calcule ainsi :

\[
\boxed{\tau_c=\operatorname{tr}(MS)
=\operatorname{tr}(MK^{-1})
-\operatorname{tr}\!\left[H^{-1}
(K^{-1}B)^\mathsf TM(K^{-1}B)\right].}
\]

Les inverses dans ces expressions désignent des résolutions. La première
trace est celle déjà étudiée par inverse sélectionnée. Le second terme
demande \(s\) résolutions, une projection massique et un système \(s\times s\).
Aucune inverse globale complète ni base dense \(Z\) n'est nécessaire à
cette identité de trace. Le remplissage des facteurs et le coût des
résolutions restent dépendants du graphe ; ce n'est pas une preuve de
complexité linéaire pour un mécanisme quelconque.

## 3. Borne sur le complément effectivement choisi

Les valeurs propres généralisées \(\mu_j>0\) de
\((Z^\mathsf TKZ,Z^\mathsf TMZ)\) vérifient

\[
\tau_c=\sum_{j=1}^{n-s}\frac1{\mu_j},
\qquad
\boxed{\mu_{\min}\ge\frac1{\tau_c}.}
\]

En effet, la trace ci-dessus est celle de
\((Z^\mathsf TKZ)^{-1}(Z^\mathsf TMZ)\), dont les valeurs propres sont
\(1/\mu_j\). Chaque terme positif est inférieur ou égal à leur somme.
Pour \(s=0\), on retrouve la borne actuelle. Pour un complément de
dimension un, la borne est exacte.

**Il n'est pas nécessaire que les colonnes choisies soient des modes propres
exacts.** Le résultat porte sur le noyau du \(B\) réellement utilisé.
Des directions proches des modes bas peuvent améliorer la borne ; elles
ne dispensent jamais de calculer et vérifier cette borne. Une valeur propre
numérique omise ne constitue pas, à elle seule, un minorant certifié.

Si \(K_D=D^\mathsf TD\) et une factorisation approchée \(K_Q\) satisfont
\(K_D\succeq(1-\eta)K_Q\), \(0\le\eta<1\), cette inégalité se restreint
au même complément. Par conséquent,

\[
\lambda_{\min}(Z^\mathsf TK_DZ,Z^\mathsf TMZ)
\ge \frac{1-\eta}{\tau_{c,Q}},
\]

où \(\tau_{c,Q}\) est calculée avec \(K_Q\), \(M\) et **le même** \(B\).
Cela raccorde la formule à la perturbation énergétique du
[certificat actuel](INVERSE_SELECTIONNEE_PREUVES.md).

## 4. Ce qui doit rester dans le système assemblé

Avec des directions retenues \(\Phi\) et un complément \(Z\), la congruence
par \([\Phi,Z]\) doit conserver **tous** les blocs de \(K-\omega^2M\).
Si l'on choisit exactement \(B=M\Phi\), le couplage massique
\(\Phi^\mathsf TMZ\) est nul. Le couplage de raideur ne l'est généralement
pas pour des modes approchés : il doit rester dans la condensation.

La condition \(\omega^2<(1-\eta)/\tau_{c,Q}\) justifie l'élimination du
complément. Elle ne garantit ni l'inversibilité du Schur restant, ni la
précision de ses réponses près des résonances globales. Les directions
retenues ajoutent \(s\) inconnues au problème couplé ; les retirer de la
physique pour obtenir une matrice positive serait un autre modèle.

Si \(B\) est stocké en binary64, il faut annoncer si la contrainte cible
est celle de ses valeurs exactes ou celle du produit mathématique \(M\Phi\).
Arrondir ce produit puis employer les deux conventions indistinctement
invaliderait une preuve de composition.

## 5. Obstacle numérique : soustraire deux grandes traces

Prenons \(K=\operatorname{diag}(10^{-30},1,2,3)\), \(M=I\), \(B=e_1\).
La trace totale vaut \(10^{30}+11/6\) et le terme retiré vaut \(10^{30}\).
La trace complémentaire exacte est \(11/6\), donc le minorant vaut \(6/11\).
En binary64, soustraire les deux grandes traces arrondies donne zéro.
Cela pourrait produire une fausse borne infinie.

La mise en production demande donc des **encadrements dirigés** : si
\(\tau\in[\tau_-,\tau_+]\) et le terme retiré
\(\chi\in[\chi_-,\chi_+]\), alors
\(\tau_c\in[\tau_- -\chi_+,\tau_+ -\chi_-]\). Il faut vérifier le rang,
la positivité des petits systèmes et une borne supérieure strictement
positive de \(\tau_c\), puis arrondir le quotient final vers le bas.
Une cancellation qui rend l'encadrement trop large impose davantage de
précision ou un refus ; elle n'autorise ni troncature à zéro ni réparation
arbitraire du signe.

Les [contre-épreuves en fractions exactes](../ci/test_trace_complement.py)
vérifient l'identité pour des contraintes obliques, une masse couplée,
la borne de coercivité, sa restriction après perturbation énergétique,
la perte par soustraction et le refus d'une contrainte redondante.
Elles n'établissent pas encore un algorithme machine certifié à grande
échelle.

La prochaine expérience doit conserver les facteurs matériels, mesurer
le coût de ces résolutions supplémentaires, puis confronter les champs
sur la bande 0–40 Hz avec le même oracle et les mêmes six charges.
Ce programme emploie des identités classiques de contraintes et de Schur ;
aucune nouveauté théorique mondiale ni avance sur un solveur complet
n'est revendiquée par ce carnet.
