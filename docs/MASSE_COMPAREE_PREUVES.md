# Masses couplées : transporter les normes sans modifier la physique

8 septembre 2026. Prototype de recherche, hors API publique 0.11.0.

Le contrôleur admet maintenant une masse intérieure dont les composantes
connexes dépassent six coordonnées. Il garde la matrice complète dans les
équations et prouve une comparaison avec une métrique auxiliaire diagonale.
Sur trois consoles à masse consistante, les 24 processus passent le seuil
physique de 10⁻⁶ ; les 12 processus avec contrôle passent leurs majorants.
C'est un progrès de couverture du prototype et un avantage de temps sur
le grand cas. Aucune supériorité généraliste sur Exudyn, MBDyn ou Simpack
n'est établie par cette expérience.

## 1. L'obstacle et l'invariant physique

`ControleFacteurs` utilisait `RacineMasseBlocs`, limité aux composantes
connexes de taille au plus six. Cela convient à des inerties de corps
individuels, mais refuse les couplages cinétiques entre sections d'une
structure flexible. La réduction de Krylov conservait déjà ces couplages ;
le calcul des normes auxiliaires bloquait leur admission.

Remplacer M par diag(M) dans K−ω²M changerait les fréquences et les champs.
Par exemple, avec K=diag(4,9), M=[[2,1],[1,2]], ω=1 et f=(1,0),
la réponse exacte est (7/13,1/13). La masse diagonalisée donne (1/2,0).
Le test rationnel conserve ce contre-exemple.

Le nouveau contrôleur garde M dans la sélection, la base, la condensation,
les produits dynamiques, les résidus, les réparations bilinéaires, la
reconstruction et les normes physiques retournées. Seuls des majorants
auxiliaires changent.

## 2. Proposition de transport sur un quotient contraint

Soit M=Mᵀ>0, L=Lᵀ>0 et 0<α<β tels que

\[
\alpha L\preceq M\preceq\beta L.
\]

Soit B de rang plein et Z une base quelconque de ker(Bᵀ). Définissons

\[
S_M=Z(Z^TMZ)^{-1}Z^T,\qquad S_L=Z(Z^TLZ)^{-1}Z^T.
\]

La congruence par Z conserve l'ordre ; l'inversion des matrices positives
le renverse. Une dernière congruence donne donc

\[
\beta^{-1}S_L\preceq S_M\preceq\alpha^{-1}S_L.
\]

Ainsi, pour tout résidu f et tout déplacement u,

\[
\sqrt{f^TS_M f}\le\alpha^{-1/2}\sqrt{f^TS_L f},\qquad
\|u\|_M\le\sqrt\beta\,\|u\|_L.
\]

La première quantité est la norme duale de f restreint à ker(Bᵀ).
Elle est inchangée si l'on ajoute Bν à f : les réactions contraintes
s'annulent exactement. Aucun projecteur euclidien ne peut être substitué
sans justification au quotient dans sa métrique.

**Composition avec le contrôleur.** Les trois images duales du polynôme de
résidu sont multipliées par α⁻¹ᐟ² avant l'enveloppe de Bernstein. Les images
massiques de W et du correctif H sont multipliées par √β. Les identités
bilinéaires du défaut fonctionnel utilisent toujours M. Par conséquent,
les preuves du contrôleur précédent restent valables en arithmétique exacte,
avec des constantes éventuellement moins serrées. Un α trop faible peut
faire refuser la borne relative, même si le champ numérique est précis.

`ControleMasseComparee` réalise le quotient de L par Householder, sans
construire Z dense. L'inégalité d'ordre est classique ; la contribution
présente est sa composition vérifiée avec ce contrôleur et son coût mesuré.
Elle n'est pas présentée comme un nouveau théorème fondamental.

## 3. Ce qui est réellement prouvé par la machine

Le candidat L=diag(M) et les nombres binary64 α,β sont des propositions.
Le certificat vérifie séparément M−αL>0 et βL−M>0, par congruences creuses
à intervalles binary64 compilées. Les produits αL et βL sont encadrés,
et ne sont pas remplacés par leurs seules valeurs arrondies.

La cible est la matrice représentée exactement par les coefficients stockés.
Les contrôles de plateforme, les refus de pivot ambigu et les budgets du
[noyau d'intervalles](INERTIE_BINAIRE_SEPARATEURS.md) restent actifs.
L'égalité semi-définie à une frontière est refusée : l'implémentation exige
une comparaison stricte. Le coût des deux preuves est dans la préparation.
Une mauvaise comparaison n'est jamais corrigée silencieusement.

Cela ne certifie pas tous les arrondis des champs et majorants : QR,
Bernstein, SVD et résolutions du contrôleur restent numériques. La propriété
`certification_machine` de ses réponses demeure fausse. Une preuve sur M
n'est ni une preuve de précision du modèle continu ni une validation de ses
paramètres physiques.

## 4. Nouveau modèle cinétique et nouvelle référence

La rigidité matérielle provient des trois modèles natifs déjà identifiés
par SHA-256. La masse est une nouvelle discrétisation : on interpole
linéairement les vitesses de translation et de rotation des sections.
Pour un segment de masse mₑ et des inerties quadratiques de section J,

\[
M_e={m_e\over6}\begin{bmatrix}2&1\\1&2\end{bmatrix}
\otimes\operatorname{diag}(1,1,1,J_x,J_y,J_z).
\]

Les coordonnées encastrées sont supprimées après assemblage. Les termes
d'inertie longitudinale d'un parallélépipède rigide employés dans l'ancienne
masse concentrée n'appartiennent pas à cette intégrale des vitesses de section.
Cette expérience change donc bien le modèle cinétique du banc : les anciens
champs et ratios Exudyn ne sont pas réutilisés.

Par élément, les valeurs propres relatives à diag(Mₑ) sont 1/2 et 3/2.
La somme des formes et la restriction aux coordonnées libres préservent
ces inégalités. Sur les matrices stockées, le calcul vérifie explicitement
les propositions α=0,49 et β=1,51. L'identité d'intégration est testée sur
un champ de vitesses arbitraire et l'encastrement.

La sélection d'une direction utilise le problème généralisé symétrique
L⁻¹ᐟ² M L⁻¹ᐟ² v = μ L⁻¹ᐟ² K_Q L⁻¹ᐟ² v, puis Φ=L⁻¹ᐟ²v et B=MΦ.
La valeur approchée ne certifie aucune bande. Le certificat d'inertie cible
ensuite le B stocké exact et vérifie le complément au seuil de 80 Hz.

`OracleChampsCouples` résout directement le système par blocs en Decimal,
à 70 puis 90 chiffres. Il conserve les couplages de masse sur les blocs
hors diagonale. Sur les 257 fréquences et six charges de chacun des trois
cas, les champs convertis en doubles sont identiques entre ces précisions.
Ce test de convergence n'est pas un encadrement de l'erreur Decimal.
L'oracle est limité aux chaînes à couplages entre blocs voisins ; cette
restriction de validation ne vient pas du contrôleur.

### Audit indépendant du facteur du juge

L'ancien juge supposait lui aussi une masse diagonale. Le nouveau calcule
un facteur de Cholesky en bande R, puis **exactement en rationnels**
E=RᵀR−M. Une congruence diagonale dyadique T donne

\[
\delta=\max_i\sum_j|(TET)_{ij}|,\quad
\ell_0=\min_i (TLT)_{ii},\quad \rho=\delta/\ell_0.
\]

La borne de norme d'une matrice symétrique par sa somme absolue maximale
implique |uᵀEu|≤ρuᵀLu. Une nouvelle comparaison αⱼL<M, effectuée par le
juge sur la masse physique complète avec αⱼ=1/4, entraîne

\[
(1-\varepsilon)M\preceq R^TR\preceq(1+\varepsilon)M,
\qquad\varepsilon=\rho/\alpha_j<1.
\]

Le juge multiplie donc l'erreur relative mesurée dans R par un arrondi
supérieur de √((αⱼ+ρ)/(αⱼ−ρ)). La correction vaut au moins 1. Le cas E=0
est valable et testé séparément du cas non exact. Les opérations restantes
sur les champs, QR et SVD restent flottantes : le verdict global est une
qualification numérique, avec défaut du facteur explicitement traité.
Le juge contrôle chaque colonne et toutes les combinaisons réelles des six
charges, sans former leur petit Gram susceptible d'amplifier l'annulation.

## 5. Campagne à paramètres fixes

Deux variantes, trois maillages, une chauffe et trois mesures chacune :
24 processus frais, CPU 8, un fil, mêmes Python 3.14.7 / NumPy 2.5.3 /
SciPy 1.18.1. Quatre blocs Krylov, une direction conservée, profondeur
Bernstein 8, intervalle 0–40 Hz, 257 fréquences, six charges. La LU est
équilibrée et reçoit deux corrections du résidu calculé avec D et M.

Les durées incluent préparation, sélection, permutations, certificats,
contrôle, reconstruction des champs et normes physiques. Compilation,
oracles, sauvegardes et jugement indépendant sont hors chronomètre pour
les deux variantes. L'import des bibliothèques n'est pas chronométré.

| Éléments | Vinkulum préparation | Vinkulum réponses | Vinkulum total | LU total | Écart du total |
|---:|---:|---:|---:|---:|---:|
| 32 | 0,04669 s | 0,10135 s | 0,14772 s | 0,10124 s | +45,9 % |
| 128 | 0,09710 s | 0,13112 s | 0,22798 s | 0,22962 s | −0,7 % |
| 512 | 0,32627 s | 0,23305 s | 0,55881 s | 0,72911 s | −23,4 % |

Chaque colonne est la médiane de ses trois mesures ; la somme des médianes
de phases n'est pas nécessairement la médiane du total. Les étendues du
total à 512 sont [0,55536 ; 0,56043] s pour Vinkulum et
[0,70617 ; 0,73332] s pour la LU. Trois mesures ne caractérisent pas la
variabilité entre machines ou sessions. Le résultat à 128 est une quasi-égalité.

Les erreurs physiques maximales, colonnes et opérateurs confondus, sont
7,59×10⁻¹¹, 5,67×10⁻¹⁰ et 2,23×10⁻⁹ pour Vinkulum. Les plus grands
majorants relatifs retournés sont respectivement 1,91×10⁻⁹, 1,75×10⁻⁹ et
1,82×10⁻⁸ : ils restent sous 10⁻⁶. Les 24 champs et les 12 contrôles sont
admis. La confrontation numérique des majorants autorise la même marge
64ε fois la norme de référence que les campagnes précédentes ; elle ne
constitue pas une certification machine des réponses.

## 6. Sources, apport et prochaine limite mathématique

[Li et Zikatanov, 2021, prépublication révisée en août 2020](https://arxiv.org/pdf/2002.06697),
sections 3.1–3.2, relient équivalence spectrale et estimateurs bilatéraux
via des préconditionneurs et une orthogonalité de Galerkin. Ce cadre
étaye l'intérêt de métriques auxiliaires calculables. Notre transport sur
ker(Bᵀ) est démontré ici ; leur résultat elliptique ne prouve pas à lui
seul le contrôle d'un système dynamique indéfini.

[Grunert, Fehr et Haasdonk, 2020](https://onlinelibrary.wiley.com/doi/abs/10.1002/zamm.201900186)
étudient la réduction du coût de préparation d'estimateurs pour systèmes
mécaniques linéaires du second ordre. Le résumé consulté annonce un gain
d'un ordre polynomial en temps et stockage grâce à des propriétés spectrales
et des développements en série. Aucun facteur d'accélération de cet article
n'est attribué à Vinkulum ; notre preuve repose sur les formules explicites
ci-dessus et les contrôles du dépôt. Les recherches ciblées sont closes
pour ce lot : les hypothèses du transport et la provenance des mesures
sont établies ; une recherche plus large ne validerait pas davantage ces champs.

La diagonale peut donner de mauvaises constantes pour d'autres masses.
Un prolongement mathématique précis est de ne prouver la minoration de M
que sur ker(Bᵀ), et la majoration que sur les images réellement interrogées
par le contrôleur. Cela pourrait resserrer les bornes sans factorisation
globale. Ce prolongement demande ses propres certificats et mesures ; il
n'est pas implémenté ici. Les trajectoires non linéaires, contacts,
assemblages généraux et pertes d'indépendance des contraintes restent
hors du résultat. Le transfert de cette chaîne dans l'API publique est ouvert.

## 7. Conservation et reproduction

L'[archive compacte](bancs/masse-comparee-2026/manifest.json) conserve les
sources, binaires, commandes, modèles, bases, intervalles, verdicts et
journaux. Les grands NPY de champs et références sont exclus, avec leurs
SHA-256 et tailles : l'archive seule ne permet pas de rejuger les champs.
Le générateur permet de recalculer des références avec des modèles source
identiques. Un rejeu de référence n'est pas une nouvelle mesure de vitesse.

Le premier générateur a échoué au jugement, car il appelait encore le
juge diagonal. Ses fichiers et son journal sont conservés séparément.
Les tests initiaux refusés sont aussi archivés : un test exigeait à tort
une marge positive pour une base insuffisante ; un autre exigeait un défaut
strictement positif pour un Cholesky exact. Les tests corrigés vérifient
explicitement les refus et les deux types de défaut. Aucune campagne
interrompue n'a été reprise et aucun essai formel n'a été supprimé.

```sh
python ci/test_comparaison_masse.py
python ci/test_controle_masse_comparee.py
python ci/test_oracle_champs_couples.py
python ci/test_juge_masse_couplee.py
python ci/test_modeles_masse_consistante.py
python ci/test_archive_masse_comparee.py
```

Le dernier test requalifie les diagnostics conservés et recompile uniquement
le code courant connu pour rejouer trois paires de certificats : inertie
du complément et comparaison massique. Il n'exécute pas de code archivé.
Pour de nouvelles mesures, les deux scripts `ci/modeles_masse_consistante.py`
et `ci/experience_masse_comparee.py` exigent des dossiers de sortie neufs.
