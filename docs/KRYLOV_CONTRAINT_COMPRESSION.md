# Comprimer le contrôle résiduel sans supposer le QR exact

Programme de recherche pour le cycle suivant, 8 septembre 2026. Cette note
ne modifie aucun code de campagne. Les propriétés algébriques ci-dessous
sont démontrées en arithmétique exacte ; leur réalisation avec encadrement
des arrondis reste à construire. **Aucun gain de temps n'a été mesuré pour
cette compression.** Le [profil mesuré](KRYLOV_CONTRAINT_PROTOTYPE.md)
place le certificat Decimal devant la construction de l'enveloppe de
contrôle. Accélérer uniquement cette construction ne suffira pas à
dépasser la LU en temps total sur le cas à 512 poutres ; la section 7
quantifie cette limite.

Depuis cette mesure, le [certificat d'inertie](INERTIE_CONTRAINTE_PROTOTYPE.md)
a réduit le coût spectral, puis le [contrôle par facteurs](CONTROLE_FACTEURS_PROTOTYPE.md)
a accéléré les normes d'opérateur et les images de réparation par fréquence.
Ce dernier transfert conserve les champs physiques et n'implémente pas
la compression commune de l'enveloppe Bernstein décrite ci-dessous.
Les coûts historiques de la section 7 restent ceux du premier prototype.

Les conventions sont celles de la
[preuve du Krylov contraint](KRYLOV_CONTRAINT_PREUVES.md) et de sa
[variante anisotrope](KRYLOV_CONTRAINT_ANISOTROPIE.md). Le modèle physique
reste celui des valeurs stockées de D et M ; le sous-espace éliminé reste
le noyau du B effectivement certifié.

## 1. Une image commune à tous les résidus

Noter m le nombre de coordonnées conservées, k le nombre de directions du
Krylov, et n_c la dimension du complément. Pour les matrices préparées de
la preuve précédente, le résidu du champ réparé idéal est

\[
R(\mu)=E_0+\mu C_1+\mu^2 F(I-\mu\Theta)^{-1}D_*,
\qquad 0\le\mu\le\rho.
\]

Ici D_* est une petite matrice k×m, distincte du facteur matériel D.
Fixer un facteur dual exact H_M de taille n_c×n tel que
H_M^T H_M=S_M, où S_M est l'inverse massique contrainte. Poser

\[
U_0=H_M E_0,\quad U_1=H_M C_1,\quad U_F=H_M F,
\qquad U=[U_0,U_1,U_F].
\]

U a n_c lignes et \(\ell=2m+k\) colonnes. Tout le résidu transformé
appartient à son image :

\[
H_MR(\mu)=U\,C(\mu),\qquad
C(\mu)=
\begin{bmatrix}
I_m\\ \mu I_m\\ \mu^2(I-\mu\Theta)^{-1}D_*
\end{bmatrix}.
\]

Supposons d'abord un QR mince **exact**
\(U=QT\), avec \(Q^TQ=I\). On peut garder
\(a=\min(n_c,\ell)\) colonnes de Q sans décision de rang ; si U est de
rang inférieur, les zéros du facteur T restent présents. Aucune direction
n'est supprimée au moyen d'un seuil numérique. Partitionner
\(T=[T_0,T_1,T_F]\) suivant les mêmes colonnes. Alors

\[
H_MR(\mu)=Q\,
\underbrace{\left[T_0+\mu T_1+
\mu^2T_F(I-\mu\Theta)^{-1}D_*\right]}_{\mathcal R(\mu)},
\]

\[
\boxed{\|H_MR(\mu)\|_2=\|\mathcal R(\mu)\|_2,\qquad
\|H_MR(\mu)v\|_2=\|\mathcal R(\mu)v\|_2.}
\]

Toutes les normes portent désormais sur des matrices à a lignes. Un QR
avec pivotage est également admissible, mais sa permutation doit être
réinjectée dans T **avant** cette partition en trois blocs. Oublier cette
permutation changerait les colonnes et donc le résidu physique.

Cette opération comprime l'image des résidus. Elle ne construit ni une
base globale Z du complément, ni une nouvelle réduction physique du modèle.
Les bases physiques nécessaires à la reconstruction des champs subsistent.

## 2. Bernstein et actions avec compensations

Soit \(t\ge\|\Theta\|_2\), avec \(\gamma=1-\rho t>0\).
Pour une profondeur h≥2, former les coefficients en puissances

\[
A_0=T_0,\quad A_1=T_1,\quad
A_j=T_F\Theta^{j-2}D_*\quad(2\le j\le h-1).
\]

Les coefficients Bernstein sur [0,ρ] sont

\[
\mathcal B_{h,j}=\sum_{i=0}^{j}
\frac{\binom ji}{\binom{h-1}i}\rho^i A_i,
\qquad 0\le j\le h-1.
\]

La preuve de convexité et la queue de résolvante donnent

\[
b_h=\max_j\|\mathcal B_{h,j}\|_2+
\frac{\rho^h\|T_F\Theta^{h-2}\|_2\,\|D_*\|_2}{\gamma},
\qquad \sup_\mu\|\mathcal R(\mu)\|_2\le b_h.
\]

Pour toute direction v des coordonnées conservées, la variante suivante
conserve les compensations entre ses composantes :

\[
\boxed{b_h(v)=\max_j\|\mathcal B_{h,j}v\|_2+
\frac{\rho^h\|T_F\Theta^{h-2}\|_2\,\|D_*v\|_2}{\gamma}.}
\]

Il faut appliquer chaque petite matrice à v **avant** la norme. La somme
des normes de colonnes pondérées par |v_j| reste une autre majoration,
mais peut perdre ces compensations. Un Gram préassemblé suivi du calcul
de \(v^T Gv\) peut également perdre une petite réponse en flottants.

Exemple exact : \(\mathcal R=(1,-1)\), sans terme fréquentiel. Pour
\(v=(1,1)^T\), l'action est nulle ; la somme des normes de colonnes vaut 2.
La représentation comprimée conserve zéro en arithmétique exacte. Cela
n'autorise pas à annoncer zéro si les défauts de compression ou les arrondis
restants sont non nuls.

## 3. Deux défauts à distinguer

Un QR flottant ne fournit pas l'identité exacte précédente. Noter
\(\widehat U\) les coefficients effectivement calculés, et
\(\widetilde Q,\widetilde T\) les facteurs stockés, interprétés comme des
objets mathématiques précis. Décomposer

\[
\underbrace{U-\widetilde Q\widetilde T}_{\Delta}
=\underbrace{U-\widehat U}_{\Delta_{\rm entree}}
+\underbrace{\widehat U-\widetilde Q\widetilde T}_{\Delta_{\rm QR}}.
\]

Le premier défaut comprend la formation des résidus à partir de D et M,
la réparation des contraintes, les résolutions massiques et les projections
duales nécessaires à U. Le second concerne uniquement la compression de
la matrice calculée. **Contrôler le second ne contrôle pas le premier.**
Un QR parfait d'un résidu préalablement arrondi à zéro ne prouve pas que
le résidu du modèle D original soit nul.

Supposer établis des nombres

\[
\|\widetilde Q\|_2\le\kappa_Q,\qquad
\|\Delta_0\|_2\le\varepsilon_0,\quad
\|\Delta_1\|_2\le\varepsilon_1,\quad
\|\Delta_F\|_2\le\varepsilon_F,
\]

où \(\Delta=[\Delta_0,\Delta_1,\Delta_F]\) suit la partition de U.
Une manière suffisante est d'additionner, par bloc, les majorants établis
des défauts d'entrée et du QR. Aucune troncature des petites valeurs
singulières de U n'est nécessaire. Si une troncature est ajoutée plus tard,
son défaut doit entrer explicitement dans ces ε.

Avec les blocs de \(\widetilde T\), refaire les petits calculs de la
section 2 et obtenir \(\widetilde b_h\) et \(\widetilde b_h(v)\). Alors

\[
\boxed{\delta_h=\kappa_Q\widetilde b_h+
\varepsilon_0+\rho\varepsilon_1+
\frac{\rho^2\varepsilon_F\|D_*\|_2}{\gamma}}
\]

majore uniformément \(\|H_MR(\mu)\|_2\). De même,

\[
\boxed{\eta_h(v)=\kappa_Q\widetilde b_h(v)+
(\varepsilon_0+\rho\varepsilon_1)\|v\|_2+
\frac{\rho^2\varepsilon_F\|D_*v\|_2}{\gamma}}
\]

majore uniformément \(\|H_MR(\mu)v\|_2\).

**Preuve.** L'identité
\(H_MR=\widetilde Q\widetilde{\mathcal R}+\Delta C\) est exacte.
Appliquer la sous-multiplicativité au premier terme, puis la somme
triangulaire aux trois blocs de \(\Delta C\), avec
\(\|(I-\mu\Theta)^{-1}\|_2\le1/\gamma\). Les normes
\(\|D_*v\|\) peuvent rester directionnelles dans le dernier terme.

On peut retenir simultanément

\[
\delta=\min_h\delta_h,\qquad
\eta(v)=\min\left(\delta\|v\|_2,\ \min_h\eta_h(v),\
\sum_j\delta_j|v_j|\right),
\]

si les δ_j sont eux aussi des majorants établis. Le choix de h peut
dépendre de v : chaque candidat est uniforme en μ et valable pour tous v.
Il reste donc permis de prendre v égal à la réponse calculée à une fréquence.

Les ε introduisent un plancher d'erreur. Lorsqu'une compensation rend
l'action comprimée très petite, ce plancher doit subsister. Par exemple,
\(U=(1,-1+2^{-30})\), \(\widetilde Q=1\),
\(\widetilde T=(1,-1)\) et \(v=(1,1)^T\) donnent une action comprimée
nulle mais une action réelle \(2^{-30}\). Ce cas exact impose le terme
additif ; il interdit de qualifier une compression uniquement par une
erreur relative à l'action comprimée.

## 4. Ne pas présumer l'orthogonalité des facteurs stockés

Si l'on établit
\(\|\widetilde Q^T\widetilde Q-I\|_2\le\varepsilon_Q\),
alors \(\kappa_Q=\sqrt{1+\varepsilon_Q}\) convient. Une valeur
\(\varepsilon_Q<1\) fournit aussi une borne inférieure sur les valeurs
singulières, mais n'est pas nécessaire à la seule majoration ci-dessus.
La racine et la norme doivent elles-mêmes être encadrées pour constituer
un certificat machine.

Exemple exact :
\(U=\widetilde Q=\operatorname{diag}(9/8,1)\),
\(\widetilde T=I\). Le défaut de factorisation est nul, mais remplacer
\(\|Uv\|\) par \(\|\widetilde Tv\|\) sous-estime l'action sur e1.
Ici \(\varepsilon_Q=17/64\) et \(\kappa_Q=9/8\).

Une autre voie utilise directement les réflecteurs stockés du QR, sans
former de Q dense. Définir exactement

\[
H_j=I-\tau_jv_jv_j^T,\qquad \nu_j=v_j^Tv_j,
\qquad \widetilde Q=H_1\cdots H_a J_a,
\]

où J_a injecte les a premières coordonnées, et où le 1 implicite de chaque
vecteur Householder fait partie de v_j. Chaque réflecteur réel est symétrique
et possède les valeurs propres 1 et \(1-\tau_j\nu_j\). Par conséquent,

\[
\boxed{\|\widetilde Q\|_2\le
\kappa_{\rm ref}:=
\prod_{j=1}^{a}\max\left(1,|1-\tau_j\nu_j|\right).}
\]

Un majorant dirigé du produit κ_ref fournit une valeur utilisable de κ_Q.
Cela ne suppose ni \(\tau_j\nu_j=2\), ni une application LAPACK exacte.
La matrice mathématique définie par les réflecteurs est distincte de
l'arrondi d'une application `ormqr` ; ce dernier doit être inclus dans
l'audit du défaut QR lorsqu'il sert à le calculer.

Pour une variante fondée sur le défaut d'orthogonalité, chaque
\(\epsilon_j=|(1-\tau_j\nu_j)^2-1|\) majore
\(\|H_j^TH_j-I\|_2\). Si tous les ε_j sont inférieurs à 1,

\[
\varepsilon_Q\le\max\left(
\prod_j(1+\epsilon_j)-1,\quad
1-\prod_j(1-\epsilon_j)\right)
\]

convient également. Les sommes de carrés définissant ν_j peuvent être
encadrées sans construire le Gram global de Q.

## 5. Auditer le défaut sans quitter le modèle d'entrée

Une réalisation concrète doit déclarer ce qui est prouvé à chaque étage :

1. **Données physiques.** E0, C1 et F visent D original, M original et le
   même B que le certificat spectral. Ils ne deviennent pas des résidus
   du modèle Q à cause du facteur utilisé pour construire le Krylov.
2. **Préparation duale.** Établir un encadrement des erreurs d'assemblage,
   de réparation et de transformation massique donnant les contributions
   à \(\Delta_{\rm entree}\). Tant que cet encadrement manque, la compression
   peut seulement préserver le statut non certifié du contrôleur actuel.
3. **Factorisation comprimée.** Contrôler
   \(\widehat U-\widetilde Q\widetilde T\) avec les valeurs stockées,
   la permutation et les réflecteurs effectivement retenus. Un résidu
   calculé en doubles sans majorant de son propre arrondi est un diagnostic.
4. **Petits calculs.** Encadrer la formation des coefficients Bernstein,
   leurs normes, les produits avec v et la borne t de Θ. Réduire leur
   dimension facilite cette étape ; cela ne la rend pas automatiquement exacte.

Pour le troisième étage, une option directe est d'appliquer les réflecteurs
à T avec arithmétique dirigée, puis de soustraire \(\widehat U\) avec le
même encadrement. Une autre option est un résidu flottant accompagné d'un
majorant absolu démontré des produits et soustractions utilisés. Le simple
fait qu'un QR soit réputé stable ne fournit pas les constantes nécessaires
pour une bibliothèque et un domaine flottant donnés.

Une borne de Frobenius des intervalles du défaut suffit à majorer sa norme
d'opérateur. Des bornes par bloc, par colonne ou des facteurs du défaut
peuvent être moins pessimistes. Elles doivent conserver leurs erreurs
d'évaluation ; recomprimer un défaut ne l'efface pas.

Les petites matrices doivent aussi être cohérentes avec l'identité
résiduelle. En particulier, D_* désigne dans cette identité le produit mathématique
\(\Theta a_0+a_1\), et non automatiquement son arrondi stocké. Si les
formules comprimées utilisent \(\widehat D_*\), le défaut
\(D_*-\widehat D_*\) n'est pas contenu dans le seul défaut du QR de U.
La résolution approchée de Y et la formation de D_* sont deux erreurs
différentes.

Plus généralement, si la petite résolvante comprimée emploie
\(\widehat\Theta,\widehat D_*\), supposer établis

\[
\|\Theta-\widehat\Theta\|\le\epsilon_\Theta,\quad
\|D_*-\widehat D_*\|\le\epsilon_D,\quad
\gamma=1-\rho t>0,\quad
\widehat\gamma=1-\rho\widehat t>0.
\]

Avec \(t\ge\|\Theta\|\), \(\widehat t\ge\|\widehat\Theta\|\) et
\(N_F\ge\|U_F\|\), la différence additionnelle de résidu est majorée par

\[
\epsilon_{\rm petit}=\rho^2 N_F\left(
\frac{\epsilon_D}{\gamma}+
\frac{\rho\epsilon_\Theta\|\widehat D_*\|}
     {\gamma\widehat\gamma}\right).
\]

Cela suit de
\((I-\mu\Theta)^{-1}-(I-\mu\widehat\Theta)^{-1}
=\mu(I-\mu\Theta)^{-1}(\Theta-\widehat\Theta)
(I-\mu\widehat\Theta)^{-1}\).
La variante directionnelle remplace ε_D par un majorant de
\(\|(D_*-\widehat D_*)v\|\), et
\(\|\widehat D_*\|\) par \(\|\widehat D_*v\|\).
Les bornes de compression de la section 3 s'appliquent alors à la
résolvante avec chapeaux, avec \(\widehat\gamma\), puis ce défaut
additionnel est ajouté. On peut prendre
\(N_F=\kappa_Q\|\widetilde T_F\|+\varepsilon_F\).
Si Θ est conservée comme même objet mathématique des deux côtés,
ε_Θ vaut zéro ; cela n'annule pas automatiquement ε_D.

Si la transformation duale emploie un facteur massique approché
\(M_Q=R_M^TR_M\), une voie supplémentaire consiste à démontrer

\[
M\succeq(1-\eta_M)M_Q\quad\text{sur }\ker B^T,
\qquad 0\le\eta_M<1.
\]

La monotonie de l'inverse sur ce sous-espace donne alors

\[
S_M\preceq\frac{S_{M_Q}}{1-\eta_M}.
\]

Un contrôle résiduel établi dans la métrique duale exacte de M_Q peut donc
être multiplié par \((1-\eta_M)^{-1/2}\) pour viser M original.
Il faut encore encadrer les applications du projecteur contraint de M_Q ;
un contrôle du seul facteur de masse ne prouve pas ce projecteur exact.
Les petits blocs massiques actuels rendent cette voie envisageable,
sans démontrer ici son coût ni son efficacité pratique.

## 6. Raccord aux champs et à l'action du Schur

Une fois δ et η(v) établis **dans la métrique originale**, reprendre
\(\alpha_* = \lambda_c-\Omega^2>0\). Pour le champ réparé idéal,
noter F_M son erreur pondérée par la racine de la masse complète et F_D
son erreur après application du facteur matériel D. Alors

\[
\|F_Mv\|\le\eta(v)/\alpha_*,\qquad
\|F_Dv\|\le\sqrt{\lambda_c}\,\eta(v)/\alpha_*.
\]

Pour son défaut de Schur positif A,

\[
\boxed{\|Av\|\le\delta\,\eta(v)/\alpha_*.}
\]

Le facteur global δ reste nécessaire. Une action résiduelle faible ne
permet pas de remplacer ce produit par \(\eta(v)^2/\alpha_*\).
Les défauts de réparation et du petit solve s'ajoutent comme dans la note
anisotrope. La marge du Schur demeure globale, de même que les constantes
d'extension multipliant l'erreur de coordonnées. La compression ne change
ni ces exigences ni les résonances physiques du système complet.

Une norme massique complète PSD, avec des ports sans masse, reste admise.
Seule la masse intérieure requise par l'inverse duale doit être SPD.

## 7. Programme conditionné par le profil

Les étapes suivantes sont proposées après achèvement et archivage du lot
en cours ; elles ne justifient aucune modification de ses sources gelées.

1. Traiter en priorité le certificat spectral : le profil final confirme
   qu'il domine la préparation. Étudier son remplacement par un
   encadrement creux du seuil, avec contrôle des arrondis, puis mesurer
   le service complet. La compression décrite ici demeure une piste
   distincte pour réduire le coût des normes.
2. Éliminer d'abord les normes identiques recalculées dans la boucle
   Bernstein : norme de la puissance FΘ^j utilisée deux fois, normes
   constantes de D_*. Revalider sans changer le contrat mathématique.
3. Prototyper le QR conjoint et comparer, sur les mêmes coefficients,
   toutes les normes, actions et majorants avec et sans compression.
   La dimension a peut être proche de n_c sur un petit problème ; dans
   ce cas, aucun bénéfice de dimension n'est promis.
4. Évaluer séparément le coût du contrôle des défauts. Une validation
   dirigée plus coûteuse que les SVD évitées peut rendre cette compression
   impropre à l'objectif de temps total, même si l'identité est correcte.
5. Si les résultats le justifient, étudier ensuite les actions avec
   compensations et le maintien de H=CΔ sous forme de faible rang, puis
   une nouvelle campagne complète avec préparation et services identiques.

Pour \(n_c\gg\ell\), un QR conjoint dense coûte classiquement
O(n_c ℓ²) et remplace les grandes normes répétées par des opérations à
a lignes. Il exige temporairement O(n_c ℓ) coefficients ; conserver
ensuite T demande O(aℓ). Ces comptes concernent le contrôle, pas les
facteurs creux ni les bases de reconstruction. Le coût d'encadrement
des défauts est additionnel et dépend de la méthode choisie.

Une borne simple évite de surinterpréter ce programme : si une phase
représente la fraction f du temps total, la supprimer entièrement ne
donne au mieux qu'un facteur \(1/(1-f)\), toutes les autres phases étant
inchangées. Sur le cas à 512 poutres, les médianes mesurées sont
**1,3026 s** pour le service contrôlé, contre **0,7254 s** pour la LU.
La construction de l'enveloppe de contrôle prend **0,0832 s**, tandis que
le certificat prend **0,6848 s**. Supprimer entièrement cette construction
laisserait un ordre de grandeur de **1,2194 s** : environ 6,4 % de temps
en moins, soit un facteur maximal d'environ 1,068 pour cette seule phase.
Le gain reste insuffisant face à la LU, même avec une compression gratuite.

Ce calcul optimiste utilise des médianes de phases et de totaux ; leur
soustraction n'est pas la médiane d'un temps contrefactuel par essai.
Les 0,0832 s concernent seulement la **préparation** du contrôleur.
Ils ne comprennent pas les contrôles par fréquence inclus dans les
0,3817 s de réponses et normes. Une éventuelle accélération de ces
contrôles doit être évaluée séparément. Aucun de ces temps ne mesure
une implémentation de la compression. Le certificat représente déjà
environ 53 % du total : réduire aussi ce poste reste la priorité du
prochain cycle, en conservant sa garantie mathématique.

Les identités exactes de cette note constituent des propriétés de
composition. La robustesse flottante, le coût du contrôle des défauts et
l'avantage mesuré sur les solveurs de référence restent des questions
expérimentales ouvertes.
