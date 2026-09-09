# Facteurs d'énergie : ce que le noyau peut fournir

Le noyau possède déjà les six déformations constitutives de chaque poutre.
L'extraction ajoutée dans [tangent.rs](../src/tangent.rs) fournit leur dérivée
pondérée directement, sans assembler ni factoriser un Gram. Initialement
livrée comme brique Rust locale, elle alimente depuis la version 0.10.0
l'API publique `Noyau.facteurs_materiels_poutres()` et la réduction matérielle
par ports. L'extraction ne modifie aucune trajectoire et ne remplace pas la
tangente générale du système. Le [guide de réduction](REDUCTION_PORTS.md)
donne un exemple natif exécutable et le contrat des réponses physiques.

## Inventaire vérifié

| Composant existant | Données et calcul réellement disponibles | Conséquence pour un facteur D |
|---|---|---|
| `Poutre`, dans [lib.rs](../src/lib.rs) | Deux corps, longueur de référence, trois rigidités de translation et trois de rotation, repères matériels locaux, choix milieu/intégrée | Les six déformations et leurs poids physiques sont conservés ; aucune récupération depuis K n'est nécessaire. |
| `Poutre::deform_t`, `coefficients`, `forces_analytiques_t`, dans [tangent.rs](../src/tangent.rs) | Cinématique générique sur nombres duaux ; déformations γ,κ ; correction de flexibilité intégrée ; gradient analytique et sa tangente | Un passage de dérivées premières suffit pour extraire six lignes locales. Il réutilise les conventions exactes des efforts. |
| `Modele::facteurs_materiels_poutres`, dans [analyse.rs](../src/analyse.rs), et son binding `Noyau` | Empilement direct des blocs locaux en CSC ; z, masse, contraintes et provenance séparés | D public porte seulement l'énergie matérielle des poutres ; aucun K, Z ni multiplicateur initialisé n'est nécessaire à son extraction. |
| `Modele::linearisation_creuse`, dans [analyse.rs](../src/analyse.rs) | K,C,M,G assemblés en CSC ; précontrainte des liaisons, effets inertiels et contributions locales du champ de forces | Cette K n'est pas, en général, un Gram positif. L'API ne conserve pas les facteurs élémentaires d'origine. |
| `Superelement`, dans [lib.rs](../src/lib.rs) et [superelement.rs](../src/superelement.rs) | K dense importée, nœuds, poses relatives de référence, amortissement β ; énergie et efforts corotationnels | Aucun facteur énergétique d'origine ni matrice de masse réduite n'est conservé. Une extension doit compléter ce type existant. |
| `craig_bampton`, dans [reduction.py](../python/vinkulum/reduction.py) | K_r,M_r,T calculés dans une réduction linéaire dense | Cette fonction ne crée pas de coordonnées modales internes natives ni d'élément portant M_r couplée. |
| `Modele::energie_invariante`, dans [lib.rs](../src/lib.rs) | Refus explicite lorsque des poutres ou superéléments sont présents : leur potentiel n'est pas compté dans ce diagnostic global | L'extraction locale ajoutée ici ne rend pas ce diagnostic d'énergie générale complet. |

La masse d'analyse reste assemblée depuis les corps : bloc mI pour les
translations et RJRᵀ pour les rotations. Les contraintes G demeurent séparées
des contributions élastiques. Elles ne doivent pas être transformées en
raideurs de pénalité pour fabriquer artificiellement un facteur positif.

## Contrat Rust livré

```rust
Poutre::facteur_materiel(&self, corps: &[Corps])
    -> Result<FacteurMaterielPoutre, String>
```

Le résultat `pub(crate)` contient :

- `noeuds: [usize; 2]`, les indices physiques A et B ;
- `deformation_ponderee: SVector<f64, 6>`, le vecteur z ;
- `d: SMatrix<f64, 6, 12>`, la dérivée D en colonnes
  `(dr_A,dtheta_A,dr_B,dtheta_B)`, dans le repère monde.

Les rotations sont perturbées à gauche, comme dans l'analyse et Newton :
R(δθ)=exp([δθ]×)R. Les lignes suivent les composantes matérielles
`(γ₁,γ₂,γ₃,κ₁,κ₂,κ₃)`. Les translations physiques sont des longueurs et
les rotations des angles ; aucune normalisation en masse ni base admissible
dense n'est introduite.

Pour une poutre de longueur L, on définit

\[
 z(q)=\begin{pmatrix}
 \sqrt{LC_{N,\mathrm{eff}}}\,\gamma(q)\\
 \sqrt{LC_M}\,\kappa(q)
 \end{pmatrix},\qquad D(q)=\frac{\partial z(q\oplus\delta)}{\partial\delta}\bigg|_0.
\]

Les matrices constitutives sont diagonales. La formulation milieu utilise
les rigidités déclarées. La formulation intégrée réutilise exactement
`Poutre::coefficients`, notamment les deux associations anisotropes :

\[
 C_{N,\mathrm{eff},2}^{-1}=GA_y^{-1}+L^2/(12EI_z),\qquad
 C_{N,\mathrm{eff},3}^{-1}=GA_z^{-1}+L^2/(12EI_y).
\]

Cette correction décrit la flexion interne condensée du modèle intégré ;
elle ne réécrit pas les modules physiques stockés dans `cn`. La dérivée de
`deform_t` est évaluée avec un seul `Dual<f64,12>`. Les différences de
positions conservent la partie basse employée par les efforts natifs.
L'extraction utilise six racines scalaires, aucune Cholesky, aucun spectre
et aucun produit DᵀD.

Le vecteur z n'est jamais mis à zéro par un seuil. Les nœuds invalides,
longueurs ou rigidités non positives/non finies, poids hors domaine flottant,
valeurs ou dérivées non finies et coupure de rotation relative à π sont
refusés. Les validations ne sont pas un certificat d'arrondi.

Le type et sa méthode restent internes au crate. Ils sont désormais appelés
par l'extraction publique ; les annotations `allow(dead_code)` du prototype
initial ont été supprimées. L'intégrateur natif conserve son chemin de
forces et de tangentes existant.

## Extraction Python publique en 0.10.0

`noyau.facteurs_materiels_poutres(t=None)` renvoie un dictionnaire dont les
tableaux possèdent leurs données :

| Clé | Contenu |
|---|---|
| `d` | D global CSC, de taille `6*nombre_poutres` × `6*nombre_corps`. |
| `deformations` | z, six composantes pondérées par poutre dans le même ordre que les lignes de D. |
| `masse` | Masse spatiale CSC de tous les corps ; chaque bloc RJRᵀ est symétrisé coefficient par coefficient dans cette nouvelle API. |
| `contraintes`, `phi` | G CSC et valeurs des contraintes à la date demandée. |
| `poutres` | Liste `(nom, indice_A, indice_B)` dans l'ordre d'empilement des blocs. |
| `t` | Date d'évaluation de G et phi, date du noyau par défaut. |
| `infos_domaine` | Diagnostics de contributions et d'états hors du périmètre matériel. Cette liste n'est pas un test d'équilibre ni de positivité. |
| `etat` | Copie descriptive des corps, vitesses, multiplicateurs, contraintes, charges et diagnostics/états internes de contact. Ce n'est pas une sauvegarde de redémarrage. |

Chaque matrice CSC est le tuple
`(n_lignes, n_colonnes, indptr, indices, donnees)`, compatible avec
`scipy.sparse.csc_matrix((donnees, indices, indptr), shape=(n_lignes, n_colonnes))`.
Les colonnes globales sont `(tx,ty,tz,rx,ry,rz)` par corps. Les corps
encastrés sont conservés : leur suppression appartient à l'adaptateur ou
à l'appelant, pas à l'extraction brute.

La nouvelle masse applique aux paires hors diagonale le même coefficient
`0.5*J_s[i,j] + 0.5*J_s[j,i]`, afin d'obtenir une symétrie exacte des données
flottantes retournées. L'ancienne `k_c_m_g_creux` conserve ses arrondis
historiques. Les lignes de contact complémentaire sont nulles dans G et
phi, comme dans l'analyse statique ; le masque d'activité dynamique n'est
pas appliqué. `infos_domaine` et `etat` permettent de reconnaître ces cas.

L'appel ne calcule aucune réaction et n'avance ni poses, ni vitesses,
ni états de contact. Un `t` différent évalue seulement les lois de
contraintes à cette date. L'extraction reste partielle lorsque d'autres
éléments sont présents ; elle ne les incorpore pas silencieusement dans D.
Les données retournées peuvent être utilisées pour définir explicitement
un modèle matériel linéaire `K=DᵀD`, y compris à une pose précontrainte,
sans confondre ce choix avec la linéarisation complète du noyau.

## Identités prouvées et limite de la précontrainte

Dans la branche différentiable de la cinématique native,

\[
 U(q)=\tfrac12 z(q)^Tz(q),\qquad f_{\mathrm{poutre}}(q)=-D(q)^Tz(q).
\]

Ces deux identités concernent l'énergie constitutive de cette poutre et son
gradient. Le produit DᵀD décrit la variation quadratique des déformations
linéarisées. À déformation nulle, il égale la raideur locale élastique :

\[
 z(q)=0\quad\Longrightarrow\quad
 K_{\mathrm{poutre}}(q)=-\partial f_{\mathrm{poutre}}/\partial\delta=D(q)^TD(q).
\]

Le signe est essentiel : `Poutre::tangente` renvoie la dérivée des **forces**,
et l'assemblage d'analyse la prend avec un signe moins.

Dans une carte de configuration fixe, la Hessienne de l'énergie comprend
en général un second terme :

\[
 \nabla^2U=D^TD+\sum_a z_a\nabla^2z_a.
\]

Il n'est pas licite de supprimer ce terme en déclarant la tangente égale à
DᵀD. La dérivée native du champ de moments sous perturbations successives
à gauche comporte aussi le transport dû à ces rotations non commutatives.
Elle peut être non symétrique hors équilibre, comme les tests préexistants
de `tangent.rs` le vérifient déjà. Un Gram positif ne peut pas reproduire
cette matrice générale.

De même, la K complète de `k_c_m_g_creux` ajoute les réactions
`∂(Gᵀλ)/∂q`, les contributions des couples et contacts, l'aérodynamique et
les effets gyroscopiques pertinents. Même si toutes les poutres sont au
repos, le facteur des seules poutres n'est pas automatiquement un facteur
de cette K complète. Les modes rigides, contraintes et contributions
manquantes doivent rester explicitement représentés.

La première utilisation raisonnable est donc l'analyse linéaire d'une
sous-structure élastique non précontrainte, ou l'étude d'un préconditionneur
matériel dont les termes omis sont traités séparément. Cette livraison
n'ajoute ni intégrateur non linéaire réduit ni assertion de stabilité des
trajectoires.

## Ne pas créer un deuxième superélément

Le `Superelement` existant calcule déjà les déplacements relatifs dans le
repère du premier nœud, puis les efforts par travail virtuel. La première
coordonnée relative vaut zéro ; les réactions sur ce nœud assurent les
résultantes. La tangente différencie cette cinématique corotationnelle.
Il existe donc déjà des termes provenant de la variation de ce repère.
Cela ne fournit pas pour autant un modèle supplémentaire de raideur d'un
continuum précontraint ou de sa réduction sous centrifugation.

L'import actuel exige K finie et symétrique, mais ne vérifie pas K positive
semi-définie. Toute K acceptée n'admet donc pas nécessairement un facteur
réel D tel que K=DᵀD. Il faut préserver ce contrat au lieu d'imposer
silencieusement la positivité à l'ancienne API.

Pour les superéléments effectivement énergétiques, la suite cohérente serait
de conserver **en option** un facteur fourni par le producteur EF, avec sa
correspondance des coordonnées de référence, dans la représentation
existante. Ce facteur pourrait être utilisé avant que les petits termes
d'énergie ne soient perdus par l'assemblage flottant de K. Refaire une
Cholesky de K déjà arrondie ne restitue ni ces termes ni leur origine
physique. Le contre-exemple d'annulation de
[test_ports_releves.py](../ci/test_ports_releves.py) matérialise ce problème.

Le stockage natif actuel du superélément ne comprend ni M_r couplée entre
nœuds, ni coordonnées internes modales, ni opérateur fréquentiel S(ω).
Une matrice issue de la réduction par ports ne peut donc pas être branchée
comme une simple K constante en prétendant transporter également sa
dynamique. L'ajout d'un facteur énergétique et l'intégration d'un modèle
réduit dynamique sont deux travaux distincts.

## Vérification et limites de l'intégration

Quatre tests Rust nommés `facteur_materiel_*` contrôlent :

1. Les six déformations élémentaires d'une poutre droite anisotrope, pour
   les deux formulations, et l'égalité DᵀD=-tangente au repos.
2. U=½||z||², le travail virtuel f=-Dᵀz, l'objectivité sous transformation
   rigide et l'annulation des directions de translation/rotation rigides.
3. Un contre-exemple précontraint où le gradient reste correct mais la
   tangente diffère substantiellement du Gram matériel et n'est pas symétrique.
4. Les refus de paramètres invalides, de valeurs non finies et de la
   coupure de rotation à π.

Cinq tests Rust `facteurs_poutres_*` vérifient en plus l'assemblage global,
la masse et G, la conservation de l'état, l'absence de factorisation ou
d'initialisation des multiplicateurs, la précontrainte, les contacts à
ligne nulle et les entrées invalides. Les tests publics de
[test_reduction_ports.py](../python/vinkulum/test_reduction_ports.py)
exercent notamment l'API native sous rotation de repère, les branches 3D,
le travail des forces, les copies et l'exclusion des superéléments.

Les contrôles ciblés sont `cargo test --release facteur_materiel`,
`cargo test --release facteurs_poutres`, `cargo fmt --check` et
`cargo clippy --release --all-targets -- -D warnings`. Depuis une roue
installée avec l'extra `reduction`, le contrat Python s'exécute avec
`python -m unittest vinkulum.test_reduction_ports`.

La version 0.10.0 raccorde cette extraction à
`ReductionMaterielle.depuis_noyau`. Cet adaptateur impose une sélection
de coordonnées admissible pour G ; la masse, les ports, leur métrique
et le minorant intérieur restent distincts. Son certificat automatique
concerne seulement la minoration spectrale des valeurs binary64 de D et
M interprétées exactement. Les bornes de champ sont évaluées en doubles,
sans certificat machine complet. Aucune trajectoire non linéaire réduite,
aucun remplacement du `Superelement` existant et aucun gain universel de
performance ne sont revendiqués par cette livraison.
