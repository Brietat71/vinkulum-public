# Composer les réductions contraintes en un système mécanique

8 septembre 2026. Prototype de recherche, hors API publique 0.11.0.

`AssemblageComplements` relie plusieurs réductions rapides avec leurs
contrôles d'erreur. Les interfaces sont partagées selon les raccordements
physiques ; les directions intérieures retenues restent propres à chaque
pièce. Le contrôle porte sur le système complet. Une pièce isolée singulière
ne provoque plus à elle seule un refus de l'assemblage.

Sur trois branches flexibles avec masses consistantes et jonction rigide
linéarisée, les 24 calculs et les 12 contrôles passent à 10⁻⁶. À 512 éléments
par branche, le total médian est de 1,665 s contre 2,288 s pour la LU corrigée.
Ce résultat étend la couverture du prototype ; il ne classe pas les moteurs
multicorps complets et ne constitue pas une livraison d'API.

## 1. Coordonnées, travail et masses

Pour la pièce j, soit Wⱼ sa normalisation de port, Xⱼ(ω) son relèvement
exact du complément, et qⱼ=(yⱼ,aⱼ) les ports normalisés et coordonnées
intérieures conservées. Soit u le déplacement des interfaces globales,
et Aⱼ l'application physique qui impose u_port,j=Aⱼu. On choisit une
métrique globale H=W⁻ᵀW⁻¹ et u=Wy.

Le vecteur global est q=(y,a₁,…,a_N). L'application Eⱼ satisfait

\[
q_j=E_jq,\qquad
(E_jq)_{\mathrm{port}}=W_j^{-1}A_jWy,
\qquad(E_jq)_{\mathrm{retenu}}=a_j.
\]

Les blocs aⱼ occupent des colonnes distinctes. Les identifier entre pièces
imposerait une liaison intérieure inexistante. Le second membre vaut
(Wᵀf,0,…,0) pour des forces physiques f appliquées aux interfaces globales.
La relation Aⱼᵀ transfère les efforts : elle conserve le travail virtuel.
Les charges intérieures ne sont pas couvertes par ce prototype.

Pour une jonction rigide linéarisée, avec décalage rⱼ et rotation de repère Rⱼ,

\[
A_j=\begin{pmatrix}R_j&-R_j[r_j]_\times\\0&R_j\end{pmatrix}.
\]

La translation du point attaché est v+θ×rⱼ. Le transfert d'effort ajoute
le moment rⱼ×F. Ces deux identités et le travail sont vérifiés exactement
sur les données dyadiques du banc.

La masse et l'énergie élastique du système sont la somme des contributions
des pièces et des termes externes déclarés, ajoutés une fois. Aucun corps
supplémentaire n'est créé automatiquement à une interface. Les intérieurs
des pièces sont privés : les couplages directs entre intérieurs de pièces
différentes ne font pas partie de ce modèle.

## 2. Transport correct des erreurs locales

Supposons les Schur exacts et approchés symétriques, avec

\[
\|\widehat S_j-S_j\|_2\le b_j.
\]

Le bloc global conservé est

\[
S=S_e+\sum_j E_j^TS_jE_j,\qquad
\widehat S=S_e+\sum_j E_j^T\widehat S_jE_j.
\]

Les congruences préservent l'ordre de Loewner. Donc

\[
-G\preceq\widehat S-S\preceq G,\qquad
G=\sum_j b_jE_j^TE_j,
\quad b=\lambda_{\max}(G)\ge\|\widehat S-S\|_2.
\]

Les défauts locaux peuvent être indéfinis. Ajouter simplement les nombres
bⱼ sans leurs applications est faux si un raccordement amplifie des
coordonnées : une contre-épreuve rationnelle produit une erreur globale
qui dépasse cette somme naïve.

### Une métrique globale qui évite une perte artificielle

Si Hⱼ=Wⱼ⁻ᵀWⱼ⁻¹ et H=ΣAⱼᵀHⱼAⱼ>0, alors

\[
\sum_j E_j^TE_j=I.
\]

Le bloc des ports est l'identité par normalisation énergétique ; chaque
mode privé apparaît exactement une fois. Il s'ensuit

\[
G\preceq (\max_j b_j)I.
\]

Le nombre de pièces ne multiplie donc pas nécessairement la borne globale.
C'est une conséquence de la métrique et de la séparation des modes privés,
pas une hypothèse d'indépendance statistique des erreurs. Le test rationnel
utilise trois pièces et une normalisation 1/3 pour vérifier cette identité
sans racine carrée approchée.

Le constructeur accepte d'autres métriques et calcule toujours G. Dans le
banc, H est assemblée puis stockée en doubles : l'identité ci-dessus n'est
pas supposée exacte pour ses arrondis, et la valeur de G est effectivement
calculée. Aucun nouveau théorème de Loewner n'est revendiqué.

## 3. Marge globale et reconstruction

Le contrôle local fournit aussi aⱼ(qⱼ), un majorant de l'action
(Ŝⱼ−Sⱼ)qⱼ. En notant r=f̂−Ŝq̂, on prend

\[
a(\widehat q)=\min\left(b\|\widehat q\|,
\sum_j\|E_j\|a_j(E_j\widehat q)\right),
\qquad
b_q={\|r\|+a(\widehat q)\over\sigma_{\min}(\widehat S)-b}.
\]

Cette expression majore l'erreur de coordonnées quand son dénominateur
est positif. La valeur singulière concerne **tout le bloc global conservé**,
y compris les modes privés. L'implémentation ne factorise ni n'inverse les
Schur des pièces isolées.

Pour une norme physique ν (masse ou déformation), soient ℓⱼ,ν(q̂ⱼ) une
borne locale du défaut de relèvement et Lⱼ,ν un majorant de ||Xⱼ||ν←2.
Le contrôleur précédent donne ces deux quantités, en incluant la réparation
des contraintes et le résidu de sa résolution auxiliaire. Alors

\[
\|u_j-\widehat u_j\|_\nu
\le B_{j,\nu}:=\ell_{j,\nu}(E_j\widehat q)
+L_{j,\nu}\|E_j\|b_q.
\]

Pour une énergie externe sur les interfaces, de norme d'extension Lₑ,ν,
la somme physique des énergies donne

\[
B_\nu=\sqrt{\sum_j B_{j,\nu}^2+L_{e,\nu}^2b_q^2}.
\]

La borne relative est Bν/(||û||ν−Bν) lorsque le dénominateur est positif.
Les normes candidates sont calculées sur les champs reconstruits et les
matrices physiques. Les petites combinaisons de charges ne sont pas évaluées
par un petit Gram projeté.

Ces démonstrations sont en arithmétique exacte, sous les hypothèses locales
K_II>0, M_II>0, masse globale semi-définie positive et bande sous le minorant
du complément. La [comparaison de masses](MASSE_COMPAREE_PREUVES.md) conserve
les couplages cinétiques. Les normalisations, applications, Schur, SVD,
résolutions et évaluations finales restent flottants : les réponses portent
`certification_machine=False`. Le certificat spectral local ne certifie pas
les arrondis de l'assemblage ni la conformité exacte des ports reconstruits.

## 4. Singularités : deux contre-exemples indispensables

Avec D=diag(2,3,1), M=I et ω=1, la pièce isolée a un déplacement terminal
libre à cette fréquence : son Schur terminal est singulier. Un ressort externe
avec facteur 2 ajoute une raideur 4. L'assemblage est régulier et une force
terminale 1 donne exactement u_port=1/4. Le service isolé refuse ; le nouveau
service et l'oracle assemblé trouvent cette réponse.

Inversement, D=diag(2,3,2), M=I et ω=1 donnent une pièce isolée régulière.
Ajouter une masse terminale 3 rend l'assemblage singulier. Il est refusé,
sans pseudo-inverse ni amortissement ajouté. Une marge non positive sur un
bloc pourtant inversible rend les bornes indisponibles ; elle ne transforme
pas le champ numérique en réponse qualifiée.

Les autres contre-épreuves confrontent les champs et les bornes à l'inverse
rationnelle du système physique complet : deux pièces différentes, masses
couplées, raccordements obliques, ports partiellement partagés, changements
d'unités, énergies externes, forces nulles et faibles. L'enrichissement d'une
pièce invalide l'assemblage déjà construit.

## 5. Trois branches physiques et référence indépendante

Le cas mesuré contient trois branches de 32, 128 ou 512 éléments chacune,
issues des masses consistantes P1 du lot précédent. Les facteurs D sont
multipliés par 1, 1,125 et 1,25 ; les masses par 1, 0,875 et 1,125. Les branches
ont donc des propriétés distinctes. Deux rotations de 90° et des bras de
levier dyadiques raccordent les extrémités à une jonction à six coordonnées.
Ses masse et ressorts externes sont explicitement conservés dans les entrées.

Il y a 564, 2 292 et 9 204 coordonnées physiques indépendantes. Après les
résolutions auxiliaires locales, le système global conserve neuf coordonnées :
six de jonction et trois modes privés. Les résolvantes Krylov locales restent
calculées ; le coût de toutes ces opérations figure dans les mesures.

`OracleDirichlet` assemble DᵀD en Decimal et condense les blocs intérieurs
physiques de chaque branche. Son Schur terminal n'est pas inversé.
`OracleAssemblage` raccorde ces Schur par les applications physiques, résout
l'équilibre global et reconstruit les champs par substitution arrière.
Aucune base candidate, direction retenue ou valeur de certificat ne lui est
fournie. Il admet donc lui aussi la singularité terminale d'une pièce isolée.
Ses pivots intérieurs par blocs doivent toutefois être inversibles : cette
restriction appartient à la référence de chaîne, pas à l'assembleur candidat.

Les calculs à 70 et 90 chiffres donnent les mêmes tableaux après conversion
en doubles sur les trois cas. Cela éprouve leur convergence, sans encadrer
rigoureusement les arrondis Decimal. Le juge vérifie les colonnes et toutes
les combinaisons réelles des six charges en masse, déformation et port. Il
emploie le facteur massique à défaut audité du lot précédent.

### Résultats à 0–40 Hz, 257 fréquences, six charges

Une chauffe et trois mesures, deux variantes, trois maillages : 24 processus
frais. Même CPU 8, un fil, Python 3.14.7, NumPy 2.5.3 et SciPy 1.18.1.
Chaque pièce reconstruit son QR, sa direction, ses certificats et ses quatre
blocs Krylov dans le processus mesuré. Le seuil du complément est 80 Hz,
la profondeur Bernstein 8, la comparaison massique [0,49 ; 1,51]. La LU
assemble le système physique complet, équilibre ses matrices et effectue
deux corrections du résidu avec les facteurs et raccordements originaux.

| Éléments par branche | Préparation contrôlée | Réponses contrôlées | Total contrôlé | Total LU corrigée | Écart total |
|---:|---:|---:|---:|---:|---:|
| 32 | 0,12954 s | 0,29061 s | 0,42094 s | 0,31752 s | +32,6 % |
| 128 | 0,28873 s | 0,36912 s | 0,65685 s | 0,72219 s | −9,0 % |
| 512 | 0,98208 s | 0,68308 s | 1,66471 s | 2,28758 s | −27,2 % |

Les colonnes sont des médianes séparées. Les étendues du total à 512 sont
[1,65139 ; 1,66643] s pour Vinkulum et [2,26992 ; 2,30818] s pour la LU.
La compilation, les imports, l'oracle, le jugement et les sauvegardes sont
hors chronomètre ; la reconstruction et les normes physiques sont comprises.
Trois mesures ne caractérisent pas la variabilité entre sessions ou machines.
Aucun nouveau ratio Exudyn, MBDyn ou Simpack n'est calculé.

Les erreurs maximales du candidat, colonnes et opérateurs confondus, sont
3,18×10⁻¹¹, 6,53×10⁻¹¹ et 3,25×10⁻¹⁰. Les plus grands majorants relatifs
sont 1,06×10⁻⁹, 2,42×10⁻⁹ et 1,95×10⁻⁸. Les 24 champs et les 12 contrôles
passent le seuil demandé de 10⁻⁶. La confrontation des majorants conserve
la marge numérique 64ε fois la norme de référence utilisée précédemment.

L'audit supplémentaire des raccordements mesure u_port,j−Aⱼu_global dans
la somme des métriques de port, par charge et sur leurs combinaisons via QR.
Sur les 24 champs enregistrés, son maximum relatif d'opérateur est
8,45×10⁻¹⁴. Cet audit est numérique et hors mesure ; il ne change ni les
champs ni les bornes et ne prouve pas une conformité machine exacte.

## 6. Antériorité et limites

[Huynh, Knezevic et Patera, ESAIM M2AN 2013](https://www.numdam.org/item/10.1051/m2an/2012022.pdf),
section 5, composent des erreurs locales de bases réduites avec une analyse
de perturbation du Schur global. Leur contexte est elliptique coercif ;
ils excluent explicitement les problèmes dynamiques de leur étude, page 215.
Cette antériorité motive la composition, sans justifier à elle seule notre
traitement des blocs conservés indéfinis. Les [preuves de rétention](RETENTION_INTERIEURE_PREUVES.md)
et les dérivations ci-dessus apportent les hypothèses nécessaires à ce lot.
Les recherches ciblées ont vérifié cette distinction ; aucun résultat
industriel de l'article n'est attribué au prototype.

Le bloc global est actuellement dense et plafonné à 2 048 coordonnées
conservées. La référence mesurée concerne une jonction à six coordonnées et
des branches en chaîne ; les tests rationnels couvrent aussi un réseau de
ports partiellement partagés. Les raccordements sont constants. Grandes
rotations, liaisons mobiles, contacts, amortissement et charges intérieures
ne sont pas inclus. La bibliothèque publique reste en 0.11.0 : intégrer ce
service et sa compilation native dans une roue est encore à réaliser.
La prochaine étape ne peut pas être un simple alias du prototype : elle doit
traiter les hypothèses d'entrée, les budgets, les statuts des certificats et
la composition avec les interfaces publiques existantes.

## 7. Archives et reproduction

L'[archive](bancs/assemblage-complements-2026/manifest.json) conserve les
entrées de chaque branche et jonction, les sources, bibliothèques compilées,
commandes, bases, preuves, diagnostics par fréquence, audit des raccordements
et journaux. Les NPY des grands champs et références sont exclus avec leurs
empreintes et tailles. L'archive seule ne permet donc pas de rejuger les
champs physiques. Le vérificateur requalifie les diagnostics et recompile
uniquement le code courant connu pour rejouer les neuf paires de certificats.
Il n'exécute ni source ni bibliothèque capturée.

```sh
python ci/test_assemblage_complements_preuves.py
python ci/test_assemblage_complements.py
python ci/test_oracle_assemblage_complements.py
python ci/test_archive_assemblage_complements.py
```

Pour recalculer références ou mesures, les scripts
`ci/references_assemblage_complements.py` et
`ci/experience_assemblage_complements.py` exigent des sorties neuves.
Le pilote à six processus et la campagne à 24 processus sont tous deux
conservés. Aucun essai n'a été repris, supprimé ou remplacé.
