# Attribution des écarts de repères de l'intégrateur énergie–moment

Les neuf attributions du candidat 0.18.0 sont admises et archivées. Le
diagnostic des différences entre repères est quantifié sur ce corpus.
La CI étendue et la qualification de la roue finale passent ; voir le
[dossier de livraison](bancs/version-0.18.0/qualification.json).

## Objet exact de la garantie

Le calcul décompose les différences entre deux **trajectoires discrètes
stockées**, pour une toupie libre à inertie constante, sur une grille commune.
Les valeurs binary64 deviennent des rationnels exacts. À chaque pas, les
identités de différence sont vérifiées exactement et chaque contribution
est encadrée avec arrondi vers l'extérieur. Le vecteur observé doit appartenir
à la somme des enclosures à chaque instant enregistré.

Cette attribution ne borne pas l'erreur par rapport à la solution continue.
Elle ne certifie ni toutes les trajectoires voisines ni les autres
intégrateurs, contraintes ou contacts. Les contributions dépendent de la
décomposition sécante déclarée ci-dessous : ce ne sont pas des expériences
contrefactuelles où une source aurait été supprimée. Leurs normes peuvent
se compenser et ne définissent pas des pourcentages causaux additifs.

## Traçabilité du calcul natif

`Noyau._trace_em(t_end, h)` exécute une copie du modèle et expose les moments
matériels effectivement portés par l'intégrateur, les rotations et les
vitesses spatiales, au point initial et à chaque pas. Les inerties exportées
sont celles réellement stockées, après la symétrisation du constructeur.
La méthode est un instrument de diagnostic privé, sans certificat implicite.
Le modèle source reste inchangé, y compris lorsqu'un refus est levé.

L'instrumentation ne change aucune équation ni opération flottante du pas.
La qualification confronte ses sorties aux sorties ordinaires de la même
roue et à une roue publiée 0.17.0, composante binary64 par composante, y
compris les bits de signe. Les traces complètes sont produites par
`ci/trace_moments_em.py` : quatre repères, trois pas (0,001 ; 0,0005 ;
0,00025 s), horizon 20 s. La structure, la finitude, l'inertie symétrique
définie positive, l'horizon et chaque date sont contrôlés avant attribution.
La grille suit le contrat `numerique::fin_pas`, dont le dernier pas peut
absorber un reliquat jusqu'à 1,5 fois le pas demandé.

## Identité de différence des moments

Pour la trace canonique, soit $f_A(x)=x\times J_A^{-1}x$. Pour l'autre
trace, $p_B$ est le moment matériel natif et $C$ la matrice matérielle
stockée ; on pose $y=Cp_B$. Aucune orthogonalité exacte de $C$ n'est supposée.
Les barres désignent les moyennes des deux extrémités du pas réel $h$.

Les défauts natifs sont évalués exactement :

$$d_A=x_1-x_0-hf_A(\bar x),\qquad
  d_B=p_{B1}-p_{B0}-h f_B(\bar p_B).$$

Puisque le champ est quadratique, l'identité suivante est exacte, sans reste
de linéarisation :

$$f_A(\bar y)-f_A(\bar x)
 =Df_A((\bar x+\bar y)/2)(\bar y-\bar x).$$

En posant $B=Df_A((\bar x+\bar y)/2)$ et $\delta=y-x$, on obtient

$$(I-hB/2)\delta_1=(I+hB/2)\delta_0
 +h[Cf_B(\bar p_B)-f_A(\bar y)]+Cd_B-d_A.$$

Les trois contributions sont donc : différence initiale ; données
transformées ; défauts algébriques natifs. L'inversibilité de chaque matrice
$I-hB/2$ est contrôlée par élimination rationnelle. Un refus d'inversion
interdit de produire l'attribution.

La propagation emploie une base fondamentale proposée en Decimal à
70 chiffres puis quantifiée sur 160 bits dyadiques. Sa matrice de transition
normalisée est recalculée en rationnels exacts ; l'erreur de proposition
n'est jamais ignorée. Les enclosures sont propagées par planchers et
plafonds rationnels. Decimal ne fait pas partie des opérations de preuve.

## Rotation et vitesse spatiale

Les rotations comparées sont $R_A$ et $R_B=M^T R_{B,natif}C^T$, avec $M$
la transformation spatiale stockée. On définit

$$C_A=\operatorname{cay}(h[J_A^{-1}\bar x]_\times),\qquad
 C_B=C^{-T}\operatorname{cay}(h[J_B^{-1}\bar p_B]_\times)C^T.$$

Pour $X_A=h[J_A^{-1}\bar x]_\times$ et
$X_B=C^{-T}h[J_B^{-1}\bar p_B]_\times C^T$, l'identité de résolvante

$$C_B-C_A=\tfrac12(C_B+I)(X_B-X_A)\tfrac12(C_A+I)$$

est vérifiée exactement à chaque pas. Elle permet de séparer la contribution
des moments de celle des données transformées dans

$$\delta R_1=\delta R_0 C_A+R_{B0}(C_B-C_A)+\delta E,$$

où $\delta E$ est la différence transformée des défauts natifs de mise à
jour de rotation. Ces défauts comprennent les produits flottants et la
projection polaire effectivement exécutée. Une base proposée à droite
stabilise les enclosures ; ses transitions sont également vérifiées en
rationnels exacts.

Enfin, avec $H_B=C^{-T}J_B^{-1}C^{-1}$,

$$\delta\omega=\delta R J_A^{-1}x+R_BJ_A^{-1}\delta p
 +R_B(H_B-J_A^{-1})y+\delta d_\omega.$$

Les erreurs du post-traitement NumPy historique (produits de matrices puis
soustraction) forment une quatrième contribution explicite. Leurs valeurs
flottantes sont relues en rationnels et comparées au résultat algébrique
exact. Les normes de Frobenius des rotations et euclidiennes des vitesses
sont encadrées au moyen de racines carrées entières dirigées.

## Critères fixés avant la campagne

- Défauts locaux normalisés de moment, rotation et vitesse : au plus
  $64\,2^{-52}$. Les échelles et normes sont celles du code du vérificateur.
- Défaut d'orthogonalité en norme infinie des rotations natives et des
  transformations : au plus $64\,2^{-52}$ ; déterminants strictement positifs.
- Largeur maximale des enclosures : $10^{-20}$ pour le moment et la rotation,
  $10^{-18}$ pour la vitesse spatiale.
- Recomposition encadrée à chaque pas ; horizon complet obligatoire pour
  une attribution globale admise. `--limite` produit un diagnostic partiel
  qui ne peut pas être déclaré complet.

Les contrôles locaux d'énergie et de norme du moment vérifient aussi les
identités télescopiques exactes
$\Delta E=(J^{-1}\bar p)\cdot d$ et
$\Delta\lVert p\rVert^2=2\bar p\cdot d$. Ils ne remplacent pas la
propagation des différences entre repères.

## Contre-épreuves et reproduction

### Résultats complets

Les neuf comparaisons passent : 420 000 pas comparés. Les maxima ci-dessous portent
sur tous les pas, alors que le rapport historique n'en échantillonnait
qu'environ 200. Les maxima historiques sont aussi retrouvés à l'identique
pour les douze couples repère/pas de la 0.14.1.

| Repère | Pas, s | Rotation, Frobenius | Vitesse, rad/s |
|---|---:|---:|---:|
| spatial | 0.001 | 1.774911e-07 | 7.246260e-07 |
| materiel | 0.001 | 9.135476e-08 | 3.729618e-07 |
| combine | 0.001 | 3.120817e-08 | 1.274107e-07 |
| spatial | 0.0005 | 3.514053e-07 | 1.434612e-06 |
| materiel | 0.0005 | 6.838528e-08 | 2.791830e-07 |
| combine | 0.0005 | 9.823193e-08 | 4.010317e-07 |
| spatial | 0.00025 | 1.692101e-07 | 6.908011e-07 |
| materiel | 0.00025 | 8.145099e-08 | 3.325236e-07 |
| combine | 0.00025 | 8.101220e-08 | 3.307320e-07 |

Au maximum de vitesse du cas spatial au pas de 0,001 s, les normes des contributions sont
environ 5,76e-11 pour l'état initial, zéro pour les données matérielles
transformées, 7,246836e-7 pour les défauts algébriques propagés et 1,17e-15
pour le post-traitement. La somme vectorielle, avec ses compensations,
encadre la valeur observée à chaque pas. Pour chacun des neuf pics de vitesse, la borne inférieure de la norme
de la contribution des défauts dépasse la somme des bornes supérieures des
trois autres contributions. Cette comparaison est faite en rationnels exacts.
Les largeurs maximales des enclosures restent inférieures à 1,08e-37 pour
le moment, 7,93e-35 pour la rotation et 2,83e-33 pour la vitesse.

L’[archive complète](bancs/attribution-reperes-0.18.0/manifest.json) conserve
les douze traces, les neuf résultats, les sources de calcul et les journaux.
Le [bilan](bancs/attribution-reperes-0.18.0/bilan.json) contient des approximations
d’affichage ; les bornes normatives restent les rationnels des résultats.

Les tests couvrent les arrondis dirigés et tous les coins d'une boîte,
l'identité de résolvante sous transformation non orthogonale, une petite
perturbation dont l'attribution doit passer, de grands défauts de moment,
rotation ou vitesse dont l'attribution doit échouer, ainsi que les traces
tronquées, mal dimensionnées, non finies et à pas manquant.

```bash
python ci/test_trace_moments_em.py
python ci/test_attribution_moments_em.py
python ci/test_attribution_rotation_em.py
python ci/attribue_rotation_em.py initial-0.00100.json.gz \
  materiel-0.00100.json.gz --sortie attribution.json
```

Le rejeu rationnel complet est plus coûteux que la simulation. Un test
synthétique court ou un contrôle de hachage ne remplace pas ce rejeu. La
base de confiance comprend l'interpréteur Python, ses entiers et Fraction,
le code du vérificateur et la fidélité de la capture native. Ce calcul
n'est pas une preuve de l'implémentation du vérificateur dans Lean.


Contrôler l’archive et ses 36 pics :

```bash
python ci/archive_attribution_reperes.py \
  --verifier docs/bancs/attribution-reperes-0.18.0 --observations
python ci/test_archive_attribution_reperes.py
```

Rejouer intégralement les neuf attributions avec leurs sources archivées
(durée de l’ordre de plusieurs dizaines de minutes, selon la machine) :

```bash
python ci/archive_attribution_reperes.py \
  --rejouer docs/bancs/attribution-reperes-0.18.0
```

Ce rejeu strict exige le même environnement NumPy pour reproduire également
les observations flottantes et leurs métadonnées. Les opérations de preuve
restent rationnelles ; aucune précision de NumPy n’est supposée.
