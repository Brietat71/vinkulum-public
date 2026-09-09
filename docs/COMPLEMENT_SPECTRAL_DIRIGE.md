# Complément spectral dirigé et composition à coordonnées retenues

Expérience du 8 septembre 2026, sur la roue figée 0.11.0. Une direction
intérieure conservée suffit à certifier la coercivité du complément sur
0–40 Hz pour les trois consoles natives testées. Le témoin mécanique
respecte le seuil relatif 10⁻⁶ aux 257 fréquences, pour les six charges
et toutes leurs combinaisons réelles. Il factorise encore un KKT complet
à chaque fréquence : la réduction rapide reste à construire.

Le certificat spectral porte sur les coefficients binary64 interprétés
exactement. Les champs et leur contrôle contre une référence Decimal
restent des calculs numériques sans encadrement complet des arrondis.
Cette expérience n'étend pas l'API `ReductionMaterielle` de la 0.11.0.

## Du théorème au calcul dirigé

Avec \(K_D=D_i^TD_i\), la factorisation approchée définit
\(K_Q=E^{-1}R^TRE^{-1}\). Les contraintes sont les valeurs **stockées**
de \(B\), même lorsqu'elles proviennent d'un produit flottant \(M\Phi\).
Le sous-espace éliminé est \(N=\ker B^T\). Posons

\[
U=K_Q^{-1}B,\quad H=B^TU,\quad J=U^TMU,\qquad
\tau_c=\operatorname{tr}(MK_Q^{-1})-\operatorname{tr}(H^{-1}J).
\]

L'[identité de trace contrainte](TRACE_COMPLEMENT_SPECTRAL_PREUVES.md)
donne, pour toute base \(Z\) de \(N\),

\[
\tau_c=\operatorname{tr}\big[(Z^TK_QZ)^{-1}Z^TMZ\big]>0,
\qquad
\lambda_{\min}(K_Q|_N,M|_N)\ge\tau_c^{-1}.
\]

Le [certificat expérimental](../ci/trace_complement_dirigee.py) procède ainsi :

1. Refaire une inverse sélectionnée à partir des entrées, puis certifier
   la trace totale, la positivité de la masse et
   \(K_D\succeq(1-\bar\eta)K_Q\), avec \(\bar\eta<1\).
2. Encadrer les résolutions triangulaires donnant \(U\), puis \(H,J\).
   Une LDL dirigée de \(H\) exige des pivots strictement positifs : le
   rang de \(B\) est démontré, sans seuil modal arbitraire.
3. Encadrer \(\chi=\operatorname{tr}(H^{-1}J)\), puis soustraire les
   intervalles. Une borne inférieure négative de \(\tau_c\) reste
   admissible : sa positivité exacte découle déjà des hypothèses.
4. Renvoyer \(\underline\lambda_c\le(1-\bar\eta)/\overline\tau_c\),
   avec division et conversion binary64 vers le bas. Le majorant
   \(\overline\tau_c\) doit être strictement positif.

Les opérations Decimal utilisent des contextes locaux dirigés. Le budget
couvre le certificat initial et la correction ; les limites du
prétraitement flottant et du stockage rectangulaire sont explicites.
Le coût supplémentaire porte sur des tableaux \(n\times r\) et des
matrices \(r\times r\), sans inverse complète ni base dense du complément.
Le remplissage des facteurs dépend du graphe : aucune complexité linéaire
générale n'est établie.

Les contre-épreuves couvrent les contraintes obliques, les masses
couplées, les annulations, les rangs indémontrables et les changements
d'échelle. Pour \(K=\operatorname{diag}(2^{-100},1,4,9)\), \(M=I\),
\(B=e_1\), les traces flottantes s'annulent alors que la trace restante
vaut exactement \(49/36\). Pour deux contraintes séparées par \(2^{-52}\),
16 chiffres refusent le rang, 40 peuvent donner une borne sûre mais
inutilisable, et 80 retrouvent ici une borne précise. Un certificat
mathématiquement valide n'implique donc pas l'acceptation de la bande.

## Composition mécanique

Retenir \(r\) coordonnées avec les \(p\) ports change le contrat de
condensation. Le [carnet de composition](RETENTION_INTERIEURE_PREUVES.md)
dérive le relèvement \(T\), l'inverse contrainte \(R_z\) et

\[
R_z=(I-zSM)^{-1}S,\qquad
S_c=T^T\mathcal A(z)T-G^TR_zG,\qquad
G=E_i^T\mathcal A(z)T.
\]

L'inverse contrainte existe sous la borne du complément, même si
l'intérieur complet est singulier. Tous les couplages de masse et de
raideur restent présents. La marge globale doit porter sur les
\(p+r\) coordonnées : les véritables résonances globales persistent.
Les tests rationnels distinguent explicitement un pôle intérieur
traversable et un pôle global singulier.

Le [témoin KKT](../ci/retention_complement.py) utilise un système creux
contraint, équilibré, avec deux corrections par le résidu évalué via
\(D_i\). Son petit bloc est calculé à partir de l'énergie du champ
reconstruit. Résidus, défauts de contrainte, écarts de fonctionnel et
valeurs singulières restent archivés ; ils ne sont pas des bornes dirigées.
Le témoin ne calcule pas encore les bornes de champ et de Schur dérivées
dans le carnet, ni une marge certifiée de type \(\sigma-b\). Sa LU seule
ne garantit pas la détection de tout pôle global après arrondis.

## Résultats conservés

Les entrées sont celles de la [confrontation harmonique](CONFRONTATION_PORTS_EXUDYN_0.10.0.md) :
console aluminium de 1 m, section 40 × 6 mm, 32/128/512 poutres natives,
six coordonnées par nœud, six charges terminales, 257 fréquences de
0 à 40 Hz. Les directions sont choisies par itération spectrale sur le
modèle QR ; leur sélection est approchée. Le certificat vise ensuite le
modèle \(D\) et le \(B\) stocké. La sélection exploratoire extrait douze
directions, dont seules les une ou deux premières servent ici ; son
coût n'est pas inclus dans celui du certificat ci-dessous.

| Poutres | Limite, 1 direction | Limite, 2 directions | Opérations Decimal, 1 / 2 directions |
|---:|---:|---:|---:|
| 32 | 65,4639 Hz | 100,3827 Hz | 67 820 / 90 000 |
| 128 | 65,7001 Hz | 100,9976 Hz | 282 284 / 373 584 |
| 512 | 65,7149 Hz | 101,0366 Hz | 1 140 140 / 1 507 920 |

La limite sans direction retenue était d'environ 28,3 Hz. Ces conversions
en Hz sont indicatives ; la décision de bande utilise les fractions
exactes du minorant et de la pulsation binary64 maximale stockés :
\(\underline\lambda_c>\Omega_{64}^2\).

Le témoin de champ conserve **une** direction et six ports. Les références
indépendantes à 70 et 90 chiffres concordent après conversion binary64.
Le juge mesure l'erreur relative de chaque charge et la norme d'opérateur
sur toutes leurs combinaisons ; les maxima d'opérateur sont :

| Poutres | Masse | Déformation | Port | Seuil 10⁻⁶ |
|---:|---:|---:|---:|:---:|
| 32 | 1,57 × 10⁻¹² | 5,85 × 10⁻¹¹ | 1,71 × 10⁻¹¹ | accepté |
| 128 | 6,07 × 10⁻¹³ | 3,34 × 10⁻¹⁰ | 3,40 × 10⁻¹¹ | accepté |
| 512 | 2,77 × 10⁻¹² | 2,22 × 10⁻⁹ | 1,72 × 10⁻¹⁰ | accepté |

Les valeurs affichées sont arrondies vers le haut. La qualification
physique concerne ces **points échantillonnés** ; la coercivité du
complément est établie sur la bande entière. Il n'y a pas encore
d'enveloppe de réponse certifiée sur toute la bande.

Une seule exécution a été réalisée, sans classement de vitesse. À 512
poutres, les certificats prennent respectivement environ 0,863 et
1,043 s ; le témoin prend 0,089 s de préparation puis 4,695 s de réponses.
Son KKT est de taille 3067 malgré un bloc conservé de taille 7. Ces durées
ne constituent pas un temps total de réduction, puisqu'elles excluent la
sélection des directions, et ne se comparent pas aux médianes des roues.

## Suite mesurée : Krylov contraint

Le [prototype Krylov contraint](KRYLOV_CONTRAINT_PROTOTYPE.md) construit
maintenant cette réduction à partir de l'inverse \(S\), avec les colonnes
des ports et des directions conservées. Sa norme duale retire les
réactions avant les normes ; une réparation après normalisation limite
le défaut des contraintes, ensuite conservé dans l'audit. Les
[preuves d'enveloppe](KRYLOV_CONTRAINT_PREUVES.md) et leur
[variante par direction](KRYLOV_CONTRAINT_ANISOTROPIE.md) composent le résidu
avec la marge de tout le bloc conservé.

La campagne de 60 essais qualifie les 24 essais Krylov à 40 Hz sur les trois
maillages, ainsi que les majorants des 12 essais avec contrôle, avec des
bases totales de 31, 39 et 47 directions. À 512 poutres,
les champs seuls prennent 0,937 s au total contre 0,725 s pour la LU
corrigée ; leur phase de réponse est 7,3 fois plus rapide. Le certificat
reste le poste de préparation dominant. Les réponses et leurs majorants
ne sont toujours pas certifiés machine ; cette campagne reste distincte
du témoin KKT unique décrit ci-dessus.

La correction énergétique actuelle emploie encore un défaut global
\(\bar\eta\), parfois pessimiste sur le complément. Le
[critère par inertie](INERTIE_COMPLEMENT_PREUVES.md) offre une piste pour
tester directement un seuil ; l'[audit de 2026](VERIFICATION_CREUSE_2026.md)
précise les méthodes récentes examinées et un énoncé écarté après
contre-épreuve exacte. Ni l'une ni l'autre de ces extensions n'est une
nouvelle capacité publique ou un gain mesuré face aux solveurs complets.

Les [données et sources mesurées](bancs/complement-spectral-2026/README.md)
et les [sondes et audits adversariaux](bancs/complement-spectral-sondes-2026/README.md)
sont conservés avec empreintes. Les tests de certificats, composition,
inertie, produits flottants et intégrité des résultats entrent dans la CI.
