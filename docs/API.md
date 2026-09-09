# Référence d'API — vinkulum

**Ce fichier est GÉNÉRÉ** par `python -m vinkulum.doc --ecrire`.
Toute édition à la main est écrasée et fait rougir la CI — une
documentation recopiée finit par décrire une autre bibliothèque.

Version : `0.19.0`

## Le noyau — `vinkulum.Noyau`

Le modèle vu de Python : on déclare, on simule, on lit

| méthode | ce qu'elle fait |
|---|---|
| `adapt_stats` | (pas rejoués, dernier pas retenu, max de R_half, max du seuil, max de |Φ| à mi-pas relatif) |
| `aero` | (nom, poussée le long de l'axe d'inflow, couple autour de cet axe) par pale, et (v_i, poussée) par inflow |
| `appariement` | Active l'appariement automatique des sphères déclarées, avec une loi de |
| `appariement_stats` | (paires actives au dernier pas, secondes cumulées d'appariement) |
| `assemble` | Résout les contraintes de position holonomes, puis G·u + ∂Φ/∂t = 0 |
| `audit_jacobien` | AUDIT : jacobien AD de Newton contre différence finie du résidu, à |
| `bilan_stabilite` | Bilan de la linéarisation mécanique locale, sans garantie de stabilité globale |
| `cardan` | CARDAN (croix de Hooke) : les axes `axe_a` (repère de a) et `axe_b` |
| `chronos` | chronos cumulés en secondes : (jacobien, résolution, résidus, fin de |
| `contact` | CONTACT par pénalité : une sphère portée par `corps` (centre `p0` en |
| `contacts` | (nom, enfoncement, force normale) de chaque contact |
| `controle` | CONTRÔLE DU MODÈLE — « est-il bien posé ? », la question qu'un |
| `corps` | Corps rigide ; `j` = tenseur 3×3 en ligne (axes corps, au CdM), `rot` = R en ligne |
| `couple` | Couple à loi entre a et b autour de `axe` (repère de a ; monde si bâti) |
| `couples` | (angle relatif déroulé, vitesse relative, couple) de chaque couple à loi, par nom |
| `d_raideur_poutre` | ∂K/∂(raideur de poutre), sans différence finie : coefficients |
| `d_residu_poutre` | ∂(résidu statique)/∂(raideur de poutre), sans différence finie |
| `d_residu_poutres_transpose` | Produits (∂R/∂p_i)ᵀ poids pour toutes les poutres, dans leur ordre |
| `diag_modal` | DIAGNOSTIC : normes des blocs réduits (Kr, Cr, Mr, T) — pour comparer |
| `distance` | Bielle à deux rotules ; `l` None = longueur initiale |
| `distance_maillage` | Distance d'un point MONDE au maillage `i`, et le point le plus proche |
| `effort` | Effort constant (F, M) au CdM du corps, repère monde |
| `effort_temporel` | Effort temporel au point `point` local au corps. F et M sont mondiaux ; |
| `energie` | Énergie mécanique totale : cinétique (translation + rotation, tenseur |
| `engrenage` | Engrenage θ_a = rapport·θ_b, axes `axe_a`/`axe_b` dans le repère du |
| `enregistre_schema` | Arme (ou désarme) l'HISTORIQUE DU SCHÉMA : (u̇, a, λ) à l'état initial |
| `etat` | Vue f64 : (t, [r], [R en ligne], [v], [ω], [v_i des inflows]) |
| `etat_precis` | État cinématique précis : les six champs de `etat()`, suivis de [r_bas] |
| `etats_inflow` | (v_i, v1s, v1c) en m/s — les trois états d'inflow |
| `etats_lb` | États de Leishman–Beddoes par station de la pale `i` : |
| `etats_wagner` | États de Wagner (z₁, z₂) par station, pour la pale `i` |
| `facteurs_materiels_poutres` | Facteurs matériels des seules poutres, à l'état courant, sans K ni Z |
| `fuite_transport` | DIAGNOSTIC de l'affirmation « le transport sort de l'espace admissible » : |
| `inflow` | Inflow uniforme : axe (monde), aire du disque, retard τ, masse volumique |
| `instationnaire` | Bascule une pale en AÉRODYNAMIQUE INSTATIONNAIRE (Wagner/Jones) |
| `k_c_m_g_creux` | (K, C, M, G) en stockage CSC, autour de l'état courant |
| `k_c_m_z` | (K, C, M, Z, G) autour de l'état courant : raideur tangente, |
| `k_m_z` | (K, M, base admissible Z) autour de l'état courant, en ligne — de quoi |
| `liaison` | LIAISON GÉNÉRIQUE — six degrés à la carte, entre `a` et `b` |
| `liaison_reperes` | Liaison définie par DEUX repères locaux indépendants. Contrairement à |
| `maillage` | MAILLAGE TRIANGULAIRE porté par un corps — la géométrie quelconque |
| `modes` | Modes propres du système contraint autour de l'état courant : |
| `modes_complexes` | Modes oscillants COMPLEXES : [(fréquence Hz, amortissement réduit ζ, σ)] |
| `modes_creux` | Modes locaux creux : valeurs propres signées proches du décalage, |
| `moment` | Moment cinétique total, au repère monde et autour de l'origine |
| `moyenne_aero` | Arme la moyenne du torseur aérodynamique à partir de `t0` (remise à zéro) |
| `pale` | Pale (théorie des tranches, polaire c81 mono-Mach). `p0` pied d'envergure |
| `pas_inflow` | UN PAS D'INFLOW SEUL, charges IMPOSÉES — le point d'entrée qui rend le |
| `phi` | Φ à l'état courant : la VIOLATION de chaque ligne de contrainte, dans |
| `phi_dot` | Φ̇ = G·u + ∂Φ/∂t à état fixé ; date courante par défaut, sans mutation |
| `pitt_peters` | Bascule un inflow en PITT–PETERS à trois états. `centre` = centre du |
| `pose` | (position, rotation en ligne) d'UN corps — le raccourci quand `etat()`, |
| `pose_effort` | Change un effort déjà déclaré. C'est ce qui rend la CONTINUATION |
| `pose_etat` | Restaure un état rendu par `etat`. Les rotations sont RÉ-ORTHONORMALISÉES |
| `pose_inflow` | Force la vitesse induite de l'inflow `i` (m/s). Sert à partir d'un état |
| `pose_inflow_carte` | Impose un PROFIL RADIAL d'inflow v(r̄) sur le disque `i` (m/s, > 0 vers |
| `pose_inflow_harmoniques` | Impose l'inflow d'un disque avec ses harmoniques 1/rev — ce qu'un |
| `pose_inflow_profil` |  |
| `pose_lam` | Pose les multiplicateurs λ de l'état courant (après `pose_etat`) : la |
| `poutre` | Poutre géométriquement exacte entre les corps a et b (nœuds, à leur pose |
| `poutres` | (γ, κ) de chaque poutre — déformation et courbure matérielles courantes |
| `reactions` | (nom, multiplicateurs λ) par élément — les EFFORTS DE LIAISON, dans le |
| `residu_statique` | Résidu STATIQUE à l'état courant : −(forces + gravité), sur les 6n |
| `schema` | L'historique armé par `enregistre_schema` : liste de (u̇, a, λ), une |
| `simule` | Intègre jusqu'à `t_end` au pas `h`. `adaptatif` = tolérance RELATIVE du |
| `simule_em` | Intégrateur ÉNERGIE-MOMENT (Simo–Wong, Cayley au point milieu) — E et |
| `simule_multirythme` | SOUS-CYCLAGE MULTI-RYTHME : les corps `rapides` à h/k, les autres à h, |
| `spectre` | Toutes les racines mécaniques locales : [(réel, imaginaire)] en s⁻¹ |
| `sphere` | Déclare une sphère candidate à l'APPARIEMENT automatique |
| `statique` | ANALYSE STATIQUE — l'équilibre, sans passer par le temps |
| `statique_info` | Rapport historique du dernier calcul statique ; None si aucun rapport |
| `stats` | (itérations de Newton cumulées, replis SVD cumulés, jacobiens calculés, contraintes redondantes neutralisées) |
| `superelement` | SUPERÉLÉMENT : la raideur d'un maillage EF quelconque, attachée à des |
| `t` | Le temps courant du modèle, en secondes |
| `torseur_moyen` | (F, M) aérodynamiques MOYENS depuis `t0`, réduits à l'origine du monde, |
| `vent` | Vitesse de l'air ambiant (repère monde). Le rotor de banc étant encastré, |
| `vis` | VIS–ÉCROU : la translation relative le long de `axe` est liée à la |

## Les modules d'analyse

### `vinkulum.certification`

Certificats a posteriori : unicité et erreur en avant des systèmes linéaires

| fonction | ce qu'elle fait |
|---|---|
| `certifier_assemblage(noyau, *, jauges, rayons)` | Certifie une pose admissible unique dans une boîte explicitement donnée |
| `certifier_initialisation(noyau, *, t=None, erreur_max=None, dimension_max=64, redondances=False)` | Borne l'erreur du système linéaire d'accélération réellement résolu par le noyau |
| `certifier_quotient_structurel(m, c, t, base, force, d, x=None, *, delta_m=None, delta_c=None, delta_t=None, delta_force=None, delta_d=None, poids=None)` | Certifie un KKT perturbé et toutes les contraintes G'=T'C', c_complet=T'd' |
| `certifier_systeme(a, b, x=None, *, inverse_approche=None, erreur_max=None, dimension_max=64)` | Certifie l'unicité de Ax=b et borne chaque composante de l'erreur de x |
| `certifier_systeme_creux(a, b, x=None, *, facteurs=None, permutations=None, delta_a=None, delta_b=None, poids=None)` | Certifie des systèmes creux perturbés sans inverse dense, jusqu'à 4096 inconnues |
| `verifier_certificat(certificat)` | Vérifie un certificat par un second calcul rationnel, sans solveur flottant |

### `vinkulum.trim`

vinkulum.trim — Newton sur les commandes, à matrice d'influence GARDÉE

| fonction | ce qu'elle fait |
|---|---|
| `demo()` | Autotest : cas algébriques à solution connue, et le PRIX de la stratégie |
| `trim(residu, x0, *, pas=0.001, tol=1e-06, iters=30, seuil_refaire=0.5, echelle=None, journal=None, broyden=True)` | Résout residu(x) = 0 par Newton à matrice d'influence gardée |

### `vinkulum.rotor`

vinkulum.rotor — rotor d'essai en articulé équivalent, et son TRIM

| fonction | ce qu'elle fait |
|---|---|
| `demo(rapide=False)` |  |
| `polaire_lineaire(a_lift=5.73, cl0=0.1, cd0=0.0079, k_cd=0.4, alpha_dec=14.0)` | Polaire c81 mono-Mach : linéaire jusqu'au décrochage, puis plateau |

### `vinkulum.floquet`

vinkulum.floquet — stabilité d'un système PÉRIODIQUE, par la monodromie

| fonction | ce qu'elle fait |
|---|---|
| `demo(rapide=False)` |  |
| `dmonodromie(construit, periode, h, base, param, p0, rel=0.0001, **kw)` | ∂Φ/∂p par différence centrée sur DEUX monodromies |
| `ecart(etat, ref)` | Perturbation de `etat` par rapport à `ref`, en 12·n composantes |
| `exposants(phi, periode, f_propre=None)` | Table triée par |λ| décroissant |
| `monodromie(construit, periode, h, base, eps=1e-06, tol_phi=1e-07, tol_per=1e-06, journal=False)` | Matrice de transition sur une période, réduite à `base` |
| `sensibilite(phi, dphi, periode)` | ∂σ/∂p de CHAQUE mode, sans re-diagonaliser et SANS APPARIEMENT |

### `vinkulum.reduction`

vinkulum.reduction — Craig–Bampton : un corps flexible réduit à ses interfaces

| fonction | ce qu'elle fait |
|---|---|
| `craig_bampton(k, m, idx_r, n_modes)` | Réduit (K, M) en gardant les DDL `idx_r` et `n_modes` modes internes |
| `demo(rapide=False)` |  |
| `frequences(k, m, combien=6)` | Fréquences propres (Hz) d'un couple (K, M) symétrique, M définie positive |

### `vinkulum.reduction_ports`

Réduction matérielle par ports : énergie d'entrée et contrôle du champ

| fonction | ce qu'elle fait |
|---|---|
| `reduire_poutres(noyau, libres, interface, metrique, omega_max, **options)` | Construit une réduction de la partie matérielle des poutres du noyau |

### `vinkulum.reduction_contrainte`

Réduction fréquentielle avec directions retenues et complément certifié

### `vinkulum.sensibilite`

vinkulum.sensibilite — le gradient par l'ADJOINT : une résolution pour tous

| fonction | ce qu'elle fait |
|---|---|
| `adjoint(construit, objectif, params, eps=None, t_end=5.0, h=0.002, exact=None)` | Gradient de `objectif` par rapport à `params`, par l'adjoint |
| `demo(rapide=False)` |  |
| `equilibre(construit, t_end=5.0, h=0.002, tol_v=1e-06)` | Amène le modèle à son équilibre statique par relaxation amortie |
| `gradient_df(construit, objectif_val, params, eps=None, t_end=5.0, h=0.002)` | Gradient par différences finies GLOBALES : 2N équilibres complets |
| `gradient_frequence(N, mode, exact)` | dλ/dp et df/dp EXACTS pour un mode propre, sans re-résolution |

### `vinkulum.convergence`

vinkulum.convergence — « mon pas de temps est-il assez fin ? »

| fonction | ce qu'elle fait |
|---|---|
| `demo()` | Autotest, en DEUX cas — et leur différence est l'enseignement |
| `etude(monte, t_end, h, facteurs=(1, 2, 4), grandeur=None, rho=0.6)` | Rejoue la MÊME simulation à h, h/2, h/4… et mesure ce qui bouge |
| `verdict(monte, t_end, h, tol, **kw)` | Accord du calcul au premier pas avec le plus fin à `tol` près |

### `vinkulum.adjoint_temps`

vinkulum.adjoint_temps — le gradient d'une TRAJECTOIRE, par l'adjoint discret

| fonction | ce qu'elle fait |
|---|---|
| `coeffs(rho)` | Les quatre paramètres de l'α-généralisé, comme le noyau les pose |
| `demo(rapide=False)` |  |
| `demo_nonlineaire()` | L'adjoint NON LINÉAIRE contre une différence finie, sur un DUFFING |
| `demo_pont(rapide=False)` | Le pont contre une différence finie CENTRÉE sur le noyau lui-même |

### `vinkulum.rotation`

vinkulum.rotation — le raidissement centrifuge d'une poutre qui tourne

| fonction | ce qu'elle fait |
|---|---|
| `demo(rapide=False)` |  |
| `frequence(omega, t_end=1.5, h=0.0002, **kw)` | Première fréquence de flexion à `omega`, une fois le régime ÉTABLI |
| `poutre_tournante(omega, n=6, L=1.0, b=0.02, hh=0.004, rho=2700.0, e_mod=70000000000.0, ei_fac=1.0, amorti=0.05)` | Poutre encastrée sur un moyeu en rotation imposée à `omega` |

### `vinkulum.contact`

vinkulum.contact — le contact par pénalité régularisée, et ce qu'il vaut

| fonction | ce qu'elle fait |
|---|---|
| `appariement(rapide=False)` | L'APPARIEMENT automatique : exactitude, puis complexité |
| `barriere_ipc(k=1000.0, d_hat=0.001, c=5.0, m=0.5, R=0.02, t_end=1.5, h=2e-05)` | Bille posée, loi de BARRIÈRE IPC : elle ne touche jamais le sol |
| `bille(k=1000000.0, expo=1.5, c=0.0, mu=0.0, m=0.5, R=0.02, z0=None, vz=0.0, nonlisse=False, e=0.0)` |  |
| `boites()` | LA BOÎTE (OBB) — troisième primitive, et la dernière du plan v0.5 |
| `capsules()` | LA CAPSULE — un segment dilaté, et la sphère en est le cas dégénéré |
| `capsules_primitives()` | LA CAPSULE contre BOÎTE, CYLINDRE et MAILLAGE — ce que le contact ne |
| `choc(m1, m2, v1, c=0.0, R=0.02, k=1000000.0, h=2e-06, t_max=0.2)` | Choc central de deux billes libres, sans pesanteur |
| `cylindres()` | LE CYLINDRE FINI — la primitive qui manquait, et ses QUATRE régimes |
| `demo(rapide=False)` |  |
| `maillages()` | GÉOMÉTRIE QUELCONQUE — maillage triangulaire et BVH |
| `multiples(h=0.0001)` | CONTACTS MULTIPLES QUI BASCULENT — ce qu'une bille seule ne teste pas |
| `nonlisse()` | CONTACT NON LISSE — multiplicateur sous complémentarité, loi d'impact en vitesse |
| `pente(th_deg, mu, t_end=0.4, h=2e-05, m=0.5, k=1000000.0, R=0.02, nonlisse=False)` | Bille sur un plan incliné PORTÉ PAR UN CORPS : accélération le long de |
| `rebond(h_libre, h_contact=None, marge=0.0, t_end=1.2, k=10000000000.0, c=0.6, R=0.02, m=0.5, z0=0.3)` | Bille lâchée qui rebondit : le vol libre est long, le choc est bref |
| `restitution(c, v0=0.5, k=1000000.0, m=0.5, R=0.02, h=2e-06)` | Rebond : vitesse après / vitesse avant, sans gravité (choc pur) |
| `tas(n, R=0.01, k=100000.0, mu=0.3, t_end=0.6, h=0.0002, m=0.01, graine=11)` | N billes lâchées en pluie sur un sol, contacts découverts et IPC |
| `tir(v0=50.0, h=0.001, ccd=False, k=100000000.0, R=0.005, m=0.01, d_hat=0.001)` | Bille lancée à grande vitesse contre une autre : traverse-t-elle ? |
| `tir_paroi(forme, ccd, v0=50.0, h=0.001, k=100000000.0, R=0.005, m=0.01, d_hat=0.001)` | Bille rapide contre une PAROI mince fixée au bâti — boîte, cylindre ou |

### `vinkulum.coque`

PLAQUE EN FLEXION — un maillage d'ÉLÉMENTS FINIS que le noyau ne sait pas

| fonction | ce qu'elle fait |
|---|---|
| `demo()` | La plaque doit retrouver la POUTRE quand elle est étroite et ν = 0 |
| `plaque(lx, ly, nx, ny, e, nu, h, rho)` | Maille une plaque rectangulaire en `nx × ny` éléments ACM |
| `pont(nx=16, ny=2, garde=9)` | LE PONT DE BOUT EN BOUT : une plaque EF devient des CORPS du multicorps |

### `vinkulum.andrews`

LE MÉCANISME D'ANDREWS — un modèle que l'auteur du solveur n'a pas conçu

| fonction | ce qu'elle fait |
|---|---|
| `demo(t_end=0.05, h=2e-06, mov=None)` | Andrews : vinkulum contre MBDyn, sur une définition qu'aucun n'a choisie |
| `monte(couple=0.033)` | Le mécanisme dans vinkulum. Rend (Noyau, index des corps) |
| `reference_mbdyn(chemin=None)` | Relit la référence depuis un `.mov` de MBDyn, si on en a un |

### `vinkulum.verification`

vinkulum.verification — la vérification APPROFONDIE du jalon 1 (hors gate, minutes)

| fonction | ce qu'elle fait |
|---|---|
| `accords()` | DEUX BRIQUES POUR LA MÊME PHYSIQUE DOIVENT S'ACCORDER |
| `amortissement()` | ρ∞ = 1 (trapèze) est INSTABLE en index 3 — mesuré : Newton diverge à |
| `arbre_tournant()` | rotatingshaft — « Stability of a rotating shaft », le second banc publié |
| `assemblage()` | ASSEMBLER UN MÉCANISME depuis une pose APPROCHÉE — Φ(q) = 0, puis Φ̇ = 0 |
| `audit_jacobien()` | LE JACOBIEN DE NEWTON CONTRE UNE DIFFÉRENCE FINIE DU RÉSIDU, élément par élément |
| `bennett()` | LE MÉCANISME DE BENNETT (1903) — un modèle que l'auteur du code n'a PAS |
| `convergence()` |  |
| `coulomb_pente()` | FROTTEMENT DE COULOMB — les deux régimes, et le seuil qui les sépare |
| `determinisme()` | Même calcul, RAYON_NUM_THREADS=1 puis tous les fils : identique au bit |
| `echelle()` | JUSQU'OÙ ÇA MONTE — mesuré, pas supposé |
| `gyroscopique_tangent()` | LE GYROSCOPIQUE N'EST COMPTÉ QU'UNE FOIS dans les tangentes exposées |
| `invariance_unites()` | LE MÊME MÉCANISME, COTÉ DANS CINQ SYSTÈMES D'UNITÉS |
| `invariances()` | CE QU'UN UTILISATEUR SUPPOSE SANS LE VÉRIFIER |
| `kapitza()` | LE PENDULE DE KAPITZA (1951) — second problème posé ailleurs, d'une |
| `lateral_buckling()` | Flambement latéral d'une poutre mince — le banc publié de Multibody |
| `lecteurs_casses()` | UN UTILISATEUR TIERS APPORTE AUSSI DES FICHIERS CASSÉS |
| `liaisons_unilaterales()` | BUTÉE et CÂBLE — les deux liaisons unilatérales que tout mécanisme a |
| `main()` |  |
| `modeles_mal_poses()` | MODÈLES PLAUSIBLES MAL POSÉS — pas absurdes, JUSTE FAUX |
| `non_holonome()` | ROULEMENT SANS GLISSEMENT — la brique qui manquait pour qu'un véhicule |
| `quatre_barres_flexible()` | fourbar — le quatre-barres de Bauchau, troisième banc publié de |
| `refus_gardes()` | CE QU'ON A REFUSÉ — la RAISON du refus, asserée et non recopiée |
| `robustesse()` | CE QUE LE NOYAU FAIT DEVANT UN MODÈLE ABSURDE — un solveur mature ne |
| `six_barres()` | 6barmech — le mécanisme à six barres du jeu de tests de MBDyn |
| `sparsite()` | LA STRUCTURE DU PROBLÈME, contre le remplissage réel |
| `srscm()` | srskm — le *spatial slider-crank* du jeu de tests de MBDyn |
| `statique()` | ANALYSE STATIQUE et POUTRE DE PRINCETON — le benchmark expérimental |
| `superelement()` | SUPERÉLÉMENT — la raideur d'un maillage EF quelconque, dans le multicorps |
| `toupie()` |  |
| `treillis_de_boucles()` | multibarmech — le dernier cas du jeu de tests MBDyn, et le plus GROS : |
| `validation_scores()` | Une hypothèse qui passe ne doit pas masquer sa voisine qui échoue |
| `vis_ecrou()` | VIS–ÉCROU (glissière hélicoïdale) — translation liée à la rotation |

### `vinkulum.bancs`

vinkulum.bancs — LA SUITE DE BANCS : courbes travail–précision, le juge du noyau

| fonction | ce qu'elle fait |
|---|---|
| `adaptatif()` | CE QUE LE PAS ADAPTATIF COÛTE — un résultat NÉGATIF, mesuré et gardé |
| `amortis()` | Modes COMPLEXES contre trois solutions exactes : oscillateur amorti |
| `bielle(h)` |  |
| `bifurcation()` | LE MASQUE DE REDONDANCE N'EST PAS FIGÉ — un mécanisme parti d'un POINT DE |
| `chaine(n_corps, h=0.001, t_end=0.2, L=0.05, m=0.02)` | Chaîne de `n_corps` maillons en rotules, lâchée à l'horizontale sous |
| `console(n_elem, h=None)` | Poutre géométriquement exacte encastrée (alu 1 m, 20 × 5 mm), chargée en |
| `dissipation()` | CE QUE ρ∞ COÛTE ET CE QU'IL ACHÈTE — le compromis, mesuré |
| `double_pendule(h)` |  |
| `echelle(tailles=(10, 30, 100))` |  |
| `energie_moment()` | SCHÉMA ÉNERGIE-MOMENT DANS LA FORMULATION — l'exclusive-ou de la projection, résolu |
| `engrenage(h)` |  |
| `ggl_boucles()` | GGL ÉCHOUE SUR TROIS MODÈLES — reproduire l'échec hors de FRELON |
| `inflow_carte()` | CARTE D'INFLOW w(r̄, ψ) imposée (7 sept.) — ce qui rend Peters–He et un |
| `instationnaire()` | AÉRODYNAMIQUE INSTATIONNAIRE DE SECTION — Wagner, jugée par Theodorsen |
| `leishman_beddoes()` | DÉCROCHAGE DYNAMIQUE DE LEISHMAN–BEDDOES — ce qu'Øye ne fait pas |
| `main(rapide=False, sortie=None)` |  |
| `modes_console(n_elem=20)` | Les TROIS premiers modes de flexion d'une console, par analyse modale du |
| `moment_projete()` | NOETHER DISCRET PAR PROJECTION — la réparation, et ses quatre réfutations |
| `multirythme()` | SOUS-CYCLAGE MULTI-RYTHME — la partition lente ne paie plus le pas de la rapide |
| `noether()` | LE MOMENT CINÉTIQUE DÉRIVE-T-IL ? — et de combien VRAIMENT |
| `oye()` | DÉCROCHAGE DYNAMIQUE D'ØYE — échelon d'incidence sur une section fixe |
| `pendule(h)` |  |
| `pitt_peters()` | INFLOW DYNAMIQUE DE PITT–PETERS (1981) — trois états, trois références |
| `poutre(h)` | Banc à pas fixe : l'erreur décroît avec le NOMBRE D'ÉLÉMENTS, pas avec h — |
| `quatre_barres(h)` |  |
| `raideur_gravite()` | La raideur tangente d'une toupie AU REPOS doit valoir exactement le |
| `seuils()` | LES SEUILS DE NEWTON SONT-ILS LA PHYSIQUE OU LE BANC ? — mesuré |
| `theodorsen_complet()` | THEODORSEN COMPLET — masse ajoutée + incidence au 3/4 de corde (6 sept.) |
| `toupie(h)` |  |
| `violation_vitesse()` | Φ̇ EST VIOLÉ À O(h²) EN INDEX 3, ET GGL LE FERME SANS PERDRE L'ORDRE — mesuré |
