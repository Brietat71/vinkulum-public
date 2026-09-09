# Revues de notes externes

Notes d'architecture transmises au projet, triées ligne par ligne contre ce que
vinkulum est déjà et ce que le métier fait. Sorties du README le 4 septembre :
ce sont des réponses datées à des textes précis, pas une description du noyau.

## Réponse à la note « Aether-Kernel » (3 septembre) — ce qu'on prend, ce qu'on refuse

Paul a transmis une note d'architecture pour un « super-noyau » CAO + FEA + MBD
en Rust. Tri, ligne par ligne, contre ce que vinkulum est déjà et ce que le
métier fait :

**On prend (et c'est la voie déjà engagée)** : Rust sans code hérité ;
portabilité par LLVM (x86, ARM, RISC-V viennent gratis) ; **couplage
monolithique** — c'est exactement notre système DAE unique, où les degrés de
liberté FEA s'ajouteront aux corps et aux contraintes comme des inconnues de
plus ; rayon (déjà là) ; **AD pour les jacobiens** ; **SoA** pour les boucles
d'ÉLÉMENTS (pas pour les douze corps d'une tête rotor — là c'est sans objet) ;
allocateur mimalloc quand l'assemblage pèsera ; **faer** plutôt que
nalgebra+sprs+PETSc : dense ET creux, SIMD portable via `pulp` sur stable,
rayon dedans — un seul socle natif.

**On refuse, avec le motif** :
- *« Fin du B-Rep, place au SDF »* — non. Un noyau qui vise Parasolid/CGM doit
  porter la géométrie EXACTE (NURBS, tolérances, STEP) : c'est ce que
  l'usinage, la cotation et l'échange exigent. Le SDF est une structure
  DÉRIVÉE, excellente pour le contact et le maillage, jamais la représentation
  maîtresse. `parry` est un détecteur de collisions (celui de Rapier), pas un
  modeleur.
- *PETSc par FFI pour « des milliards de DDL »* — contradiction avec « zéro
  code hérité », et hors échelle : un solveur creux direct natif (faer) tient
  10⁶–10⁷ DDL sur un nœud ; la distribution MPI est un sujet d'après, s'il
  vient.
- *`std::simd`, `crossbeam` lock-free, atomiques Acquire/Release pour
  l'assemblage, `hwloc`* — prématuré et en partie faux : l'assemblage parallèle
  se fait par tampons par fil puis réduction (COO → CSR), pas par atomiques ;
  `std::simd` est nightly, `pulp` est stable ; hwloc/NUMA compte à 128 cœurs,
  pas avant.
- *Enzyme comme « cœur névralgique »* — l'AD en Rust par Enzyme
  (`std::autodiff`) est expérimental/nightly à ma connaissance : on démarre
  avec des nombres duaux sur stable (`num-dual`, mode avant — il nous faut des
  jacobiens de petites fonctions par élément, c'est son cas d'usage), et on
  bascule sur Enzyme quand il sera stable.

**Ce que la note ne dit pas, et qui est tout le sujet** : aucune formulation
(quel intégrateur, quelles contraintes, quels éléments, quel contact), aucune
VÉRIFICATION (pas un cas analytique, pas une régression contre un solveur
indépendant), rien sur la licence et le clean-room, rien sur
l'aérodynamique — notre cible. C'est une liste de crates, pas une
architecture ; l'architecture est ce que ce README décrit au-dessus, brique
par brique, chacune jugée avant d'être étendue.

## Réponse à la version étayée de la note Aether (3 septembre, soir) — lecture critique

Texte relu point par point, sources vérifiées. Ce qui tient, ce qui ne tient
pas, ce qu'on en garde.

**D'où vient le texte.** Le vocabulaire « Région / Rep Router / certificats /
CutFEM + ghost penalty + Nitsche » est celui du site **FrankenSim**
(frankensim.org — « Certified Simulation & Design Kernel for Rust », MIT,
v0.0.1, 100+ crates, 160 k lignes, pas sur crates.io). La note y ajoute SBM,
GGL, port-hamiltonien, FETI-DP, Enzyme, Arrow, wgpu/egui, que FrankenSim ne
mentionne pas. Les références [1]–[4] sont citées pour tout indistinctement
(une FIGURE ResearchGate de la finite cell method est censée étayer GGL et
FETI-DP) ; le texte finit par un avertissement médical, et parle de parois
artérielles et de scans médicaux : sa matière première est de la mécanique
biomédicale, pas de la voilure tournante. C'est un texte génératif à
citations décoratives — à lire pour les idées, pas comme une preuve.

**Ce qui est vrai et bien décrit.** CutFEM + ghost penalty (Burman 2010) et
Nitsche : exact. SBM (Main & Scovazzi 2018 ; haut ordre par Taylor, 2022 ;
Gap-SBM 2025) : exact, y compris le fait que le vecteur distance vient du
gradient d'un vrai champ de distance. GGL (Gear–Gupta–Leimkuhler 1985)
comme réduction d'index 3 → 2 : exact. FETI-DP/BDDC polylogarithmiques :
exact. Enzyme et le problème d'analyse de type en Rust (travaux GSoC,
métadonnées depuis la MIR) : exact et à jour ; `ad-trait` (Yale 2025,
forward + reverse, tangentes SIMD) : existe. Arrow/PyO3 zéro-copie : exact.

**Ce qui est faux, ou passé sous silence.**
- *« L'index 3 dérive, il faut GGL »* — faux pour un intégrateur qui résout
  Φ(q) = 0 à chaque pas par Newton : vinkulum le fait, mesuré |Φ| 7e-15 après
  4 000 pas et 7e-19 en boucle fermée. La dérive frappe les réductions
  d'index 1/2 SANS projection. GGL est une alternative valable (Adams l'a en
  option), pas une nécessité ; et l'α-généralisé en index 3 a sa preuve de
  convergence (Arnold–Brüls 2007).
- *« Symplectique, plus d'amortissement numérique »* — l'amortissement
  numérique contrôlé (ρ∞) est une FONCTIONNALITÉ pour le multicorps flexible
  raide : sans lui, les hautes fréquences parasites du maillage ne meurent
  jamais. Un intégrateur variationnel conserve l'énergie, y compris celle du
  bruit.
- *SDF comme représentation maîtresse* — le min/max régularisé d'une CSG
  n'est PAS un champ de distance (seulement une borne) : le gradient n'est
  pas unitaire près des arêtes, or SBM a besoin du vrai vecteur distance —
  la note le dit elle-même et ne voit pas la contradiction ; il faut
  redistancer (eikonal). Et une CSG régularisée n'est plus la géométrie
  exacte : un cylindre cesse d'être un cylindre. Le B-Rep reste maître.
- *Éléments immergés pour une machine comme FRELON* — la coque fait 0,6 à
  1 mm, les pales sont des poutres : une grille de fond doit être bien plus
  fine que l'épaisseur, soit des millions de cellules pour ce qu'une poutre
  géométriquement exacte ou une coque conforme fait en centaines de DDL. Les
  méthodes immergées brillent sur l'optimisation topologique et le
  fluide-structure volumique, pas sur les structures minces.
- *FETI-DP, XFEM, VEM, « localisation d'Anderson »* — remplissage. Pertinent
  à 10⁷ DDL distribués ; un creux direct natif (faer) fait tout ce dont un
  hélicoptère a besoin sur un nœud.
- *Highway (C++) « ou équivalents Rust », `std::simd`* — contradiction avec
  « zéro C/C++ », et `std::simd` est nightly.
- *Enzyme « pilier »* — `std::autodiff` avance vers nightly ; le projet Rust
  a écrit qu'il **ne demandera pas de stabilisation pendant la période 2026**
  et envisage même une réimplémentation plus petite. Bâtir dessus aujourd'hui,
  c'est bâtir sur nightly.

**Ce qu'on en garde pour vinkulum.**
1. AD en mode avant par surcharge, sur stable (`ad-trait` ou `num-dual`),
   pour les jacobiens de liaisons et d'éléments ; Enzyme quand il sera stable.
2. Binding numpy zéro-copie (`rust-numpy`) — aujourd'hui `simule` rend des
   listes, c'est une copie.
3. Le formalisme port-hamiltonien comme DISCIPLINE de couplage multiphysique
   — il est vraiment utile pour ce que FRELON a et que les noyaux MBD
   ignorent : la chaîne batterie → ESC → moteur → étage → rotor, à énergie
   conservée. Pas comme solveur.
4. GGL en OPTION, le jour où un cas le justifie (mesure, pas doctrine).
5. FrankenSim : à lire (MIT), pas à dépendre (v0.0.1, périmètre démesuré).
