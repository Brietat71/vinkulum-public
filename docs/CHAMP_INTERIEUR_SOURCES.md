# Champs, observables et enrichissement : sources primaires

Lecture du 7 septembre 2026. Cette recherche ciblée complète les carnets
sur les ports ; elle ne constitue pas une nouvelle comparaison de solveurs.
Les énoncés algébriques adaptés au prototype sont démontrés séparément dans
[le carnet de preuves](CHAMP_INTERIEUR_PREUVES.md).

| Source lue | Énoncé utile et conditions | Limite pour Vinkulum |
|---|---|---|
| Rannacher, ICM 2002, prépublication 2003, *Adaptive Finite Element Methods for Partial Differential Equations*, propositions 1–2, équations (2.1), (3.3). [Texte primaire](https://arxiv.org/pdf/math/0305006) | Représentation de l'erreur d'observable par résidus primal et dual ; reste cubique sous régularité et stationnarité. Dans le cas linéaire, le reste disparaît. | Un dual approché ne donne pas seul une borne : conserver son résidu et la constante de stabilité. L'antériorité du contrôle par observable est explicite. |
| Smetana–Patera, 2016, *Optimal Local Approximation Spaces for Component-Based Static Condensation Procedures*, proposition 3.5 et théorème 4.4. [Article primaire](https://ris.utwente.nl/ws/portalfiles/portal/168414349/15m1009603.pdf) | L'espace spectral de l'opérateur de transfert compact réalise la largeur de Kolmogorov, d_n = √λ_(n+1), dans les métriques prescrites. Le spectral greedy contrôle l'ensemble fini d'apprentissage. | Traiter les mouvements rigides et les charges ; aucune extension automatique à toute fréquence, aux contraintes ponctuelles ou à une autre métrique. |
| Buhr–Smetana, 2018, *Randomized Local Model Order Reduction*, §3.3, propositions 3.7–3.8. [Article primaire](https://arxiv.org/pdf/1706.09179) | Un estimateur gaussien avec constante explicite majore la norme de l'erreur de transfert avec probabilité prescrite ; la borne d'union doit compter les étapes adaptatives. | Métriques SPD et constantes spectrales adéquates requises. Une probabilité d'échec contrôlée ne certifie pas les arrondis. Aucun échantillonnage aléatoire n'est utilisé dans le prototype actuel. |
| Buhr, 2018, *Exponential Convergence of Online Enrichment in Localized Reduced Basis Methods*, théorèmes 1 et 4. [Texte primaire](https://arxiv.org/pdf/1710.02104) | Contraction énergétique de facteur √(1−1/(N_D c_pu²)) par correcteur du résidu local maximal, pour une forme symétrique coercive et des sous-domaines recouvrants avec partition de l'unité stable. L'enrichissement couplé est optimal dans les espaces locaux considérés. | Le contraste et le recouvrement entrent dans c_pu. La preuve ne couvre pas directement notre dynamique globalement indéfinie ni nos rondes d'enrichissement Krylov. |
| Smetana–Taddei, 2023, *Localized Model Reduction for Nonlinear Elliptic Partial Differential Equations*, SISC 45(3), A1300–A1331. [Publication](https://doi.org/10.1137/22M148402X), [préprint lu](https://www.math.u-bordeaux.fr/~ttaddei/data/KSTT_arxiv.pdf) | Lemme 5.1 et proposition 5.1 du préprint : résidus locaux vers résidu global, puis théorie de Brezzi–Rappaz–Raviart. Avec stabilité β, Lipschitz L et résidu R, τ=2LR/β²<1 donne une boule contrôlée de rayon β/L·(1−√(1−τ)). | Constantes quantitatives et proximité nécessaires ; une inégalité inverse dépend du maillage. La contraction d'enrichissement prouvée concerne le cas linéaire coercif. Ce théorème ne permet pas de déclarer résolue l'extension multicorps non linéaire. |

Les recherches ont ciblé les représentations d'erreur d'observable, les
opérateurs de transfert et les théorèmes d'enrichissement. Les textes ont
été rapprochés des hypothèses du code ; une seconde lecture indépendante
a vérifié les passages décisifs de Rannacher, Buhr et Smetana–Taddei.
La recherche s'arrête ici parce que les cinq familles de claims ont une
source primaire ou une limite explicite. Aucun code concurrent n'a été lu.

**Conclusion de provenance :** l'erreur primal-dual, les espaces de transfert
optimaux et l'enrichissement résiduel ont des antériorités établies. Le travail
propre à Vinkulum est leur composition contrôlée avec le facteur énergétique
d'entrée, les masses couplées, les applications de ports et les résidus
flottants effectivement conservés. Aucune nouveauté mondiale n'est établie.
