# Vinkulum 0.9.0 — adjoints par produits implicites

La 0.9.0 réduit le coût du gradient d'une trajectoire mécanique : le pont
temporel résout directement un système transposé et applique les dérivées,
sans construire les grandes matrices de transition d'état et de paramètres.
Le noyau assemble et exporte `K,C,M,G` en stockage CSC sans base dense `Z`.

Les deux nouvelles méthodes `Noyau.k_c_m_g_creux` et
`Noyau.d_residu_poutres_transpose`, ainsi que les rappels optionnels
`d_residu_transpose` et `d_initial_transpose` de `PontNoyau`, étendent l'API.
Cette extension compatible justifie une **version mineure** selon les
[règles du projet](VERSIONNEMENT.md).

Les appels existants au pont restent valides. La sélection automatique
utilise SciPy à partir de 180 degrés physiques s'il est disponible ;
`creux=False` conserve NumPy et `creux=True` exige SciPy.
Le produit de sensibilité d'une poutre inclut les deux formulations et
les six directions de rigidité déjà disponibles.

Le Jacobien adjoint est équilibré et son erreur arrière contrôlée, avec
raffinement éventuel. L'état initial garde sa différence finie historique
à moins que l'utilisateur fournisse son produit de dérivée exact. Cette
possibilité évite deux constructions initiales du modèle par paramètre,
mais exige de connaître toutes ses dépendances.

Le pont refuse désormais aussi les candidats de contact automatique,
les éléments aérodynamiques et les liaisons non holonomes : leurs états
ou leurs équations ne sont pas dérivés par le recul mécanique actuel.
GGL, correction σ, adaptatif, événements, changements de branche et
historiques internes non enregistrés restent hors de son domaine.

La campagne finale conserve **156 processus**, échauffements compris.
À 120 poutres, le gradient d'un paramètre partagé passe de **10,846 s à
0,334 s** avec l'interface existante. Pour 120 rigidités indépendantes,
il passe de **11,532 s à 0,758 s**, ou **0,264 s** avec les produits et
l'initialisation exacte du modèle vérifié. Le pic RSS passe alors de
**335 à 77 Mio**. Le premier import de SciPy ralentit toutefois le produit
isolé à 30 poutres de **3,63×** ; tous les coûts sont conservés.

Validation de livraison : **74 tests Rust**, **8 tests du prototype**,
**87 tests Python**, **41 groupes de vérification**, **46 bancs mécaniques**,
**9 bancs de contact**, contrôle des **187 entrées d'API**, sept tests du
juge d'adjoint et rejeu des archives. La roue finale est contrôlée hors
du dépôt. [Journaux et empreintes](bancs/version-0.9.0.json).

La [dérivation, les contrats et les mesures](OPERATEURS_ADJOINTS.md)
détaillent ce transfert de la différentiation implicite modulaire.
Les modes et les sensibilités statiques/modales restent denses ; le pont
garde toute la trajectoire et refactorise chaque pas.

Ces gains sont mesurés contre la roue **Vinkulum 0.8.2**, pas contre les
adjoints d'Exudyn, MBDyn ou Simpack. La
[confrontation externe à précision commune](CONFRONTATION_EXUDYN_0.8.2.md)
reste celle de la statique 0.8.2 ; aucune supériorité généraliste n'est
établie par cette livraison.
