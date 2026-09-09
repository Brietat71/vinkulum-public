# Vinkulum 0.14.1 — fiabilité et qualification du noyau lisse

Cette livraison corrective exécute le premier lot du
[plan de fiabilité](PLAN_FIABILITE.md). Elle corrige un défaut du Newton
énergie–moment et deux défauts de suivi de la validation. Elle n'ajoute
pas de certificat géométrique ou temporel.

## Défaut mécanique reproduit et corrigé

Pour un solide libre, multiplier uniformément l'inertie J par un facteur
positif ne change pas l'équation d'Euler :

`J ω̇ = (J ω) × ω`.

La 0.14.0 arrêtait pourtant Newton lorsque le résidu du moment devenait
inférieur à `1e-17 * (1 + ||Π₀||)`. Le nombre 1 introduisait un plancher
dimensionné : une petite inertie pouvait faire accepter `Π₁ = Π₀` avant
la première correction, même hors des axes principaux. Le test de
stagnation avait aussi un plancher dimensionné.

Cas reproduit : J = diag(0,001 ; 0,002 ; 0,003), vitesse matérielle
(0,05 ; 20 ; 0,03) rad/s, durée 0,1 s et pas 0,001 s. Avec un facteur
d'inertie 1e-12, l'ancienne roue donne une erreur de vitesse matérielle
d'environ **0,0416005 rad/s** et une dérive relative du moment spatial
**0,00245323** (environ 0,245 %). La référence intègre séparément les
équations d'Euler ; elle n'utilise pas les matrices du noyau.

Le résidu est désormais `||f||∞ / max(||Π₀||∞, ||Π₁||∞)`, avec traitement
explicite du moment nul. Newton accepte au niveau epsilon ou, en cas de
stagnation, reprend son meilleur itéré uniquement si son résidu relatif
est inférieur à 32 epsilon. La norme infinie évite de mettre au carré
les petites ou grandes composantes du moment.

Après correction, le cas donne environ **2,64116e-7 rad/s**, correspondant
à l'erreur temporelle observée, avec une dérive du moment inférieure à
1e-14 dans la sonde. Les contre-épreuves couvrent huit facteurs de 1e-20
à 1e20, deux repères matériels et deux pas. Elles vérifient l'ordre deux,
les invariants et l'équivalence des mouvements. Le moment nul est testé
séparément. Aucun de ces contrôles flottants n'est une preuve globale.

## Validation : ne plus masquer une référence non exécutée

- Une exception d'Andrews était enregistrée sans marqueur d'erreur et
  pouvait laisser la commande sortir avec le code 0. Elle produit
  maintenant une erreur d'exécution identifiée et un code d'échec.
- Une exception d'une référence généraliste supprimait les références
  suivantes de la famille. Chacune est maintenant exécutée et enregistrée
  séparément ; son identité est conservée en cas d'échec.
- Les écarts physiques restent des résultats publiés. Leurs seuils ne
  changent pas, et un écart physique seul ne devient pas une exception.

Le pilote `ci/qualifie_validation.py` exécute les familles dans des
processus isolés, contrôle chaque fichier du paquet installé contre la
roue, enregistre les dépendances, les tolérances numériques et les logs,
puis rapproche chaque ligne de `docs/VALIDATION.md`. Ce rapport historique
n'est pas écrasé. Une famille interrompue ou une ligne historique absente
rend la qualification incomplète.

La reproduction complète sur les roues 0.14.0 et 0.14.1 retrouve les **62 lignes**
historiques : **44 dans leur tolérance, 18 écarts, aucune erreur d'exécution**.
Les valeurs et critères numériques sont identiques entre les deux roues.
Cette reproduction ne transforme pas ces 18 écarts en problèmes résolus.

## Corpus lisse et diagnostic des repères

Le corpus associe chute libre avec reprise, pendule et projections à
grand bras de levier, double pendule, quadrilatère articulé redondant,
plateforme de Stewart et toupie. Les critères sont identiques dans la
confrontation des deux roues : les **32 violations de critères** de la
0.14.0 dans cette campagne concernent la toupie ; la 0.14.1 est admise.
Ce nombre compte des critères, pas 32 mécanismes indépendants.

Le diagnostic à vingt secondes reste ouvert. Au pas 0,001 s, les deux
roues présentent le même maximum échantillonné : environ **1,72e-7** sur
la rotation en norme de Frobenius et **7,24e-7 rad/s** sur la vitesse
spatiale. En raffinant à 0,0005 s, ce dernier écart atteint **1,43e-6 rad/s**
sur la 0.14.1 contre **1,08e-6** sur la 0.14.0 ; à 0,00025 s il redescend.
La correction des petites inerties ne résout donc pas ce diagnostic.

Les équations variationnelles indépendantes montrent une forte sensibilité
sur ce cas, avec une norme de sensibilité de la vitesse matérielle proche
de **9,05e4 à dix secondes**. Cette mesure constitue un diagnostic ; elle
ne borne pas l'accumulation des arrondis et ne prouve pas à elle seule
l'origine de tous les écarts entre repères.

## Reproduire et vérifier

Installer les extras `validation` dans un environnement isolé, puis :

```bash
python ci/qualifie_validation.py --roue /chemin/roue.whl \
  --historique docs/VALIDATION.md --sortie /tmp/qualification-nouvelle --delai 3600
python ci/audit_fiabilite_lisse.py --sortie /tmp/lisse.json
python ci/audit_reperes_lisses.py --sortie /tmp/reperes.json
```

Les 181 tests Python de la roue isolée, 87 tests Rust, huit tests de poutre,
41 cas de vérification, 46 bancs étendus et neuf cas de contact passent.
Les trois intégrations nécessitant le binaire Exudyn restent ignorées.
Les 12 certificats linéaires, huit quotients et 13 314 contre-épreuves
exactes restent admis, sans divergence avec l'oracle.

La reconstruction des sources 0.14.0 retrouve les fichiers du paquet livré
à l'identique, extension native comprise. Les campagnes 0.14.1 portent sur
le candidat initial ; la roue finale actualise le README embarqué et le
tag de plateforme. **Tous les fichiers du paquet `vinkulum/` sont identiques
entre candidat et livraison**, ce que le manifeste vérifie octet par octet.
Les tests Python sont aussi exécutés sur la roue finale installée séparément.

Vérifier le dossier sans NumPy ni Vinkulum installé :

```bash
python3 -I -S ci/verifie_fiabilite.py docs/bancs/version-0.14.1.json
```

Les roues, dépendances, sources et résultats sont identifiés dans le
[manifeste de livraison](bancs/version-0.14.1.json). La qualification reste
limitée à Linux x86_64 et CPython 3.14. Les étapes Lean, assemblage par
intervalles, trajectoire validée et certificats creux restent à réaliser.
