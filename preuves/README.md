# Preuves du garde de vitesse — chantier de l'étape 3

Lean 4 **4.19.0**, Mathlib **v4.19.0** ; les commits des dépendances sont
figés dans `lake-manifest.json`. Depuis ce répertoire, avec le toolchain
indiqué dans `lean-toolchain` :

```sh
lake build
lake env lean Audit.lean
lake build garde
python contre_epreuves.py
```

La dernière commande requiert une installation compilée de Vinkulum et les
dépendances de sa campagne de certification. Elle compare les entiers complets
et la décision, sans conversion flottante des sorties, avec le Rust installé
et un oracle indépendant `fractions.Fraction`. Elle échoue à la première
divergence. Les mots d'entrée sont transmis comme entiers décimaux de 64 bits.

## Portée mathématique

`GardeVitesse.lean` prouve, sur des listes finies de dyadiques et des entiers
non bornés, l'équivalence exacte :

```
abs(R) * 2^46 ≤ T * 2^46 + E
    ↔
abs(dt + Σ g*u) ≤ tol + 2^-46 * (abs(dt) + Σ abs(g*u))
```

L'hypothèse est `tol ≥ 0`. L'unité entière est `2^-2148`. La preuve couvre
les produits signés, les sommes, les magnitudes, la boucle d'accumulation et
le décalage entier de la comparaison.

Les définitions exécutées sont regroupées dans `Modele.lean` et importées
par les preuves comme par `Main.lean`.

`Binary64.lean` interprète un mot fini par la significande IEEE 754 :
`f/2^52` avec exposant −1022 pour les sous-normaux, `1+f/2^52` avec exposant
`e−1023` pour les normaux. Il prouve l'égalité avec la valeur dyadique, les
bornes de mantisse et d'exposant et l'ajout du bit implicite par OR.

## Correspondance avec `src/certificat_vitesse.rs`

| Rust | Modèle Lean | Obligation |
|---|---|---|
| `f64::to_bits`, masques et décalages | `Mot64`, division et modulo | Masques et décalages naturels prouvés équivalents à division/modulo ; primitives Rust dans la base de confiance |
| `decode`, bit implicite | `decode`, `bit_implicite`, `decode_exact` | Valeur rationnelle prouvée pour tout mot décodé |
| `exposant` | `puissance − 1074` | Changement d'origine explicite |
| `u128` des mantisses | produit de naturels | Bornes des entrées finies prouvées dans `produit_u128` |
| `BigInt`, `BigUint` | entiers et naturels non bornés | Sémantique de num-bigint dans la base de confiance |
| boucle `ligne` | `accumule` | Égalité au résidu et à l'échelle prouvée |
| `accepte` | `garde_exact`, `decalage_exact` | Équivalence du critère prouvée |

## Base de confiance et limites

`Audit.lean` affiche les axiomes des théorèmes : `propext`, `Classical.choice`
et `Quot.sound` selon le théorème. Aucun `sorry`, `admit`, axiome ajouté ou
`native_decide` ne fait partie de ces preuves. Le noyau Lean, son exécution,
le matériel et la chaîne de distribution des outils restent de confiance.
Les objets `.olean` de Mathlib récupérés dans le cache sont aussi des
artefacts de confiance : cette qualification ne reconstruit pas toutes les
dépendances Lean depuis leurs sources. Le manifeste fige leurs révisions.
Les tactiques produisent des termes contrôlés par Lean ; leur succès seul
ne remplace pas ce contrôle.

L'exécutable de confrontation ajoute le compilateur Lean, son runtime et son
analyseur d'entrée. Côté Rust : compilateur, num-bigint, PyO3 et passage des
bits. La campagne apporte une vérification finie de correspondance ; elle
ne prouve pas l'équivalence de tous les programmes compilés.

Cette preuve ne porte ni sur les erreurs de géométrie, ni sur l'intégration
temporelle, ni sur l'exactitude physique des coefficients fournis au garde.
Elle ne constitue pas une certification du noyau entier. La revue de correspondance et la qualification de livraison restent requises.

## Exécution obligatoire en CI locale

`ci/local.sh` appelle `preuves/verifier.sh` après la campagne du garde Rust.
Une version différente de Lean, une preuve qui échoue, un axiome inattendu
ou une divergence interrompt la CI. Aucun téléchargement automatique ni
contournement si Lean manque. Ajouter le toolchain au `PATH` avant la CI.
Les dépendances et leurs objets peuvent être récupérés par
`lake exe cache get` lors de la première installation ; conserver le manifeste
versionné sans le régénérer. Le workflow `.github/workflows/ci.yml`
installe le toolchain vérifié par SHA-256, prépare Mathlib puis appelle
la même chaîne de vérification.

Archive Linux utilisée pour cette qualification :
`lean-4.19.0-linux.tar.zst`, distribution officielle
<https://github.com/leanprover/lean4/releases/tag/v4.19.0>, SHA-256 mesuré :
`6fe3ce97a58f44e2b3567d455b994eacec5bfe9ae7774f2a573444480ba813fe`.
Ce hash identifie l'archive téléchargée ; ce n'est pas une signature indépendante.

## Qualification conservée

`qualification.json` associe les empreintes des sources, la roue 0.14.1
isolée confrontée, les 13 audits de théorèmes et l'empreinte des entrées.
L'exécutable compilé et l'exécution interprétée ont chacun été confrontés
aux **13 314 cas**, sans divergence des entiers ni des décisions ; les
**7 entrées invalides** du protocole Lean sont refusées. Ces refus ne sont
pas un remplacement des sept contrôles invalides de l'interface Rust,
également exécutés par sa campagne existante.
