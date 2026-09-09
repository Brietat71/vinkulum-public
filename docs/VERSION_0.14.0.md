# Vinkulum 0.14.0 — certification des contraintes exactement redondantes

La 0.14.0 étend `certifier_initialisation` avec l'option
`redondances=True`. Elle certifie le rang et la compatibilité des
contraintes stockées, l'accélération unique et la réaction généralisée
unique, même lorsque les multiplicateurs de liaison ne sont pas uniques.

Il s'agit d'une nouvelle capacité facultative, d'où le numéro mineur
**0.14.0**. Sans cette option, le contrat de la 0.13.0 est conservé.
Les lois mécaniques et les critères d'acceptation du solveur ordinaire
ne sont pas modifiés. Sans demande de certification, les contributions
temporaires de G sont libérées avant factorisation ; leur durée de vie
a été corrigée lors de la relecture de cette extension.

## Garantie et preuve

Le système complet est `M a + Gᵀ λ = f`, `G a = c`. Le certificat
vérifie une factorisation exacte `G=T C` et la compatibilité exacte
`c=T c_I`, où C contient les lignes retenues par le noyau.
Il établit aussi l'inversibilité du système avec jauge effectivement
résolu, par le certificat linéaire de la 0.13.0.

Ces propriétés impliquent `rang(G)=nombre de lignes de C`, l'existence
du système complet et l'unicité de a et de Gᵀλ. Les multiplicateurs
varient dans un espace affine de direction `ker(Gᵀ)`. Leur borne concerne
le représentant qui annule les multiplicateurs des lignes inactives ;
elle n'affirme ni une répartition unique ni une norme minimale.

Le [contrat et sa démonstration](CERTIFICATION_QUOTIENT.md) précisent
les bornes, les unités, les refus et la base de confiance. Le témoin
natif capture toutes les contributions de G avant filtrage, avec le
second membre complet, le masque et le vecteur réellement calculé.
Le générateur construit T par élimination rationnelle ; le vérificateur
indépendant contrôle les produits, les structures et les inégalités.

## Exemple et contre-exemple

Un pendule avec cinq contraintes dupliquées donne **dix lignes de rang
exact cinq**. Son [certificat archivé](bancs/certificat-quotient-pendule-0.14.0.json)
borne l'erreur du représentant numérique à environ **1,39 × 10⁻¹⁶**
sur le système exporté et fournit des bornes séparées de réaction généralisée.
Il est vérifiable sans installer Vinkulum :

```bash
python3 -I -S python/vinkulum/_verification_lineaire.py docs/bancs/certificat-quotient-pendule-0.14.0.json
```

Deux contraintes séparées par un bras de levier de 2⁻⁸⁰ peuvent être
classées redondantes par le seuil numérique du QR tout en restant
exactement indépendantes. Ce cas est reproduit sur le noyau ; la
certification refuse le quotient proposé. Ce refus porte sur la preuve
de rang exact et ne modifie pas la stratégie numérique ordinaire.

Deux lignes identiques avec seconds membres 0 et 2⁻¹⁰⁷⁴ sont également
refusées : la compatibilité n'est pas une comparaison à un seuil flottant.

## Validation et livraison

- 15 nouveaux tests de quotient : référence rationnelle indépendante,
  ambiguïté des multiplicateurs, dépendance 1/3, permutations, échelles,
  rang zéro, faux rang, incompatibilité et certificats altérés.
- Huit quotients de campagne vérifiés indépendamment : six pendules à
  trois échelles, un corps libre et une liaison non holonome dupliquée.
- 12 certificats linéaires et 13 314 contre-épreuves du garde entier
  conservés ; aucune divergence avec l'oracle exact.
- 175 tests Python sur la roue installée hors dépôt, 87 tests Rust,
  huit tests du prototype de poutre, 41 cas de vérification,
  46 bancs étendus et neuf cas de contact.
- 2 304 projections et six trajectoires de pendule : critères et résultats
  strictement identiques à ceux de la roue 0.13.0.
- Trois intégrations Exudyn restent ignorées faute de binaire installé.
  Le refus du label de couverture complète conserve le code 2 attendu.

Le [manifeste de livraison](bancs/version-0.14.0.json) et l'[archive des preuves](bancs/version-0.14.0-preuves.json.gz)
identifient les sources, la roue installée, les certificats et les journaux.
Roue : `vinkulum-0.14.0-cp314-cp314-linux_x86_64.whl`, SHA-256
`44fd977563ff413f1d3f975cbc78439e34abd3fb410d0200b849eae7080989b9`.
Environnement : CPython 3.14.7, NumPy 2.5.3, SciPy 1.18.1, Rust 1.98.0,
Linux x86_64. Les autres plateformes ne sont pas qualifiées par ce lot.

## Limites maintenues

**Le noyau entier n'est pas certifié.** La limite dense reste de 128
inconnues avant retrait des redondances, avec un budget par défaut de 64.
Contacts non lisses, partitions gelées et schémas de redémarrage imposés
restent exclus. Les redondances perturbées et les erreurs des coefficients
mécaniques ne sont pas couvertes ; l'erreur globale des trajectoires reste
ouverte. Aucune preuve formelle du code ni amélioration de vitesse n'est
revendiquée. Les temps de campagne sont diagnostiques.
