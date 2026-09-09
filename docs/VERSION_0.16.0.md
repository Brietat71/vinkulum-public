# Vinkulum 0.16.0 — certificats creux et dépendances structurées

Cette version rend publiques les garanties du prototype :

- `certifier_systeme_creux` vérifie des facteurs creux et borne l'erreur pour
  toutes les perturbations déclarées, sans inverse dense.
- `certifier_quotient_structurel` vérifie un KKT perturbé et la conservation
  explicite du rang dans une famille `G'=T'C'`, avec bornes de réaction.

Le [contrat détaillé](CERTIFICATION_CREUSE.md) précise les hypothèses et
limites. La structure n'est pas déduite de contraintes proches. Les données
sont copiées ; le calcul dynamique du noyau n'est pas modifié dans ce lot.

Les formats `vinkulum.lineaire.creux.1` et `vinkulum.quotient.structurel.1`
se vérifient avec `verifier_certificat`, y compris son exécution autonome.
Les formats précédents et le prototype historique restent lisibles.
Cette extension compatible justifie l'augmentation mineure de version.

## Qualification isolée

Roue : `vinkulum-0.16.0-cp314-cp314-manylinux_2_39_x86_64.whl`.
SHA-256 : `d0779d3b2d5e0bcaf11843697e1997daccbaed779a5c9c884b16f14a64823374`.
Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, environnement isolé du checkout.

- 195 tests Python réussis depuis le paquet installé.
- Huit certificats de ce paquet vérifiés sans le module natif.
- Cinq archives creuses historiques relues par le nouveau vérificateur.
- 46/46 bancs rapides et 9/9 cas de contact réussis sur la roue isolée.
- Preuves et contre-épreuves historiques obligatoires dans la CI avant push.

Le manifeste `bancs/certification-creuse-0.16.0/qualification.json` conserve
les empreintes de la roue, des modules et des documents. La campagne physique
historique à 62 lignes n'est pas réattribuée à cette version ; ses résultats
restent ceux des roues identifiées dans le dossier 0.14.1.

## Travaux ouverts

Le chantier modal et les nouvelles confrontations MBDyn/Exudyn à précision
commune restent à qualifier. Le diagnostic de repères à temps long demeure
ouvert. Ces garanties algébriques ne certifient pas le noyau entier.
