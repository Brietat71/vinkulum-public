# Studio CAD — OCCT 8 minimum

Studio **0.4.0** ajoute une première chaîne de conception solide à l'éditeur
multicorps. Le noyau mécanique reste **Vinkulum 0.19.0**.

## Périmètre livré

Le bouton **CAD** et **Créer → Conception CAD…** ouvrent les opérations suivantes :
boîte, cylindre, sphère, extrusion de rectangle ou de disque XY, soustraction,
union, intersection, congés sur toutes les arêtes, import et export STEP.
Les modifications sont intégrées au document mécanique et à son annuler/rétablir.
Les opérations booléennes conservent le corps outil B : il reste un corps du
mécanisme tant que l'utilisateur ne le supprime pas.

Il s'agit de conception solide avec paramètres d'opération. Le journal conservé
est une provenance, pas encore un arbre de fonctions régénérable. L'esquisse
contrainte interactive, la sélection de faces/arêtes et les assemblages STEP
multi-pièces sont des étapes ultérieures. Le dialogue reste modal, avec le calcul
CAD exécuté dans un processus distinct, annulable et limité à 60 secondes.

## Installation reproductible

Dans l'environnement Studio Python 3.14 avec son noyau natif installé :

```sh
python ci/prepare_cad.py build/cad-sources
uv pip install build/cad-sources/build123d-0.11.1 \
  build/cad-sources/ocpsvg-0.6.0 './apps/studio[cad,test]'
python -m vinkulum_studio
```

Les distributions sources sont vérifiées par SHA-256 avant extraction. Les
patches sont conservés dans [ci/patches](../ci/patches). Aucun `--ignore-requires-python`,
aucune contrainte de dépendance ignorée et aucun repli vers OCCT 7 ne sont utilisés.
Les adaptations locales ont des versions distinctes :

| Composant | Version utilisée | Adaptation |
|---|---|---|
| OCCT via cadquery-ocp-novtk | 8.0.1 / 8.0.1.0.0 | Bibliothèque amont, sans recompilation ni modification |
| build123d | 0.11.1+vinkulum.occt8 | Imports des collections OCCT 8, limites de Bnd_Box, dépendances |
| ocpsvg | 0.6.0+vinkulum.occt8 | Collection de points et dépendance OCCT 8 |
| ocp_gordon | 0.3.1 | Version amont compatible OCCT 8 |

Les anciens alias TopTools/TColgp/TColStd sont remplacés par les types concrets
correspondants de `OCP.collections`. `Bnd_Box.Get()` renvoie désormais un type
non enregistré dans le binding ; les coins min/max fournissent les mêmes six
coordonnées. Les bibliothèques partagées restent séparées et remplaçables dans
le paquet. Les sources et licences amont restent identifiées dans
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

## Géométrie et mécanique

- Le BREP est conservé en millimètres. OCCT lit les unités déclarées dans STEP.
- Le repère du corps mécanique est recentré sur son centre de masse volumique.
  Les positions et le maillage de visualisation sont convertis en mètres.
- Avec une densité homogène ρ, `m = ρ V_mm³ × 10⁻⁹` et
  `I_kg·m² = (I_mm⁵ / V_mm³) × 10⁻⁶ × m`, autour du centre de masse.
  Les propriétés de masse viennent du BREP, indépendamment de la tessellation.
- Lors d'une modification CAD, les ancrages, repères de liaison et points de
  charge sont rebasés pour conserver leur position/orientation dans le monde.
- Le worker mécanique consomme le document capturé, la masse et le tenseur SI.
  Il ne charge ni OCCT ni build123d et ne calcule pas les inerties sur un maillage.
- Les projets CAD utilisent le schéma 2. Les projets antérieurs restent lisibles ;
  les projets sans CAD restent sauvegardés en schéma 1.

Un import accepte un seul solide valide, jusqu'à 8 Mo de STEP. Les données
capturées sont limitées à 2 Mo de BREP et 30 000 sommets/triangles par pièce,
4 Mo de BREP et 100 000 éléments de maillage par document. Un refus conserve
le document précédent. Ces bornes limitent cette première intégration ; elles
ne constituent pas une promesse de prise en charge de grands assemblages.

## Qualification

La qualification locale courante compte **48 tests Studio réussis**, dont les
parcours CAD et les contrôles de précision de l'inspecteur. Les journaux,
versions et mesures de la recette sont conservés dans
[le dossier Studio CAD 0.4.0](bancs/studio-cad-040/README.md).

La recette propre à Studio confronte volumes et inerties à des formules
analytiques indépendantes, vérifie rotations et centres de masse, booleans,
congés, conversions STEP mètres/millimètres, sauvegarde BREP, rebasage des
attaches, échec du worker et chaîne CAD → calcul mécanique → rendu Qt/VTK.

**96 tests amont ciblés** passent pour `test_bound_box`, `test_mass_properties`,
`test_location` et `test_build_part` de build123d 0.11.1. Cela qualifie les
parcours testés, pas l'intégralité de l'API build123d ou d'OCCT.
`uv pip check` confirme la cohérence du graphe installé. Le contrôle du paquet
exécute aussi son propre worker CAD, vérifie une pièce percée et son aller-retour
STEP avant de lancer le double pendule et le contrôle de rendu.

Les propriétés intégrales sont numériques et restent soumises aux tolérances
géométriques d'OCCT. La validation topologique et les références analytiques ne
sont pas une certification générale des trajectoires mécaniques.

Sources amont : [OCCT 8.0.1](https://github.com/Open-Cascade-SAS/OCCT/tree/V8_0_1),
[build123d](https://github.com/gumyr/build123d),
[OCP](https://github.com/CadQuery/OCP),
[ocpsvg](https://pypi.org/project/ocpsvg/0.6.0/),
[ocp_gordon](https://pypi.org/project/ocp-gordon/0.3.1/).
