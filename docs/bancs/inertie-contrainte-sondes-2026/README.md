# Sondes préalables du certificat d'inertie contraint

Archive du 8 septembre 2026. Trois étapes exploratoires sont conservées
intégralement, avec leurs refus. **Cette archive ne constitue aucun
classement formel de solveurs.** Aucun solveur n'a été relancé pour
l'archivage.

- `sondes-double/` : motifs creux, signatures binary64 à 40/50/60/70/80/100 Hz
  et choix d'échelle sur les trois maillages natifs. Les directions viennent
  d'essais antérieurs et sont offertes à la sonde. Les signatures sont
  **non certifiées** ; le rapport conserve aussi les échelles défavorables.
- `pilote-initial/resultats.json.gz` : premier essai de certificats dirigés
  à γ = fl((2π·80)²), aux précisions Decimal 32 et 16. Les trois tailles
  passent à 32 chiffres ; 32 et 128 poutres passent à 16, tandis que
  **512 refuse à 16**, avec un pivot non séparé de zéro. Les directions
  préexistantes sont offertes. **Le script et la version initiale du
  prototype n'ont pas été figés.** Ce JSON conserve des preuves et diagnostics,
  mais le code actuel, qui a changé depuis, ne reconstitue pas à lui seul
  cette exécution historique.
- `pilote-controle-v1/` : pilote avec sélection fraîche, certificat dirigé
  à 32 chiffres, puis réponses contrôlées avec trois et quatre blocs de
  Krylov. Le script, ses sept sources, les B/Phi, les trois certificats,
  les résultats détaillés et le manifeste d'origine sont figés. Il s'agit
  d'**un seul passage** ; QR, sélection et certificat sont réutilisés entre
  les variantes à trois et quatre blocs. Les sommes de phases ne définissent
  donc pas un protocole comparatif de processus frais.

Au seuil demandé de 10⁻⁶, le pilote contrôlé présente les résultats suivants :

| Poutres | Trois blocs : majorants demandés | Quatre blocs : majorants demandés |
|---:|:---:|:---:|
| 32 | refus : masse et déformation dépassent le seuil | passent |
| 128 | passent | passent |
| 512 | refus : déformation dépasse le seuil | passent |

Ces refus de majorants ne sont pas des erreurs de résolution : le juge
physique indépendant passe les six champs calculés. Le certificat d'inertie
porte sur la coercivité stricte du modèle DᵀD−γM sur ker(Bᵀ). Il ne certifie
pas les arrondis des réponses ni ceux des majorants flottants du contrôleur.

## Provenance et intégrité

Les origines exactes sont `/tmp/vinkulum-sonde-inertie-2026`,
`/tmp/vinkulum-inertie-pilot-initial.json` et
`/tmp/vinkulum-inertie-controle-pilote-v1`. Le journal sibling
`/tmp/vinkulum-inertie-controle-pilote-v1.log` et le script parent
`/tmp/vinkulum-inertie-controle-pilote.py` sont également conservés dans
le dernier dossier. Ce script parent et `pilote.py` ont des octets identiques.
Les entrées natives et références physiques externes restent celles de
`/tmp/vinkulum-confrontation-ports-0.10.0-corrigee` ; aucun de ces modèles
n'est dupliqué ici. Les scripts conservent leurs chemins absolus historiques.

Tous les JSON d'origine sont compressés sans modification de leurs octets
décompressés. Les caches `.pyc`, déjà présents dans le manifeste pilote-v1,
sont eux aussi conservés en gzip sous `captures_cache/`, hors des dossiers
ignorés `__pycache__`. Ces huit captures binaires documentent l'état historique,
ne remplacent pas les sources `.py` et ne sont jamais exécutées par l'archive.
Leurs chemins d'origine et leurs empreintes décompressées restent dans la
provenance ; le déplacement ne retire aucune entrée du manifeste historique.
Le manifeste original est
`pilote-controle-v1/manifest.json.gz` ; ses 22 empreintes ont été vérifiées
avant et après copie, en décompressant les objets concernés.

Le `manifest.json` racine donne, pour chaque autre fichier archivé, sa taille
et son SHA256, puis le chemin, la taille et le SHA256 de ses octets d'origine
lorsqu'il s'agit d'une copie. Il indique la compression appliquée. Il ne
s'inclut pas lui-même. L'intégrité des copies et du manifeste historique a
été vérifiée sans exécuter les codes archivés. Adapter les chemins dans une
copie de travail serait nécessaire pour un rejeu ; ne pas écrire de nouvelles
mesures dans cette archive.
