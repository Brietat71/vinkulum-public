# Correctif de robustesse statique — travail après 0.6.2

**Étape suivante :** [l’assemblage local des tangentes](RAIDEUR_LOCALE.md) réduit ensuite le coût de ces calculs. Les mesures ci-dessous sont celles du correctif initial.

Le chemin statique de Princeton converge désormais avec l'amortissement
actif, y compris les 50 paliers qui échouaient dans la confrontation.
Ce correctif ne démontre pas une supériorité générale sur MBDyn.

## Modification du noyau

Le premier essai conserve le critère en forces. Si ce palier échoue,
`statique()` restaure son état de départ, puis essaie un critère fondé sur
la correction des coordonnées physiques, calculée avec la factorisation
courante et l'équilibrage de Jacobi. Les réactions sont prises en compte
dans le résidu d'essai ; leur incrément n'entre pas dans la norme de pose.

Dans cette seconde stratégie, une absence de diminution des forces sur dix
itérations n'interrompt plus un chemin dont les corrections progressent.
Le budget d'itérations, les contrôles de finitude, la continuation en charge
et la restauration sur échec sont conservés. Les critères d'acceptation de
l'équilibre physique et leur plancher numérique ne sont pas relâchés.

Remplacer globalement le premier critère avait dégradé le cas de la
chaînette. Le choix de reprendre depuis l'état sauvegardé permet de conserver
ce cas et de résoudre la poutre anisotrope.

## Mesures isolées

Charge de 8,896 N à 45°, paramètres de la confrontation. Référence : position
du bout de la poutre MBDyn à 160 intervalles, archivée avec ses raffinements.
Trois répétitions après échauffement, un fil demandé. Temps de l'application
de charge et des appels statiques, **hors construction et démarrage Python**.
Ils ne doivent pas être présentés comme directement comparables aux temps
complets de l'exécutable MBDyn dans le premier rapport.

| Intervalles Vinkulum | Temps direct médian | Écart maximal au bout | Résidu maximal hors encastrement |
|---|---|---|---|
| 10 | 0.085 s | 254.87 µm | 3.5e-09 |
| 20 | 0.320 s | 63.72 µm | 7.3e-09 |
| 40 | 1.240 s | 15.93 µm | 1.6e-08 |
| 60 | 2.747 s | 7.08 µm | 2.1e-08 |

Le raffinement divise environ par quatre l'erreur lorsque le nombre
d'intervalles double. Le coût est encore proche d'une loi quadratique sur
cette plage. La précision spatiale et l'assemblage constituent donc des
limites restantes, même après résolution du défaut de robustesse.

Le chemin **50 paliers**, rejoué une fois à 10/20/40 intervalles, réussit
**3 fois sur 3** ; il retrouve les positions de l'appel direct à moins de
1e-8 m. Il ne s'agit pas de neuf répétitions ni d'une nouvelle médiane :
les neuf échecs initiaux appartiennent à l'archive de la version précédente.

## Vérification

- Le nouveau test anisotrope échoue avec le binaire antérieur et passe avec
  le correctif ; il vérifie les forces libres, les contraintes et le
  rapprochement vers MBDyn lorsque le maillage est raffiné.
- Formatage et Clippy sans avertissement ; **25 tests Rust**, **42 tests
  Python**, **41/41 groupes de vérification**, **46/46 bancs rapides** et
  **9/9 bancs de contact** passent. L'API générée reste à jour, 183 entrées.
- Le contrôle de continuation de la console a été adapté : il ne doit pas
  imposer artificiellement plusieurs paliers quand le solveur peut réussir
  en un seul. Un budget de huit itérations exerce encore la continuation,
  puis vérifie qu'elle retrouve le même équilibre.

[Mesures](bancs/statique-princeton-correctif.json) ·
[Journaux et empreintes](bancs/statique-princeton-controles.json) ·
[Objectif complet, toujours actif](OBJECTIF_MBDYN.md).

Reproduction des mesures :

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_statique.py --sortie /tmp/statique.json
```

La prochaine amélioration visée porte sur la raideur : assembler les
contributions locales au lieu de réévaluer toutes les poutres pour chaque
coordonnée du modèle. La confrontation doit ensuite être rejouée à
précision commune, sans assimiler la réparation de ce défaut à l'objectif
« considérablement mieux sur tous les tableaux ».
