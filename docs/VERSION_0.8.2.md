# Vinkulum 0.8.2 — efforts analytiques et tangente exacte des poutres

7 septembre 2026 · Version corrective, API conservée.

Les moments nodaux sont dérivés par travail virtuel à partir des deux
énergies de poutre existantes. La tangente du champ d'efforts utilise un
passage de duaux d'ordre un. Elle remplace les différences finies des six
colonnes de rotation, qui exigeaient douze réévaluations des efforts.
Les sensibilités constitutives de raideur emploient ce même calcul exact,
en tenant compte des rigidités effectives de l'option intégrée.

La géométrie du modèle est conservée. Les séries près des rotations nulles
évitent les annulations ; les petites parties de la corde sont maintenues
dans les contributions au travail. La coupure du logarithme principal à π
reste explicitement non différentiable. Les critères physiques de Newton
et les restaurations après refus restent ceux de la 0.8.1.

La [dérivation et ses limites](TANGENTE_POUTRE_ANALYTIQUE.md) explicitent
les transports de rotation et le contre-calcul indépendant par énergie.
Trois nouveaux tests Rust couvrent efforts et tangentes, les repères,
les échelles et la coupure de branche. La
[lecture critique de la stack scientifique](COMPOSABILITE_SCIENTIFIQUE.md)
précise les hypothèses de composabilité des méthodes envisagées.

Les six rampes Princeton prennent **1,98 à 2,56 fois moins de temps**
que la roue 0.8.1 figée, avec 334 ou 335 évaluations principales.
L’appel complet d’analyse `k_c_m_z` ne gagne que 1,04–1,09 fois à
120 éléments : le coût global reste une cible de travail.

La campagne contre la roue 0.8.1 conserve 51 configurations et trois
répétitions par version, avec les contrôles physiques indépendants et les
restaurations. Les deux cas à petite échelle au seuil resserré restent
des échecs numériques ; une dérivée plus précise ne suffit pas à résoudre
leur limite de précision.

Validation : 73 tests Rust, 8 du prototype, 78 tests Python, 41 groupes
de vérification, 46 bancs mécaniques et 9 bancs de contact. La roue finale
est également contrôlée dans un environnement indépendant, hors du dépôt.
[Journaux et empreintes](bancs/version-0.8.2.json).

La confrontation externe reste celle de la
[roue 0.8.1](CONFRONTATION_STATIQUE_MBDYN_0.8.1.md). Les mesures internes
de cette version ne permettent pas de recalculer directement un rapport
de temps face à MBDyn. Aucune exécution Simpack n'est archivée.
