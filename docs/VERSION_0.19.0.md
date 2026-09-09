# Vinkulum 0.19.0 — repères explicites et charges temporelles

Cette version étend le noyau pour l'éditeur rigide 3D Studio 0.2.0. L'API
historique `liaison` reste disponible. Les modèles existants ne sont pas migrés.

`liaison_reperes` reçoit les deux points et orientations locaux de la liaison.
Cela permet de conserver et diagnostiquer un montage incohérent sans reconstruire
implicitement le second repère à partir de la pose courante des corps.

`effort_temporel` applique trois lois de force et trois lois de moment mondiaux,
à un point local `p` du corps. Le moment au centre vaut
`M(t) + (R p) × F(t)`. Pour une perturbation de rotation spatiale gauche,
sa dérivée vaut `[F(t)]× [R p]×`. Le résidu dynamique et la tangente utilisent
le même instant de calcul ; la linéarisation locale inclut cette contribution.
Les lois constantes, linéaires et tabulées suivent les conventions natives.

Les analyses statiques par continuation ne prennent pas en charge cette nouvelle
famille de charges et la refusent explicitement. Les déclarations automatiques
de conservation d'énergie et de moment sont également refusées en présence de
ces efforts, y compris lorsqu'une loi particulière est constante ou nulle.

Huit contre-épreuves Python couvrent l'équivalence constante, les solutions
analytiques de translation et rotation, la force déportée avec référence DOP853,
la tangente numérique, les tables, les refus et la compatibilité des liaisons.
Une contre-épreuve Rust vérifie la dérivée du moment déporté et l'identité de
puissance. Voir la [référence d'API](API.md) et la [qualification Studio](STUDIO_3D.md).

Ces ajouts n'étendent pas les théorèmes Lean existants à tout problème comportant
ces charges. Ils n'apportent pas de certification générale du noyau ni de borne
d'erreur automatique sur les trajectoires.
