# Sondes préalables des séparateurs et de l'inertie binary64

Les refus sont conservés : Decimal 16 et binary64 en ordre naturel échouent
sur 512 poutres ; les séparateurs passent. La première sonde Decimal fait
18 appels, la sonde Python six, la sonde compilée six. B et les permutations
sont offerts à ces sondes de certificat ; elles ne classent pas les solveurs.
Leur provenance initiale est dans les chemins et manifestes capturés.

Le service pilote à 512 poutres recalcule toute sa préparation, avec le
validateur CSR accéléré ; un seul passage, distinct des 36 essais formels.
Ses sources ont été capturées après le calcul et leurs SHA comparés à ceux
que le worker a enregistrés avant et après l'exécution. Son module compilé
est celui de `natif-initial/build/`, identifié par son SHA dans le résultat.
Un champ NPY est exclu avec identité conservée ; ses octets restent dans
`/tmp/vinkulum-pilote-inertie-compilee-n512-v1`.

Les premières sondes ne capturaient pas tout l'environnement ni les sources
des dépendances. Leurs scripts historiques contiennent des chemins absolus.
Les archives et les résultats formels ne sont pas des reprises de ces sondes.
Aucun code ni binaire de cette archive n'est exécuté par les contrôles CI.
