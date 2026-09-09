# Sondes exploratoires du Krylov contraint

Cette archive conserve les fichiers disponibles de six séries exploratoires ayant précédé la campagne formelle de septembre 2026. Elle contient **43 observations** : 39 réponses et quatre essais de précision du certificat. Les résultats défavorables sont conservés : 12 réponses sont déclarées hors du seuil commun de 10⁻⁶. Aucun calcul physique n'a été relancé pour constituer cette archive.

**Les sources des implémentations ont changé entre les sondes et leurs anciennes versions n'ont pas toutes été figées. Aucune reproduction exacte ni aucun classement formel de performances ne sont revendiqués.** Les scripts présents sont des copies capturées après les expériences ; ils documentent les appels, pas l'identité complète du code exécuté à chaque étape. Aucun état actuel des dépendances n'est attribué rétroactivement aux sondes.

## Séries conservées

| Série | Observations | Résultat consigné | Provenance disponible |
|---|---:|---|---|
| Krylov initial | 9 | Deux blocs hors seuil aux trois tailles ; quatre et six blocs sous seuil | Script, journal et JSON intégral |
| Contrôle scalaire | 9 | Champs sous seuil à quatre/six blocs, mais majorants massiques trop larges ; certains indisponibles à deux blocs | Script et journal ; ancien JSON écrasé |
| Contrôle anisotrope | 9 | Majorants encore trop larges malgré des champs précis à quatre/six blocs | Journal et JSON intégral ; pas de script distinct conservé |
| Contrôle réparé | 9 | Deux blocs hors seuil ; quatre/six blocs avec champs et majorants relatifs déclarés sous seuil | Script, journal et JSON intégral |
| Sélection fraîche | 3 | Une direction recalculée et quatre blocs : trois champs déclarés sous seuil | Script, journal et JSON intégral |
| Précisions du complément | 4 | Certificats déclarés obtenus à 80/32/24/16 chiffres sur 512 poutres | Journal et JSON intégral ; script absent |

Les sondes de réponses portent sur 32/128/512 poutres, six charges et 257 fréquences jusqu'à 40 Hz. Les JSON conservés contiennent les juges complets ; leurs maxima et acceptations ont été vérifiés à partir des valeurs enregistrées, sans nouvelle résolution. Les 34 enregistrements JSON concordent avec les 34 lignes correspondantes des journaux. Les neuf lignes du contrôle scalaire ont seulement été transcrites dans la synthèse : elles ne restituent pas les tableaux perdus par fréquence et charge. Le fichier anisotrope ne remplace pas le JSON scalaire disparu.

## Enseignements et limites

Les essais initiaux orientent vers quatre blocs, avec des tailles réduites de 31/39/47 coordonnées. Deux blocs restent insuffisants au seuil retenu ; six blocs n'améliorent pas systématiquement les contrôles. Les étapes scalaire et anisotrope montrent qu'un champ précis peut coexister avec un majorant trop conservateur. L'étape réparée motive le travail sur les défauts de contrainte, sans prouver par ces fichiers seuls que tous les majorants absolus couvrent l'erreur réelle.

La sélection fraîche expose les coûts auparavant offerts : QR, direction modale et certificat sont recalculés. Le script mesure des phases séparées ; ce n'est pas le chronométrage intégral du worker formel. Son appel du certificat ne fixe pas explicitement la précision. Les essais de précision enregistrent notamment un minorant de 170485,66247046008 à 16 chiffres contre 170485,69949640392 à 80 chiffres. Ces observations ont orienté un réglage ensuite fixé avant la campagne ; elles ne remplacent jamais la décision d'un nouveau certificat.

Les temps bruts restent dans les rapports et `synthese.json`. Les premières sondes chargent gratuitement des directions et un certificat, ou excluent QR/Krylov de la préparation du contrôleur. Les champs complets et leurs normes ne sont pas chronométrés selon un contrat uniforme. Il n'y a donc ni comparaison à LU/HCB à coût égal, ni médiane de campagne, ni facteur de vitesse à tirer de ces fichiers.

Les majorants de réponse sont flottants. Les pilotes ne conservent pas les erreurs absolues et normes de référence nécessaires à une confrontation exhaustive de chaque majorant. « Sous seuil déclaré » décrit leurs nombres enregistrés ; cela ne désigne ni une certification machine des réponses ni l'acceptation plus exigeante de la campagne formelle.

## Organisation et intégrité

- `rapports/` : cinq JSON originaux compressés intégralement avec gzip déterministe ; décompression identique octet pour octet aux fichiers capturés.
- `journaux/` : six journaux complets, sans retrait des refus.
- `sources/` : quatre scripts disponibles, sans adaptation de leurs chemins. Les exécuter maintenant utiliserait d'autres états des dépendances.
- `synthese.json` : les 43 observations structurées, leur provenance et leurs limites.
- `manifest.json` : inventaire exact, tailles et SHA-256 ; SHA des originaux et transformations dans `provenance_fichiers`. Les absences sont explicites dans `fichiers_manquants`.

Le manifeste ne s'inclut pas lui-même. Ses empreintes garantissent l'intégrité de cette capture ; elles n'établissent pas la provenance d'exécution des anciennes implémentations. Cette archive complète les résultats de la campagne formelle sans en modifier les sources ni le classement.
